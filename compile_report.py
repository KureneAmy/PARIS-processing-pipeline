#!/usr/bin/env python3
"""
compile_report.py — PARIS Analysis Report Compiler
====================================================
Generates HTML and Markdown reports from PARIS pipeline outputs.

Supports both pipeline modes:
  - structure  : genome-wide RNA secondary structure via duplex groups (DGs)
  - interaction: intermolecular RNA–RNA interaction discovery

Usage:
    python compile_report.py --config config.yaml --work_dir /path/to/output
    python compile_report.py --config config.yaml --work_dir /path/to/output \\
        --html_out report_structure.html --md_out report_structure.md \\
        --inline_images
"""

import argparse
import base64
import logging
import os
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

import yaml

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Utility helpers
# ---------------------------------------------------------------------------

def _safe_float(value, default=None):
    """Convert value to float; return default on failure."""
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _safe_int(value, default=None):
    """Convert value to int; return default on failure."""
    try:
        return int(value)
    except (TypeError, ValueError):
        return default


def _fmt(value, decimals=1):
    """Format a numeric value with thousands separator; handle None / 'N/A'."""
    if value is None or value == "N/A":
        return "N/A"
    try:
        f = float(value)
        if decimals == 0:
            return f"{int(f):,}"
        return f"{f:,.{decimals}f}"
    except (TypeError, ValueError):
        return str(value)


def image_to_data_uri(filepath, mime_type="image/png"):
    """Encode an image file as a base64 data URI for inline HTML embedding."""
    if not filepath or not os.path.isfile(filepath):
        return ""
    try:
        with open(filepath, "rb") as fh:
            data = base64.b64encode(fh.read()).decode("ascii")
        return f"data:{mime_type};base64,{data}"
    except OSError as exc:
        logger.warning("Cannot encode image %s: %s", filepath, exc)
        return ""


def file_status_badge(filepath):
    """Return an HTML badge indicating whether a file exists."""
    if filepath and os.path.isfile(filepath):
        return '<span class="badge badge-pass">Found</span>'
    return '<span class="badge badge-low">Missing</span>'


# ---------------------------------------------------------------------------
# Log / file parsers
# ---------------------------------------------------------------------------

def parse_star_log(log_path):
    """
    Parse a STAR *_Log.final.out file and return alignment statistics.

    Returns a dict with keys:
        total_reads, uniquely_mapped, uniquely_mapped_pct,
        multi_mapped, multi_mapped_pct,
        chimeric_reads, chimeric_pct,
        mapping_rate (overall unique + multi)
    """
    stats = {}
    if not log_path or not os.path.isfile(log_path):
        logger.debug("STAR log not found: %s", log_path)
        return stats

    try:
        with open(log_path, encoding="utf-8", errors="replace") as fh:
            content = fh.read()
    except OSError as exc:
        logger.warning("Cannot read STAR log %s: %s", log_path, exc)
        return stats

    def _extract(pattern, text, cast=_safe_float):
        m = re.search(pattern, text, re.IGNORECASE | re.MULTILINE)
        if m:
            return cast(m.group(1))
        return None

    stats["total_reads"] = _extract(
        r"Number of input reads\s+\|\s+([0-9]+)", content, _safe_int
    )
    stats["uniquely_mapped"] = _extract(
        r"Uniquely mapped reads number\s+\|\s+([0-9]+)", content, _safe_int
    )
    stats["uniquely_mapped_pct"] = _extract(
        r"Uniquely mapped reads %\s+\|\s+([0-9.]+)%?", content
    )
    stats["multi_mapped"] = _extract(
        r"Number of reads mapped to multiple loci\s+\|\s+([0-9]+)", content, _safe_int
    )
    stats["multi_mapped_pct"] = _extract(
        r"% of reads mapped to multiple loci\s+\|\s+([0-9.]+)%?", content
    )
    stats["chimeric_reads"] = _extract(
        r"Number of chimeric reads\s+\|\s+([0-9]+)", content, _safe_int
    )
    stats["chimeric_pct"] = _extract(
        r"% of chimeric reads\s+\|\s+([0-9.]+)%?", content
    )

    # Derive overall mapping rate (unique + multi)
    u = stats.get("uniquely_mapped_pct")
    m_pct = stats.get("multi_mapped_pct")
    if u is not None:
        stats["mapping_rate"] = round(u + (m_pct or 0.0), 2)
    else:
        stats["mapping_rate"] = None

    return stats


def parse_trimmomatic_log(log_path):
    """
    Parse a Trimmomatic log for read counts.

    Returns a dict with keys: input_reads, surviving, dropped, surviving_pct
    """
    stats = {}
    if not log_path or not os.path.isfile(log_path):
        logger.debug("Trimmomatic log not found: %s", log_path)
        return stats

    try:
        with open(log_path, encoding="utf-8", errors="replace") as fh:
            content = fh.read()
    except OSError:
        return stats

    # "Input Reads: 1000000 Surviving: 900000 (90.00%) Dropped: 100000 (10.00%)"
    m = re.search(
        r"Input Reads:\s+([0-9]+)\s+Surviving:\s+([0-9]+)\s+\(([0-9.]+)%\)\s+Dropped:\s+([0-9]+)",
        content,
        re.IGNORECASE,
    )
    if m:
        stats["input_reads"] = _safe_int(m.group(1))
        stats["surviving"] = _safe_int(m.group(2))
        stats["surviving_pct"] = _safe_float(m.group(3))
        stats["dropped"] = _safe_int(m.group(4))
    return stats


def parse_dedup_log(log_path):
    """
    Parse a readCollapse deduplication log.

    Returns a dict with: input_reads, output_reads, dup_rate
    """
    stats = {}
    if not log_path or not os.path.isfile(log_path):
        logger.debug("Dedup log not found: %s", log_path)
        return stats

    try:
        with open(log_path, encoding="utf-8", errors="replace") as fh:
            content = fh.read()
    except OSError:
        return stats

    # Look for lines like "Input: N" / "Output: N"
    m_in = re.search(r"(?:Input|Total reads?)[:\s]+([0-9,]+)", content, re.IGNORECASE)
    m_out = re.search(r"(?:Output|Unique reads?)[:\s]+([0-9,]+)", content, re.IGNORECASE)

    if m_in:
        stats["input_reads"] = _safe_int(m_in.group(1).replace(",", ""))
    if m_out:
        stats["output_reads"] = _safe_int(m_out.group(1).replace(",", ""))

    if stats.get("input_reads") and stats.get("output_reads"):
        dup_count = stats["input_reads"] - stats["output_reads"]
        stats["dup_rate"] = round(100.0 * dup_count / stats["input_reads"], 2) if stats["input_reads"] else 0.0

    return stats


def count_fastq_reads(fastq_path):
    """Count reads in a FASTQ file (handles .gz). Returns None if not found."""
    if not fastq_path or not os.path.isfile(fastq_path):
        return None
    try:
        import gzip
        opener = gzip.open if fastq_path.endswith(".gz") else open
        with opener(fastq_path, "rt", encoding="utf-8", errors="replace") as fh:
            count = sum(1 for _ in fh)
        return count // 4
    except (OSError, EOFError):
        return None


def parse_dg_geometric(filepath):
    """
    Parse a *_DG.geometric file and return duplex group statistics.

    The geometric file has one DG per block separated by blank lines or
    lines starting with 'Group'. Each Group header contains metadata;
    alignment lines follow.
    """
    stats = {
        "total_dg": 0,
        "total_reads": 0,
        "max_support": 0,
        "mean_support": 0.0,
        "median_support": 0.0,
        "unique_rnas": 0,
    }
    if not filepath or not os.path.isfile(filepath):
        logger.debug("DG geometric file not found: %s", filepath)
        return stats

    group_supports = []
    rna_set = set()

    try:
        with open(filepath, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                if line.startswith("Group"):
                    # "Group 1 readnum: 5 ..."
                    m = re.search(r"readnum[:\s]+([0-9]+)", line, re.IGNORECASE)
                    if m:
                        group_supports.append(_safe_int(m.group(1), default=0))
                    # Extract RNA identifiers from the group line
                    rna_matches = re.findall(
                        r"([A-Za-z0-9_.\-]+)\([^)]*\):[0-9]+-[0-9]+", line
                    )
                    for rna in rna_matches:
                        rna_set.add(rna)
    except OSError as exc:
        logger.warning("Cannot read DG geometric file %s: %s", filepath, exc)
        return stats

    if group_supports:
        import statistics as _stat
        stats["total_dg"] = len(group_supports)
        stats["total_reads"] = sum(group_supports)
        stats["max_support"] = max(group_supports)
        stats["mean_support"] = round(_stat.mean(group_supports), 1)
        try:
            stats["median_support"] = _stat.median(group_supports)
        except Exception:
            stats["median_support"] = 0

    stats["unique_rnas"] = len(rna_set)
    return stats


def parse_interactions_filtered(filepath, top_n=10):
    """
    Parse an *_interactions.filtered file and return interaction statistics.

    Returns a dict with:
        total_pairs, unique_rna_pairs, top_pairs (list of (rna1, rna2, count)),
        rna1_counts, rna2_counts
    """
    stats = {
        "total_pairs": 0,
        "unique_rna_pairs": 0,
        "top_pairs": [],
        "top_rna1": "—",
        "top_rna2": "—",
    }
    if not filepath or not os.path.isfile(filepath):
        logger.debug("Interactions filtered file not found: %s", filepath)
        return stats

    pair_counts = Counter()
    rna1_counts = Counter()
    rna2_counts = Counter()

    try:
        with open(filepath, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                # Group line format: "Group N readnum: K ... RNA1(...):start-end|RNA2(...):start-end"
                if line.startswith("Group"):
                    rna_pairs = re.findall(
                        r"([A-Za-z0-9_.\-]+)\([^)]*\):[0-9]+-[0-9]+\|([A-Za-z0-9_.\-]+)\([^)]*\):[0-9]+-[0-9]+",
                        line,
                    )
                    for rna1, rna2 in rna_pairs:
                        key = tuple(sorted([rna1, rna2]))
                        pair_counts[key] += 1
                        rna1_counts[rna1] += 1
                        rna2_counts[rna2] += 1
    except OSError as exc:
        logger.warning("Cannot read interactions file %s: %s", filepath, exc)
        return stats

    stats["total_pairs"] = sum(pair_counts.values())
    stats["unique_rna_pairs"] = len(pair_counts)
    stats["top_pairs"] = [
        (k[0], k[1], v)
        for k, v in pair_counts.most_common(top_n)
    ]

    if rna1_counts:
        stats["top_rna1"] = rna1_counts.most_common(1)[0][0]
    if rna2_counts:
        stats["top_rna2"] = rna2_counts.most_common(1)[0][0]

    return stats


def parse_alt_file(filepath):
    """
    Parse a *_DG.alt alternative structure file.

    Returns a dict with: total_alt_structures, total_loci
    """
    stats = {"total_alt_structures": 0, "total_loci": 0}
    if not filepath or not os.path.isfile(filepath):
        logger.debug("Alt file not found: %s", filepath)
        return stats

    loci = set()
    alt_count = 0

    try:
        with open(filepath, encoding="utf-8", errors="replace") as fh:
            for line in fh:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                alt_count += 1
                # First field is typically a locus identifier
                parts = line.split("\t")
                if parts:
                    loci.add(parts[0])
    except OSError as exc:
        logger.warning("Cannot read alt file %s: %s", filepath, exc)
        return stats

    stats["total_alt_structures"] = alt_count
    stats["total_loci"] = len(loci)
    return stats


# ---------------------------------------------------------------------------
# Data collection
# ---------------------------------------------------------------------------

def collect_stats(config, work_dir, inline_images=True):
    """
    Walk pipeline output directories and collect all statistics.

    Returns a report_data dict ready for template rendering.
    """
    samples = config.get("samples", [])
    mode = config.get("mode", "structure")
    report_cfg = config.get("report", {})
    params = config.get("params", {})

    # Reference genome path — derive a display name
    ref_fa = config.get("paths", {}).get("genome_fa", "")
    reference_genome = os.path.basename(ref_fa) if ref_fa else "—"

    report_data = {
        "report_title": report_cfg.get("title", f"PARIS Analysis Report — {mode.capitalize()} Mode"),
        "author": report_cfg.get("author", ""),
        "institution": report_cfg.get("institution", ""),
        "pi_name": report_cfg.get("pi_name", ""),
        "project_id": report_cfg.get("project_id", ""),
        "report_date": datetime.now().strftime("%Y-%m-%d"),
        "mode": mode,
        "num_samples": len(samples),
        "samples_list": ", ".join(samples) if samples else "N/A",
        "reference_genome": reference_genome,
        "multiqc_report": os.path.join(work_dir, "qc", "multiqc_report.html"),
        # Pipeline params for display
        "min_len_3": params.get("minlen3", 28),
        "min_len_5": params.get("minlen5", 20),
        "chim_segment_min": params.get("chimSegmentMin", 15),
        "dg_minlen": params.get("dg_minlen", 15),
        "dg_minpair": params.get("dg_minpair", 2),
        # Aggregated KPIs (filled below)
        "avg_mapping_rate": "N/A",
        "total_dg_count": "N/A",
        "total_alt_structures": "N/A",
        "total_interaction_pairs": "N/A",
        "total_unique_rna_pairs": "N/A",
        "vis_enabled": config.get("rna_intrxn_visualization", {}).get("enabled", False),
        "per_sample": {},
    }

    mapping_rates = []
    total_dg = 0
    total_alt = 0
    total_intrxn = 0
    total_unique_pairs = 0

    for sample in samples:
        sample_dir = os.path.join(work_dir, sample)
        log_dir = os.path.join(work_dir, "logs")
        sd = {}

        # ------------------------------------------------------------------
        # STAR alignment stats
        # ------------------------------------------------------------------
        # STAR writes *_Log.final.out in the prefix directory
        if mode == "structure":
            star_log = os.path.join(
                sample_dir,
                f"{sample}_starGenome_Log.final.out",
            )
        else:
            star_log = os.path.join(
                sample_dir,
                f"{sample}_starSmallRNA_Log.final.out",
            )
        star_stats = parse_star_log(star_log)

        # ------------------------------------------------------------------
        # Trimmomatic / dedup stats
        # ------------------------------------------------------------------
        trim3_log = os.path.join(log_dir, f"{sample}_trim3.log")
        trim5_log = os.path.join(log_dir, f"{sample}_trim5.log")
        dedup_log = os.path.join(log_dir, f"{sample}_dedup.log")

        trim3_stats = parse_trimmomatic_log(trim3_log)
        trim5_stats = parse_trimmomatic_log(trim5_log)
        dedup_stats = parse_dedup_log(dedup_log)

        # Derive total reads from STAR or trim log
        raw_reads = (
            star_stats.get("total_reads")
            or trim3_stats.get("input_reads")
        )
        clean_reads = trim5_stats.get("surviving") or trim3_stats.get("surviving")

        # Duplication rate from dedup
        dup_rate = dedup_stats.get("dup_rate")

        # Mapping stats
        uniquely_mapped = star_stats.get("uniquely_mapped")
        mapping_rate = star_stats.get("mapping_rate")
        chimeric_reads = star_stats.get("chimeric_reads")
        chimeric_pct = star_stats.get("chimeric_pct")

        if mapping_rate is not None:
            mapping_rates.append(mapping_rate)
        sd["mapping_rate_pct"] = mapping_rate
        sd["qc_status"] = True

        # Determine status badge string
        if mapping_rate is not None:
            if mapping_rate >= 70:
                status_badge = '<span class="badge badge-pass">PASS</span>'
                status_md = "PASS"
            elif mapping_rate >= 40:
                status_badge = '<span class="badge badge-warn">WARN</span>'
                status_md = "WARN"
            else:
                status_badge = '<span class="badge badge-low">LOW</span>'
                status_md = "LOW"
        else:
            status_badge = '<span class="badge badge-info">N/A</span>'
            status_md = "N/A"

        # Build per-sample QC rows (booktabs table rows)
        qc_items = [
            ("Total Reads (raw)",        _fmt(raw_reads, 0)),
            ("Clean Reads (post-trim)",  _fmt(clean_reads, 0)),
            ("Duplication Rate",         f"{_fmt(dup_rate, 2)}%" if dup_rate is not None else "N/A"),
            ("Uniquely Mapped Reads",    _fmt(uniquely_mapped, 0)),
            ("Unique Mapping Rate",      f"{_fmt(star_stats.get('uniquely_mapped_pct'), 2)}%" if star_stats.get("uniquely_mapped_pct") is not None else "N/A"),
            ("Overall Mapping Rate",     f"{_fmt(mapping_rate, 2)}%" if mapping_rate is not None else "N/A"),
            ("Chimeric Reads",           _fmt(chimeric_reads, 0)),
            ("Chimeric Rate",            f"{_fmt(chimeric_pct, 2)}%" if chimeric_pct is not None else "N/A"),
        ]
        sd["qc_rows"] = "\n".join(
            f"<tr><td>{k}</td><td>{v}</td></tr>" for k, v in qc_items
        ) or "<tr><td colspan='2' class='text-muted text-center'>QC data not available</td></tr>"

        # QC summary row for the aggregated booktabs table
        sd["qc_summary_row"] = (
            f"<tr>"
            f"<td>{sample}</td>"
            f"<td>{_fmt(raw_reads, 0)}</td>"
            f"<td>{_fmt(clean_reads, 0)}</td>"
            f"<td>{_fmt(uniquely_mapped, 0)}</td>"
            f"<td>{_fmt(mapping_rate, 2)}%</td>"
            f"<td>{_fmt(chimeric_reads, 0)}</td>"
            f"<td>{_fmt(dup_rate, 2)}%</td>"
            f"<td>{status_badge}</td>"
            f"</tr>"
        )
        sd["qc_summary_md"] = (
            f"| {sample} "
            f"| {_fmt(raw_reads, 0)} "
            f"| {_fmt(clean_reads, 0)} "
            f"| {_fmt(uniquely_mapped, 0)} "
            f"| {_fmt(mapping_rate, 2)}% "
            f"| {_fmt(chimeric_reads, 0)} "
            f"| {_fmt(dup_rate, 2)}% "
            f"| {status_md} |"
        )

        # ------------------------------------------------------------------
        # Mode-specific results
        # ------------------------------------------------------------------
        if mode == "structure":
            dg_file = os.path.join(sample_dir, f"{sample}_DG.geometric")
            dg_stats = parse_dg_geometric(dg_file)

            alt_file = os.path.join(sample_dir, f"{sample}.alt")
            alt_stats = parse_alt_file(alt_file)

            sd["dg_total"] = _fmt(dg_stats["total_dg"], 0)
            sd["dg_unique_rnas"] = _fmt(dg_stats["unique_rnas"], 0)
            sd["dg_median_support"] = _fmt(dg_stats["median_support"], 1)
            sd["ng_count"] = "See NGmin SAM"

            total_dg += dg_stats["total_dg"]
            total_alt += alt_stats["total_alt_structures"]

            dg_stat_items = [
                ("Total Duplex Groups",       _fmt(dg_stats["total_dg"], 0)),
                ("Total Supporting Reads",    _fmt(dg_stats["total_reads"], 0)),
                ("Max Read Support",          _fmt(dg_stats["max_support"], 0)),
                ("Mean Read Support",         _fmt(dg_stats["mean_support"], 1)),
                ("Median Read Support",       _fmt(dg_stats["median_support"], 1)),
                ("Unique RNA Loci",           _fmt(dg_stats["unique_rnas"], 0)),
            ]
            if dg_stats["total_dg"] > 0:
                sd["dg_stats_rows"] = "\n".join(
                    f"<tr><td>{k}</td><td>{v}</td></tr>"
                    for k, v in dg_stat_items
                )
            else:
                sd["dg_stats_rows"] = ""

            alt_stat_items = [
                ("Total Alternative Structures", _fmt(alt_stats["total_alt_structures"], 0)),
                ("Unique Loci",                  _fmt(alt_stats["total_loci"], 0)),
            ]
            if alt_stats["total_alt_structures"] > 0:
                sd["alt_stats_rows"] = "\n".join(
                    f"<tr><td>{k}</td><td>{v}</td></tr>"
                    for k, v in alt_stat_items
                )
            else:
                sd["alt_stats_rows"] = ""

            # DG BED figure (optional)
            sd["dg_bed_figure"] = ""  # Would be set if a BED coverage PNG exists

        else:  # interaction mode
            intrxn_file = os.path.join(sample_dir, f"{sample}_interactions.filtered")
            intrxn_stats = parse_interactions_filtered(intrxn_file)

            sd["interaction_pairs"] = _fmt(intrxn_stats["total_pairs"], 0)
            sd["unique_rna_pairs"] = _fmt(intrxn_stats["unique_rna_pairs"], 0)
            sd["top_rna1"] = intrxn_stats.get("top_rna1", "—")
            sd["top_rna2"] = intrxn_stats.get("top_rna2", "—")

            total_intrxn += intrxn_stats["total_pairs"]
            total_unique_pairs += intrxn_stats["unique_rna_pairs"]

            intrxn_stat_items = [
                ("Total Interaction Pairs",    _fmt(intrxn_stats["total_pairs"], 0)),
                ("Unique RNA Interaction Pairs", _fmt(intrxn_stats["unique_rna_pairs"], 0)),
                ("Most Frequent RNA (arm 1)",   intrxn_stats.get("top_rna1", "—")),
                ("Most Frequent RNA (arm 2)",   intrxn_stats.get("top_rna2", "—")),
            ]
            if intrxn_stats["total_pairs"] > 0:
                sd["interaction_stats_rows"] = "\n".join(
                    f"<tr><td>{k}</td><td>{v}</td></tr>"
                    for k, v in intrxn_stat_items
                )
            else:
                sd["interaction_stats_rows"] = ""

            # Top pairs table rows
            if intrxn_stats["top_pairs"]:
                sd["top_pairs_rows"] = "\n".join(
                    f"<tr><td>{r1}</td><td>{r2}</td><td>{_fmt(cnt, 0)}</td><td>{i+1}</td></tr>"
                    for i, (r1, r2, cnt) in enumerate(intrxn_stats["top_pairs"])
                )
            else:
                sd["top_pairs_rows"] = ""

            # Arc diagram figures
            plots_dir = os.path.join(sample_dir, "plots")
            arc_figures = []
            if os.path.isdir(plots_dir):
                for fname in sorted(os.listdir(plots_dir)):
                    if fname.lower().endswith((".png", ".svg", ".jpg")):
                        fpath = os.path.join(plots_dir, fname)
                        if inline_images:
                            mime = "image/svg+xml" if fname.endswith(".svg") else "image/png"
                            src = image_to_data_uri(fpath, mime)
                        else:
                            src = os.path.relpath(fpath, start=work_dir)
                        if src:
                            arc_figures.append({
                                "src": src,
                                "caption": fname.replace("_", " ").rsplit(".", 1)[0],
                            })
            sd["arc_figures"] = arc_figures
            sd["network_figure"] = ""  # Future: network graph image

        report_data["per_sample"][sample] = sd

    # ------------------------------------------------------------------
    # Aggregated KPIs
    # ------------------------------------------------------------------
    if mapping_rates:
        avg_mr = sum(mapping_rates) / len(mapping_rates)
        report_data["avg_mapping_rate"] = f"{avg_mr:.1f}%"
    else:
        report_data["avg_mapping_rate"] = "N/A"

    if mode == "structure":
        report_data["total_dg_count"] = _fmt(total_dg, 0)
        report_data["total_alt_structures"] = _fmt(total_alt, 0)
    else:
        report_data["total_interaction_pairs"] = _fmt(total_intrxn, 0)
        report_data["total_unique_rna_pairs"] = _fmt(total_unique_pairs, 0)

    return report_data


# ---------------------------------------------------------------------------
# HTML / MD context builders
# ---------------------------------------------------------------------------

def _pipeline_params_rows(config):
    """Build HTML table rows for pipeline parameters."""
    p = config.get("params", {})
    paths = config.get("paths", {})
    rows = [
        ("Mode",                   config.get("mode", "—")),
        ("Threads",                p.get("threads", "—")),
        ("Min Length after 3′ Trim", p.get("minlen3", "—")),
        ("Min Length after 5′ Trim", p.get("minlen5", "—")),
        ("HEADCROP",               17),
        ("chimSegmentMin",         p.get("chimSegmentMin", "—")),
        ("chimJunctionOverhangMin", p.get("chimJunctionOverhangMin", "—")),
        ("DG minLen",              p.get("dg_minlen", "—")),
        ("DG minPair",             p.get("dg_minpair", "—")),
        ("dg2bed option",          p.get("dg2bed_option", "—")),
        ("Reference FASTA",        os.path.basename(paths.get("genome_fa", "—"))),
        ("GTF Annotation",         os.path.basename(paths.get("genome_gtf", "—"))),
    ]
    return "\n".join(
        f"<tr><td>{k}</td><td><code>{v}</code></td></tr>" for k, v in rows
    )


def _pipeline_params_table_md(config):
    """Build Markdown table for pipeline parameters."""
    p = config.get("params", {})
    paths = config.get("paths", {})
    rows = [
        ("Mode",                   config.get("mode", "—")),
        ("Threads",                p.get("threads", "—")),
        ("Min Length after 3′ Trim", p.get("minlen3", "—")),
        ("Min Length after 5′ Trim", p.get("minlen5", "—")),
        ("HEADCROP",               17),
        ("chimSegmentMin",         p.get("chimSegmentMin", "—")),
        ("chimJunctionOverhangMin", p.get("chimJunctionOverhangMin", "—")),
        ("DG minLen",              p.get("dg_minlen", "—")),
        ("DG minPair",             p.get("dg_minpair", "—")),
        ("dg2bed option",          p.get("dg2bed_option", "—")),
        ("Reference FASTA",        os.path.basename(paths.get("genome_fa", "—"))),
        ("GTF Annotation",         os.path.basename(paths.get("genome_gtf", "—"))),
    ]
    header = "| Parameter | Value |\n|-----------|-------|\n"
    body = "\n".join(f"| {k} | `{v}` |" for k, v in rows)
    return header + body


def _output_files_rows(config, work_dir, samples, mode):
    """Build HTML table rows for output files listing."""
    rows = []
    for sample in samples:
        sample_dir = os.path.join(work_dir, sample)
        files_to_check = [
            (os.path.join(sample_dir, f"{sample}_trim3.fastq"),
             f"{sample}_trim3.fastq", "3′-trimmed reads"),
            (os.path.join(sample_dir, f"{sample}_trim3_nodup.fastq"),
             f"{sample}_trim3_nodup.fastq", "Deduplicated reads"),
            (os.path.join(sample_dir, f"{sample}_trim5.fastq"),
             f"{sample}_trim5.fastq", "5′-trimmed reads"),
        ]
        if mode == "structure":
            files_to_check += [
                (os.path.join(sample_dir, f"{sample}_starGenome_Aligned.sortedByCoord.out.bam"),
                 f"{sample}_starGenome.bam", "STAR genome-aligned BAM"),
                (os.path.join(sample_dir, f"{sample}_starGenome_Chimeric.out.junction"),
                 f"{sample}_starGenome_Chimeric.out.junction", "Chimeric junctions"),
                (os.path.join(sample_dir, f"{sample}_DG.geometric"),
                 f"{sample}_DG.geometric", "Duplex groups (geometric format)"),
                (os.path.join(sample_dir, f"{sample}_DG.geometricsam"),
                 f"{sample}_DG.geometricsam", "Duplex groups SAM"),
                (os.path.join(sample_dir, f"{sample}_NGmin.sam"),
                 f"{sample}_NGmin.sam", "NG-assembly SAM"),
                (os.path.join(sample_dir, f"{sample}_DG.bed"),
                 f"{sample}_DG.bed", "Duplex groups (BED12)"),
                (os.path.join(sample_dir, f"{sample}.alt"),
                 f"{sample}.alt", "Alternative structure calls"),
            ]
        else:
            files_to_check += [
                (os.path.join(sample_dir, f"{sample}_starSmallRNA_Chimeric.out.junction"),
                 f"{sample}_starSmallRNA_Chimeric.out.junction", "Small RNA chimeric junctions"),
                (os.path.join(sample_dir, f"{sample}_interactions.geometric"),
                 f"{sample}_interactions.geometric", "Interaction duplex groups"),
                (os.path.join(sample_dir, f"{sample}_interactions.filtered"),
                 f"{sample}_interactions.filtered", "Filtered inter-molecular interactions"),
            ]

        for fpath, fname, desc in files_to_check:
            badge = file_status_badge(fpath)
            rows.append(f"<tr><td><code>{fname}</code></td><td>{desc}</td><td>{badge}</td></tr>")

    # QC files
    rows.append(
        f"<tr><td><code>qc/multiqc_report.html</code></td><td>MultiQC aggregated QC report</td>"
        f"<td>{file_status_badge(os.path.join(work_dir, 'qc', 'multiqc_report.html'))}</td></tr>"
    )
    return "\n".join(rows) or "<tr><td colspan='3'>No files found</td></tr>"


def _output_files_table_md(config, work_dir, samples, mode):
    """Build Markdown table rows for output files listing."""
    rows = []
    for sample in samples:
        sample_dir = os.path.join(work_dir, sample)
        files_to_check = []
        if mode == "structure":
            files_to_check = [
                (f"{sample}_DG.geometric",   "Duplex groups (geometric format)"),
                (f"{sample}_DG.geometricsam","Duplex groups SAM"),
                (f"{sample}_NGmin.sam",       "NG-assembly SAM"),
                (f"{sample}_DG.bed",          "Duplex groups (BED12)"),
                (f"{sample}.alt",             "Alternative structure calls"),
            ]
        else:
            files_to_check = [
                (f"{sample}_interactions.geometric", "Interaction duplex groups"),
                (f"{sample}_interactions.filtered",  "Filtered inter-molecular interactions"),
            ]
        for fname, desc in files_to_check:
            exists = os.path.isfile(os.path.join(sample_dir, fname))
            status = "✓ Found" if exists else "✗ Missing"
            rows.append(f"| `{fname}` | {desc} | {status} |")

    exists = os.path.isfile(os.path.join(work_dir, "qc", "multiqc_report.html"))
    rows.append(f"| `qc/multiqc_report.html` | MultiQC aggregated QC report | {'✓ Found' if exists else '✗ Missing'} |")
    return "\n".join(rows)


def build_html_context(report_data, config, work_dir):
    """Build the full context dict for HTML template rendering."""
    ctx = dict(report_data)
    samples = config.get("samples", [])
    mode = report_data["mode"]

    ctx["pipeline_params_rows"] = _pipeline_params_rows(config)
    ctx["output_files_rows"] = _output_files_rows(config, work_dir, samples, mode)

    # QC summary table rows
    ctx["qc_summary_rows"] = "\n".join(
        sd.get("qc_summary_row", "")
        for sd in report_data["per_sample"].values()
    ) or "<tr><td colspan='8' class='text-muted text-center'>No QC data available</td></tr>"

    return ctx


def build_md_context(report_data, config, work_dir):
    """Build the full context dict for Markdown template rendering."""
    ctx = dict(report_data)
    samples = config.get("samples", [])
    mode = report_data["mode"]

    ctx["pipeline_params_table"] = _pipeline_params_table_md(config)
    ctx["output_files_table"] = _output_files_table_md(config, work_dir, samples, mode)

    # QC summary table
    header = (
        "| Sample | Total Reads | Clean Reads | Uniquely Mapped | Mapping Rate | "
        "Chimeric Reads | Dup. Rate | Status |\n"
        "|--------|-------------|-------------|-----------------|--------------|"
        "---------------|-----------|--------|\n"
    )
    ctx["qc_summary_table"] = header + "\n".join(
        sd.get("qc_summary_md", "") for sd in report_data["per_sample"].values()
    )

    # Per-sample QC details section
    qc_sections = []
    for sname, sd in report_data["per_sample"].items():
        qc_items_md = re.findall(r"<td>([^<]+)</td>\s*<td>([^<]+)</td>", sd.get("qc_rows", ""))
        if qc_items_md:
            tbl = "| Metric | Value |\n|--------|-------|\n"
            tbl += "\n".join(f"| {k} | {v} |" for k, v in qc_items_md)
        else:
            tbl = "_QC data not available._"
        qc_sections.append(f"#### {sname}\n\n{tbl}\n")
    ctx["per_sample_qc_section"] = "\n".join(qc_sections) or "_No QC data available._"

    if mode == "structure":
        # DG stats section
        dg_sections = []
        for sname, sd in report_data["per_sample"].items():
            dg_items = re.findall(r"<td>([^<]+)</td>\s*<td>([^<]+)</td>", sd.get("dg_stats_rows", ""))
            if dg_items:
                tbl = "| Metric | Value |\n|--------|-------|\n"
                tbl += "\n".join(f"| {k} | {v} |" for k, v in dg_items)
            else:
                tbl = f"_Duplex group data not available for {sname}._"
            dg_sections.append(f"#### {sname}\n\n{tbl}\n")
        ctx["dg_stats_section"] = "\n".join(dg_sections) or "_No duplex group data available._"

        # Alt structure section
        alt_sections = []
        for sname, sd in report_data["per_sample"].items():
            alt_items = re.findall(r"<td>([^<]+)</td>\s*<td>([^<]+)</td>", sd.get("alt_stats_rows", ""))
            if alt_items:
                tbl = "| Metric | Value |\n|--------|-------|\n"
                tbl += "\n".join(f"| {k} | {v} |" for k, v in alt_items)
            else:
                tbl = f"_Alternative structure data not available for {sname}._"
            alt_sections.append(f"#### {sname}\n\n{tbl}\n")
        ctx["alt_structure_section"] = "\n".join(alt_sections) or "_No alternative structure data available._"

    else:  # interaction
        # Interaction stats section
        intrxn_sections = []
        for sname, sd in report_data["per_sample"].items():
            items = re.findall(r"<td>([^<]+)</td>\s*<td>([^<]+)</td>", sd.get("interaction_stats_rows", ""))
            if items:
                tbl = "| Metric | Value |\n|--------|-------|\n"
                tbl += "\n".join(f"| {k} | {v} |" for k, v in items)
            else:
                tbl = f"_Interaction data not available for {sname}._"
            intrxn_sections.append(f"#### {sname}\n\n{tbl}\n")
        ctx["interaction_stats_section"] = "\n".join(intrxn_sections) or "_No interaction data available._"

        # Top pairs section
        top_sections = []
        for sname, sd in report_data["per_sample"].items():
            pairs = re.findall(
                r"<td>([^<]+)</td>\s*<td>([^<]+)</td>\s*<td>([^<]+)</td>\s*<td>([^<]+)</td>",
                sd.get("top_pairs_rows", ""),
            )
            if pairs:
                tbl = "| RNA 1 | RNA 2 | Read Support | Rank |\n|-------|-------|--------------|------|\n"
                tbl += "\n".join(f"| {r1} | {r2} | {cnt} | {rank} |" for r1, r2, cnt, rank in pairs)
            else:
                tbl = f"_No interaction pairs found for {sname}._"
            top_sections.append(f"#### {sname}\n\n{tbl}\n")
        ctx["top_pairs_section"] = "\n".join(top_sections) or "_No interaction pair data available._"

    return ctx


# ---------------------------------------------------------------------------
# Template rendering
# ---------------------------------------------------------------------------

def _render_template(template_path, context):
    """
    Render a Jinja2 template; fall back to simple {{key}} substitution
    if Jinja2 is not available.
    """
    try:
        from jinja2 import Environment, FileSystemLoader, Undefined

        class SilentUndefined(Undefined):
            """Return empty string for missing variables instead of raising."""
            def __str__(self):
                return ""
            def __iter__(self):
                return iter([])
            def __bool__(self):
                return False

        env = Environment(
            loader=FileSystemLoader(str(Path(template_path).parent)),
            autoescape=False,
            keep_trailing_newline=True,
            undefined=SilentUndefined,
        )
        template = env.get_template(Path(template_path).name)
        return template.render(**context)
    except ImportError:
        logger.warning("Jinja2 not available; using basic {{key}} substitution.")
        with open(template_path, encoding="utf-8") as fh:
            content = fh.read()
        for key, value in context.items():
            if isinstance(value, (str, int, float)):
                content = content.replace("{{" + key + "}}", str(value))
                content = content.replace("{{ " + key + " }}", str(value))
        return content


def generate_reports(config, work_dir, html_template, md_template,
                     html_out, md_out, inline_images=True):
    """
    Collect statistics, build contexts, render templates, and write output files.
    """
    report_data = collect_stats(config, work_dir, inline_images=inline_images)

    html_ctx = build_html_context(report_data, config, work_dir)
    md_ctx = build_md_context(report_data, config, work_dir)

    # Write HTML
    if html_template and os.path.isfile(html_template):
        html_content = _render_template(html_template, html_ctx)
        Path(html_out).parent.mkdir(parents=True, exist_ok=True)
        with open(html_out, "w", encoding="utf-8") as fh:
            fh.write(html_content)
        logger.info("HTML report written to %s", html_out)
    else:
        logger.warning("HTML template not found: %s", html_template)

    # Write Markdown
    if md_template and os.path.isfile(md_template):
        md_content = _render_template(md_template, md_ctx)
        Path(md_out).parent.mkdir(parents=True, exist_ok=True)
        with open(md_out, "w", encoding="utf-8") as fh:
            fh.write(md_content)
        logger.info("Markdown report written to %s", md_out)
    else:
        logger.warning("Markdown template not found: %s", md_template)

    logger.info("Report generation complete (mode=%s).", report_data["mode"])
    return html_out, md_out


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Compile PARIS analysis reports (HTML and Markdown)."
    )
    parser.add_argument(
        "--config", default="config.yaml",
        help="Path to pipeline config.yaml (default: config.yaml).",
    )
    parser.add_argument(
        "--work_dir", default=None,
        help="Path to pipeline work/output directory. Overrides config paths.work_dir.",
    )
    parser.add_argument(
        "--html_template", default=None,
        help="Path to HTML Jinja2 template. Auto-detected from mode if not given.",
    )
    parser.add_argument(
        "--md_template", default=None,
        help="Path to Markdown template. Auto-detected from mode if not given.",
    )
    parser.add_argument(
        "--html_out", default=None,
        help="Output HTML file path. Defaults to <work_dir>/report_<mode>.html.",
    )
    parser.add_argument(
        "--md_out", default=None,
        help="Output Markdown file path. Defaults to <work_dir>/report_<mode>.md.",
    )
    parser.add_argument(
        "--inline_images", action="store_true", default=False,
        help="Embed images as base64 data URIs in the HTML for offline viewing.",
    )
    parser.add_argument(
        "--mode", default=None,
        help="Override pipeline mode (structure or interaction).",
    )
    return parser.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    # Load config
    config_path = args.config
    if not os.path.isfile(config_path):
        logger.error("Config file not found: %s", config_path)
        sys.exit(1)

    with open(config_path, encoding="utf-8") as fh:
        config = yaml.safe_load(fh)

    # Override mode if requested
    if args.mode:
        config["mode"] = args.mode

    # Override work_dir if requested
    work_dir = args.work_dir or config.get("paths", {}).get("work_dir", ".")

    mode = config.get("mode", "structure")
    report_cfg = config.get("report", {})

    # Resolve template paths — check config first, then script-relative templates/
    script_dir = Path(__file__).parent.resolve()
    templates_dir = script_dir / "templates"

    def _resolve_template(arg_val, cfg_key, default_name):
        if arg_val:
            return arg_val
        cfg_val = report_cfg.get(cfg_key)
        if cfg_val and os.path.isfile(cfg_val):
            return cfg_val
        candidate = templates_dir / default_name
        if candidate.is_file():
            return str(candidate)
        logger.warning("Template not found: %s", default_name)
        return None

    html_template = _resolve_template(
        args.html_template, "html_template",
        f"PARIS_{mode}_Report.html",
    )
    md_template = _resolve_template(
        args.md_template, "md_template",
        f"PARIS_{mode}_Report.md",
    )

    # Determine inline_images: CLI flag or config key
    inline_images = args.inline_images or report_cfg.get("inline_images", False)

    # Output paths
    html_out = (
        args.html_out
        or report_cfg.get("output_html")
        or os.path.join(work_dir, f"report_{mode}.html")
    )
    md_out = (
        args.md_out
        or report_cfg.get("output_md")
        or os.path.join(work_dir, f"report_{mode}.md")
    )

    generate_reports(
        config=config,
        work_dir=work_dir,
        html_template=html_template,
        md_template=md_template,
        html_out=html_out,
        md_out=md_out,
        inline_images=inline_images,
    )


if __name__ == "__main__":
    main()
