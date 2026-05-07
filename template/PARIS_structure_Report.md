# {{ report_title }}

---

| Field | Value |
|-------|-------|
| **Pipeline Mode** | Structure |
| **Project ID** | {{ project_id }} |
| **Report Date** | {{ report_date }} |
| **Author** | {{ author }} |
| **Institution** | {{ institution }} |
| **PI / Investigator** | {{ pi_name }} |
| **Pipeline** | EasyOmics PARIS v1.0 |
| **Number of Samples** | {{ num_samples }} |
| **Samples** | {{ samples_list }} |
| **Reference Genome** | {{ reference_genome }} |

---

## Table of Contents

1. [Executive Summary](#executive-summary)
2. [Pipeline Parameters](#pipeline-parameters)
3. [Data Quality Assessment](#data-quality-assessment)
4. [Duplex Group Analysis](#duplex-group-analysis)
5. [Alternative Structures](#alternative-structures)
6. [Output Files](#output-files)
7. [Methods](#methods)

---

## Executive Summary

This report summarises the results of PARIS (Psoralen Analysis of RNA Interactions and
Structures) sequencing data processed through the **EasyOmics PARIS pipeline** in
**structure** mode.
The analysis identifies RNA duplex groups (DGs) from chimeric reads to characterise
intramolecular and intermolecular RNA secondary structures at transcriptome scale.

**Key highlights:**

- **{{ num_samples }}** sample(s) analysed
- Pipeline mode: **Structure** (genome-wide RNA secondary structure)
- Average mapping rate: **{{ avg_mapping_rate }}**
- Total duplex groups identified: **{{ total_dg_count }}**
- Total alternative structures: **{{ total_alt_structures }}**

---

## Pipeline Parameters

{{ pipeline_params_table }}

---

## Data Quality Assessment

### QC Overview

{{ qc_summary_table }}

> **Note:** Mapping rate ≥ 70%: PASS; 40–70%: WARN; < 40%: LOW.
> Chimeric read rate reflects the proportion of reads supporting RNA duplexes.

### Per-Sample QC Details

{{ per_sample_qc_section }}

{% if multiqc_report %}
### MultiQC Report

A comprehensive MultiQC report aggregating FastQC results is available at:

```
{{ multiqc_report }}
```

{% endif %}

---

## Duplex Group Analysis

Duplex groups (DGs) represent clusters of chimeric reads supporting a specific
RNA base-pairing interaction. Each DG defines a structural element in the
transcriptome revealed by PARIS crosslinking.

### Per-Sample Duplex Group Statistics

{{ dg_stats_section }}

---

## Alternative Structures

Alternative (competing) RNA conformations are identified by `alternativestructure.py`
from overlapping duplex groups, revealing structural heterogeneity at individual loci.

### Per-Sample Alternative Structure Summary

{{ alt_structure_section }}

---

## Output Files

| File | Description | Status |
|------|-------------|--------|
{{ output_files_table }}

---

## Methods

### Pipeline Overview

PARIS data were processed using the **EasyOmics PARIS Processing Pipeline**
(structure mode) implemented in Snakemake. Steps:

1. **Adapter Trimming** – 3′ adapter removal with Trimmomatic (minLen={{ min_len_3 }})
2. **Read Deduplication** – PCR duplicate collapse with `readCollapse`
3. **5′ Trimming** – HEADCROP (17 nt) + length filter (minLen={{ min_len_5 }}) with Trimmomatic
4. **Quality Control** – FastQC + MultiQC aggregation
5. **Genome Alignment** – STAR with chimeric read detection (chimSegmentMin={{ chim_segment_min }})
6. **Duplex Group Assembly** – `samPairingCalling` (minLen={{ dg_minlen }}, minPair={{ dg_minpair }})
7. **NG Assembly** – Near-gapless alignment extraction with `sam2ngmin.py`
8. **BED Conversion** – DGs to BED12 format with `dg2bed.py`
9. **Alternative Structure Analysis** – `alternativestructure.py`

### Software Versions

| Software | Version | Purpose |
|----------|---------|---------|
| Snakemake | — | Workflow management |
| FastQC | — | Read-level quality control |
| MultiQC | — | QC report aggregation |
| Trimmomatic | 0.39 | Adapter & quality trimming |
| STAR | — | Genome alignment with chimeric detection |
| SAMtools | — | BAM indexing & filtering |
| samPairingCalling | — | Duplex group assembly |
| sam2ngmin | — | NG alignment extraction |
| dg2bed | — | BED format conversion |
| alternativestructure | — | Alternative structure detection |

### References

1. Lu Z, et al. (2016). RNA Duplex Map in Living Cells Reveals Higher-Order
   Transcriptome Structure. *Cell* 165(5):1267–1279.

2. Dobin A, et al. (2013). STAR: ultrafast universal RNA-seq aligner.
   *Bioinformatics* 29(1):15–21.

3. Ewels P, et al. (2016). MultiQC: summarize analysis results for multiple
   tools and samples in a single report. *Bioinformatics* 32(19):3047–3048.

4. Bolger AM, et al. (2014). Trimmomatic: a flexible trimmer for Illumina
   sequence data. *Bioinformatics* 30(15):2114–2120.

---

## Contact Information

For questions regarding this report or the EasyOmics PARIS pipeline, please contact:

- **Institution:** {{ institution }}
- **PI / Investigator:** {{ pi_name }}
- **Author:** {{ author }}
- **Project ID:** {{ project_id }}

---

*Report generated by EasyOmics PARIS Pipeline (structure mode) on {{ report_date }}.*
