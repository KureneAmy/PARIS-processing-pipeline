# {{ report_title }}

---

| Field | Value |
|-------|-------|
| **Pipeline Mode** | Interaction |
| **Project ID** | {{ project_id }} |
| **Report Date** | {{ report_date }} |
| **Author** | {{ author }} |
| **Institution** | {{ institution }} |
| **PI / Investigator** | {{ pi_name }} |
| **Pipeline** | EasyOmics PARIS v1.0 |
| **Number of Samples** | {{ num_samples }} |
| **Samples** | {{ samples_list }} |
| **Reference (Small RNA)** | {{ reference_genome }} |

---

## Table of Contents

1. [Executive Summary](#executive-summary)
2. [Pipeline Parameters](#pipeline-parameters)
3. [Data Quality Assessment](#data-quality-assessment)
4. [Interaction Analysis](#interaction-analysis)
5. [Top RNA Interaction Pairs](#top-rna-interaction-pairs)
6. [Arc Diagram Visualization](#arc-diagram-visualization)
7. [Output Files](#output-files)
8. [Methods](#methods)

---

## Executive Summary

This report summarises the results of PARIS (Psoralen Analysis of RNA Interactions and
Structures) sequencing data processed through the **EasyOmics PARIS pipeline** in
**interaction** mode.
The analysis identifies intermolecular RNA–RNA interactions from chimeric reads by
mapping to a small RNA reference, enabling discovery of RNA interaction networks.

**Key highlights:**

- **{{ num_samples }}** sample(s) analysed
- Pipeline mode: **Interaction** (RNA–RNA interaction discovery)
- Average mapping rate: **{{ avg_mapping_rate }}**
- Total interaction pairs identified: **{{ total_interaction_pairs }}**
- Unique RNA interaction pairs: **{{ total_unique_rna_pairs }}**

---

## Pipeline Parameters

{{ pipeline_params_table }}

---

## Data Quality Assessment

### QC Overview

{{ qc_summary_table }}

> **Note:** Mapping rate ≥ 70%: PASS; 40–70%: WARN; < 40%: LOW.
> In interaction mode, unique mapping (outFilterMultimapNmax=1) is used to
> ensure high-confidence chimeric read assignment.

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

## Interaction Analysis

In interaction mode, PARIS captures intermolecular RNA–RNA contacts through psoralen
crosslinking. Only trans-interactions (different RNA species on each arm) are retained
after filtering.

### Per-Sample Interaction Statistics

{{ interaction_stats_section }}

---

## Top RNA Interaction Pairs

{{ top_pairs_section }}

> **Note:** Interactions are ranked by the number of supporting chimeric read pairs.
> Only the top pairs per sample are shown here; the full list is in the output files.

---

## Arc Diagram Visualization

{% if vis_enabled %}
Arc diagrams for configured RNA pairs are available in the output directory:

```
{work_dir}/{sample}/plots/
```

Each arc diagram shows the distribution of interaction sites along RNA molecules.
Configure RNA pairs in `config.yaml` under `rna_intrxn_visualization.rna_pairs`.
{% else %}
Arc diagram visualization was not enabled in this run.
To enable, set `rna_intrxn_visualization.enabled: true` and configure `rna_pairs`
in `config.yaml`.
{% endif %}

---

## Output Files

| File | Description | Status |
|------|-------------|--------|
{{ output_files_table }}

---

## Methods

### Pipeline Overview

PARIS data were processed using the **EasyOmics PARIS Processing Pipeline**
(interaction mode) implemented in Snakemake. Steps:

1. **Adapter Trimming** – 3′ adapter removal with Trimmomatic (minLen={{ min_len_3 }})
2. **Read Deduplication** – PCR duplicate collapse with `readCollapse`
3. **5′ Trimming** – HEADCROP (17 nt) + length filter (minLen={{ min_len_5 }}) with Trimmomatic
4. **Quality Control** – FastQC + MultiQC aggregation
5. **Small RNA Alignment** – STAR with unique mapping (outFilterMultimapNmax=1) and chimeric read detection (chimSegmentMin={{ chim_segment_min }})
6. **Interaction Assembly** – `samPairingCalling` (minLen={{ dg_minlen }}, minPair={{ dg_minpair }})
7. **Interaction Filtering** – Only inter-molecular events retained (different RNA on each chimeric arm)
{% if vis_enabled %}
8. **Interaction Visualization** – Arc diagrams for configured RNA pairs with `intrxn_specificity_edited.py`
{% endif %}

### Software Versions

| Software | Version | Purpose |
|----------|---------|---------|
| Snakemake | — | Workflow management |
| FastQC | — | Read-level quality control |
| MultiQC | — | QC report aggregation |
| Trimmomatic | 0.39 | Adapter & quality trimming |
| STAR | — | Small RNA alignment with chimeric detection |
| SAMtools | — | BAM processing |
| samPairingCalling | — | Interaction duplex group assembly |
| intrxn_specificity | — | Interaction arc diagram visualization |

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

*Report generated by EasyOmics PARIS Pipeline (interaction mode) on {{ report_date }}.*
