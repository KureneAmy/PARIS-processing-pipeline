import os

configfile: "config.yaml"

SAMPLES = config["samples"]
WORK_DIR = config["paths"]["work_dir"]
MODE = config["mode"]

# 设置输出目录结构
LOG_DIR = f"{WORK_DIR}/logs"

def raw_fastq(wc):
    return config["raw_files"][wc.sample]

############################################
# rule all (direct dependencies)
############################################
if MODE == "structure":
    rule all:
        input:
            f"{WORK_DIR}/.setup_done",
            expand(f"{WORK_DIR}/{{sample}}/{{sample}}_starGenome_Aligned.sortedByCoord.out.bam", sample=SAMPLES),
            expand(f"{WORK_DIR}/{{sample}}/{{sample}}_starGenome_Aligned.sortedByCoord.out.bam.bai", sample=SAMPLES),
            expand(f"{WORK_DIR}/{{sample}}/{{sample}}_NGmin.sam", sample=SAMPLES),
            expand(f"{WORK_DIR}/{{sample}}/{{sample}}_DG.bed", sample=SAMPLES),
            expand(f"{WORK_DIR}/{{sample}}/{{sample}}.alt", sample=SAMPLES) if config["alternative_structure"]["enabled"] else [],
            f"{WORK_DIR}/qc/multiqc_report.html",
else:  # interaction
    rule all:
        input:
            expand(f"{WORK_DIR}/{{sample}}/{{sample}}_interactions.filtered", sample=SAMPLES),
            expand(f"{WORK_DIR}/{{sample}}/{{sample}}_intrxn_visual.done", sample=SAMPLES) if config["rna_intrxn_visualization"]["enabled"] else [],
            f"{WORK_DIR}/qc/multiqc_report.html",

############################################
# 0) Setup directories
############################################
rule setup_directories:
    output:
        touch(f"{WORK_DIR}/.setup_done")
    shell:
        r"""
        mkdir -p {WORK_DIR}/logs
        mkdir -p {WORK_DIR}/qc
        """

############################################
# 1) Preprocessing
############################################
rule fastqc_raw:
    input:
        raw_fastq
    output:
        f"{WORK_DIR}/qc/{{sample}}/{{sample}}_fastqc.html"
    log:
        f"{LOG_DIR}/{{sample}}_fastqc_raw.log"
    container:
        config["container"]
    shell:
        r"""
        mkdir -p {WORK_DIR}/qc/{wildcards.sample}/
        fastqc -o {WORK_DIR}/qc/{wildcards.sample}/ {input} 2>{log}
        """

rule trim3:
    input:
        raw_fastq
    output:
        f"{WORK_DIR}/{{sample}}/{{sample}}_trim3.fastq"
    log:
        f"{LOG_DIR}/{{sample}}_trim3.log"
    params:
        adapters=config["paths"]["adapters_3p"],
        threads=config["params"]["threads"],
        minlen=config["params"]["minlen3"]
    container:
        config["container"]
    shell:
        r"""
        mkdir -p {WORK_DIR}/{wildcards.sample}
        java -jar /mnt/software/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads {params.threads} -phred33 \
          {input} {output} \
          ILLUMINACLIP:{params.adapters}:3:20:10 \
          SLIDINGWINDOW:10:30 MINLEN:{params.minlen} 2>{log}
        """

rule dedup:
    input:
        f"{WORK_DIR}/{{sample}}/{{sample}}_trim3.fastq"
    output:
        f"{WORK_DIR}/{{sample}}/{{sample}}_trim3_nodup.fastq"
    log:
        f"{LOG_DIR}/{{sample}}_dedup.log"
    params:
         script=config["paths"]["read_collapse"]
    container:
        config["container"]
    shell:
        r"""
        {params.script} {input} {output} 2>{log}
        """

rule trim5:
    input:
        f"{WORK_DIR}/{{sample}}/{{sample}}_trim3_nodup.fastq"
    output:
        f"{WORK_DIR}/{{sample}}/{{sample}}_trim5.fastq"
    log:
        f"{LOG_DIR}/{{sample}}_trim5.log"
    params:
        threads=config["params"]["threads"],
        minlen=config["params"]["minlen5"]
    container:
        config["container"]
    shell:
        r"""
        java -jar /mnt/software/Trimmomatic-0.39/trimmomatic-0.39.jar SE -threads {params.threads} -phred33 \
          {input} {output} \
          HEADCROP:17 MINLEN:{params.minlen} 2>{log}
        """

rule fastqc:
    input:
        f"{WORK_DIR}/{{sample}}/{{sample}}_trim5.fastq"
    output:
        f"{WORK_DIR}/qc/{{sample}}/{{sample}}_trim5_fastqc.html"
    log:
        f"{LOG_DIR}/{{sample}}_fastqc_trim5.log"
    container:
        config["container"]
    shell:
        r"""
        fastqc -o {WORK_DIR}/qc/{wildcards.sample}/ {input} 2>{log}
        """

############################################
# 2) STAR mapping
############################################

def star_common_params():
    p = config["params"]
    return {
        "chimSegmentMin": p["chimSegmentMin"],
        "chimJunctionOverhangMin": p["chimJunctionOverhangMin"],
        "threads": p["star_threads"],
    }

if MODE == "structure":
    rule star_map_genome_structure:
        input:
            f"{WORK_DIR}/{{sample}}/{{sample}}_trim5.fastq"
        output:
            aligned=f"{WORK_DIR}/{{sample}}/{{sample}}_starGenome_Aligned.sortedByCoord.out.bam",
            chim=f"{WORK_DIR}/{{sample}}/{{sample}}_starGenome_Chimeric.out.sam",
            junc=f"{WORK_DIR}/{{sample}}/{{sample}}_starGenome_Chimeric.out.junction"
        log:
            f"{LOG_DIR}/{{sample}}_star_genome.log"
        params:
            index=config["paths"]["star_genome_index"],
            multimap=config["params"]["outFilterMultimapNmax_structure"],
            **star_common_params()
        container:
            config["container"]
        shell:
            r"""
            STAR --runMode alignReads --genomeDir {params.index} \
            --readFilesIn {input} \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix {WORK_DIR}/{wildcards.sample}/{wildcards.sample}_starGenome_ \
            --outReadsUnmapped Fastq \
            --outSAMattributes All \
            --outFilterMultimapNmax {params.multimap} \
            --alignIntronMin 1 \
            --scoreGapNoncan -4 --scoreGapATAC -4 \
            --chimSegmentMin {params.chimSegmentMin} \
            --chimOutType SeparateSAMold Junctions \
            --chimJunctionOverhangMin {params.chimJunctionOverhangMin} \
            --runThreadN {params.threads} 2>{log}
            """

    rule starGenome_bam_index:
        input:
            bam=f"{WORK_DIR}/{{sample}}/{{sample}}_starGenome_Aligned.sortedByCoord.out.bam"
        output:
            bai=f"{WORK_DIR}/{{sample}}/{{sample}}_starGenome_Aligned.sortedByCoord.out.bam.bai"
        log:
            f"{LOG_DIR}/{{sample}}_bam_index.log"
        container:
            config["container"]
        shell:
            r"""
            samtools index {input.bam} 2>{log}
            """

else:  # interaction mode
    rule star_map_smallrna_interaction:
        input:
            f"{WORK_DIR}/{{sample}}/{{sample}}_trim5.fastq"
        output:
            chim=f"{WORK_DIR}/{{sample}}/{{sample}}_starSmallRNA_Chimeric.out.sam",
            junc=f"{WORK_DIR}/{{sample}}/{{sample}}_starSmallRNA_Chimeric.out.junction"
        log:
            f"{LOG_DIR}/{{sample}}_star_smallrna.log"
        params:
            index=config["paths"]["star_smallrna_index"],
            multimap=config["params"]["outFilterMultimapNmax_interaction"],
            **star_common_params()
        container:
            config["container"]
        shell:
            r"""
            STAR --runMode alignReads --genomeDir {params.index} \
            --readFilesIn {input} \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix {WORK_DIR}/{wildcards.sample}/{wildcards.sample}_starSmallRNA_ \
            --outReadsUnmapped Fastq \
            --outSAMattributes All \
            --outFilterMultimapNmax {params.multimap} \
            --alignIntronMin 1 \
            --scoreGapNoncan -4 --scoreGapATAC -4 \
            --chimSegmentMin {params.chimSegmentMin} \
            --chimJunctionOverhangMin {params.chimJunctionOverhangMin} \
            --chimOutType SeparateSAMold Junctions \
            --runThreadN {params.threads} 2>{log}
            """

############################################
# 3) Structure branch
############################################
if MODE == "structure":
    rule primary_gapped_sam:
        input:
            f"{WORK_DIR}/{{sample}}/{{sample}}_starGenome_Aligned.sortedByCoord.out.bam"
        output:
            f"{WORK_DIR}/{{sample}}/{{sample}}_Aligned_prim_N.sam"
        log:
            f"{LOG_DIR}/{{sample}}_primary_gapped.log"
        container:
            config["container"]
        shell:
            r"""
            samtools view -h -F 0x900 {input} | awk '($1 ~ /^@/) || ($6 ~ /N/)' > {output} 2>{log}
            """

    rule dg_assembly_structure:
        input:
            aligned=f"{WORK_DIR}/{{sample}}/{{sample}}_Aligned_prim_N.sam",
            chim=f"{WORK_DIR}/{{sample}}/{{sample}}_starGenome_Chimeric.out.sam",
            junc=f"{WORK_DIR}/{{sample}}/{{sample}}_starGenome_Chimeric.out.junction"
        output:
            dg=f"{WORK_DIR}/{{sample}}/{{sample}}_DG.geometric",
            dgsam=f"{WORK_DIR}/{{sample}}/{{sample}}_DG.geometricsam"
        log:
            f"{LOG_DIR}/{{sample}}_dg_assembly.log"
        params:
            script=config["paths"]["samPairingCalling"],
            module_dir=config["paths"]["PARIS_module"],
            bin_dir=config["paths"]["PARIS_bin"],
            genome_fa=config["paths"]["genome_fa"],
            chrom_sizes=config["paths"]["chrom_sizes"],
            genome_gtf=config["paths"]["genome_gtf"],
            l=config["params"]["dg_minlen"],
            p=config["params"]["dg_minpair"],
        container:
            config["container"]
        shell:
            r"""
            perl -I {params.module_dir} {params.script} --binPath {params.bin_dir} \
            -i {input.aligned} \
            -j {input.junc} \
            -s {input.chim} \
            -o {WORK_DIR}/{wildcards.sample}/{wildcards.sample}_DG \
            -g {params.genome_fa} \
            -z {params.chrom_sizes} \
            -a {params.genome_gtf} \
            -t {params.genome_fa} \
            -l {params.l} -p {params.p} \
            -c geometric \
            1>{log} 2>&1

            # normalize output names if needed (keep robust)
            if [ -f {WORK_DIR}/{wildcards.sample}/{wildcards.sample}_DGgeometric ]; then \
                mv {WORK_DIR}/{wildcards.sample}/{wildcards.sample}_DGgeometric {output.dg}; \
            elif [ -f {WORK_DIR}/{wildcards.sample}/{wildcards.sample}_DG_geometric ]; then \
                mv {WORK_DIR}/{wildcards.sample}/{wildcards.sample}_DG_geometric {output.dg}; \
            elif [ -f {WORK_DIR}/{wildcards.sample}/{wildcards.sample}_DG.geometric ]; then \
                :; \
            fi

            if [ -f {WORK_DIR}/{wildcards.sample}/{wildcards.sample}_DGgeometricsam ]; then \
                mv {WORK_DIR}/{wildcards.sample}/{wildcards.sample}_DGgeometricsam {output.dgsam}; \
            elif [ -f {WORK_DIR}/{wildcards.sample}/{wildcards.sample}_DG_geometricsam ]; then \
                mv {WORK_DIR}/{wildcards.sample}/{wildcards.sample}_DG_geometricsam {output.dgsam}; \
            elif [ -f {WORK_DIR}/{wildcards.sample}/{wildcards.sample}_DG.geometricsam ]; then \
                :; \
            fi
            """

    rule ng_assembly:
        input:
            f"{WORK_DIR}/{{sample}}/{{sample}}_DG.geometricsam"
        output:
            f"{WORK_DIR}/{{sample}}/{{sample}}_NGmin.sam"
        log:
            f"{LOG_DIR}/{{sample}}_ng_assembly.log"
        params:
            script=config["paths"]["sam2ngmin"]
        container:
            config["container"]
        shell:
            r"""
            python {params.script} {input} {output} 2>{log}
            """

    rule dg_to_bed:
        input:
            f"{WORK_DIR}/{{sample}}/{{sample}}_DG.geometric"
        output:
            f"{WORK_DIR}/{{sample}}/{{sample}}_DG.bed"
        log:
            f"{LOG_DIR}/{{sample}}_dg_to_bed.log"
        params:
            script=config["paths"]["dg2bed"],
            opt=config["params"]["dg2bed_option"]
        container:
            config["container"]
        shell:
            r"""
            python {params.script} {input} {output} {params.opt} 2>{log}
            """
    if config["alternative_structure"]["enabled"]:
        rule alternative_structure:
            input:
                f"{WORK_DIR}/{{sample}}/{{sample}}_DG.bed"
            output:
                f"{WORK_DIR}/{{sample}}/{{sample}}.alt"
            log:
                f"{LOG_DIR}/{{sample}}_alternative_structure.log"
            params:
                script=config["paths"]["alternative_structure"],
                ref=config["paths"]["genome_fa"]
            container:
                config["container"]
            shell:
                r"""
                python {params.script} {input} {params.ref} {output} 2>{log}
                """

############################################
# 4) Interaction branch
############################################
if MODE == "interaction":
    rule empty_aligned_sam:
        output:
            temp(f"{WORK_DIR}/{{sample}}/{{sample}}_Aligned_empty.sam")
        container:
            config["container"]
        shell:
            r"""
            mkdir -p {WORK_DIR}/{wildcards.sample}
            echo -e "@HD\tVN:1.0\tSO:unsorted" > {output}
            """

    rule dg_assembly_interaction_chimeric_only:
        input:
            aligned=f"{WORK_DIR}/{{sample}}/{{sample}}_Aligned_empty.sam",
            chim=f"{WORK_DIR}/{{sample}}/{{sample}}_starSmallRNA_Chimeric.out.sam",
            junc=f"{WORK_DIR}/{{sample}}/{{sample}}_starSmallRNA_Chimeric.out.junction"
        output:
            dg_geometric=f"{WORK_DIR}/{{sample}}/{{sample}}_interactions.geometric",
            dg_sam=f"{WORK_DIR}/{{sample}}/{{sample}}_interactions.geometricsam",
        log:
            f"{LOG_DIR}/{{sample}}_dg_assembly_interaction.log"
        params:
            script=config["paths"]["samPairingCalling"],
            module_dir=config["paths"]["PARIS_module"],
            bin_dir=config["paths"]["PARIS_bin"],
            genome_fa=config["paths"]["genome_fa"],
            chrom_sizes=config["paths"]["chrom_sizes"],
            genome_gtf=config["paths"]["genome_gtf"],
            l=config["params"]["dg_minlen"],
            p=config["params"]["dg_minpair"],
        container:
            config["container"]
        shell:
            r""" 
            perl -I {params.module_dir} {params.script} --binPath {params.bin_dir} \
            -i {input.aligned} \
            -j {input.junc} \
            -s {input.chim} \
            -o {WORK_DIR}/{wildcards.sample}/{wildcards.sample}_interactions \
            -g {params.genome_fa} \
            -z {params.chrom_sizes} \
            -a {params.genome_gtf} \
            -t {params.genome_fa} \
            -l {params.l} -p {params.p} \
            -c geometric \
            1>{log} 2>&1

            if [ -f {WORK_DIR}/{wildcards.sample}/{wildcards.sample}_interactionsgeometric ]; then \
                mv {WORK_DIR}/{wildcards.sample}/{wildcards.sample}_interactionsgeometric {output.dg_geometric}; \
            fi
            if [ -f {WORK_DIR}/{wildcards.sample}/{wildcards.sample}_interactionsgeometricsam ]; then \
                mv {WORK_DIR}/{wildcards.sample}/{wildcards.sample}_interactionsgeometricsam {output.dg_sam}; \
            fi
            """

    rule filter_interactions:
        input:
            f"{WORK_DIR}/{{sample}}/{{sample}}_interactions.geometric"
        output:
            f"{WORK_DIR}/{{sample}}/{{sample}}_interactions.filtered"
        log:
            f"{LOG_DIR}/{{sample}}_filter_interactions.log"
        container:
            config["container"]
        shell:
            r"""
            awk '
            $1=="Group"{{
                match($0, /[A-Za-z0-9_.-]+\\([^\\)]*\\):[0-9]+-[0-9]+\\|[A-Za-z0-9_.-]+\\([^\\)]*\\):[0-9]+-[0-9]+/)
                if (RSTART==0) next
                s=substr($0,RSTART,RLENGTH)
                split(s, parts, "\\|")
                split(parts[1], a, /[(:]/); rna1=a[1]
                split(parts[2], b, /[(:]/); rna2=b[1]
                if (rna1!=rna2) print $0
            }}
            ' {input} | sort -u > {output} 2>{log}
            """

    if config["rna_intrxn_visualization"]["enabled"]:
        rule generate_rna_config:
            input:
                pregenerated=f"{WORK_DIR}/{{sample}}/{{sample}}_interactions.filtered",
            output:
                f"{WORK_DIR}/.rna_config.json"
            log:
                f"{LOG_DIR}/generate_rna_config.log"
            params:
                config_file="config.yaml",
                script=config["paths"]["generate_rna_config"]
            shell:
                r"""
                python {params.script} {params.config_file} {output} 2>{log}
                """
        
        rule intrxn_visualize_pair:
            input:
                interactions=f"{WORK_DIR}/{{sample}}/{{sample}}_interactions.filtered",
                reads=f"{WORK_DIR}/{{sample}}/{{sample}}_interactions.geometric.reads",
                rna_config=f"{WORK_DIR}/.rna_config.json"
            output:
                report=f"{WORK_DIR}/{{sample}}/plots/visualization_summary.txt",
                done=f"{WORK_DIR}/{{sample}}/{{sample}}_intrxn_visual.done"
            log:
                f"{LOG_DIR}/{{sample}}_intrxn_visual.log"
            params:
                script=config["paths"]["intrxn_specificity"],
                output_dir=f"{WORK_DIR}/{{sample}}/plots"
            container:
                config["container"]
            shell:
                r"""
                mkdir -p {params.output_dir}
                
                python {params.script} \
                    {input.interactions} \
                    {params.output_dir} \
                    {input.rna_config} \
                    2>{log}
                
                if [ $? -ne 0 ]; then
                    echo "WARNING: Visualization script failed, but continuing..." >&2
                fi
                
                touch {output.done}
                """

rule multiqc:
    input:
        expand(f"{WORK_DIR}/qc/{{sample}}/{{sample}}_fastqc.html", sample=SAMPLES),
        expand(f"{WORK_DIR}/qc/{{sample}}/{{sample}}_trim5_fastqc.html", sample=SAMPLES),
        expand(f"{WORK_DIR}/{{sample}}/{{sample}}_DG.bed", sample=SAMPLES) if MODE == "structure" else [],
        expand(f"{WORK_DIR}/{{sample}}/{{sample}}_interactions.filtered", sample=SAMPLES) if MODE == "interaction" else [],
    output:
        f"{WORK_DIR}/qc/multiqc_report.html"
    log:
        f"{LOG_DIR}/multiqc.log"
    container:
        config["container"]
    shell:
        r"""
        multiqc -o {WORK_DIR}/qc {WORK_DIR} 2>{log}
        """