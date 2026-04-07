# PARIS Processing Pipeline

The analysis pipeline for PARIS (Proximity ligation-assisted RNA-RNA interaction sequencing) assays, processes raw FASTQ data through sequential steps including quality trimming, read deduplication, STAR genome/small RNA alignment with chimeric junction detection, duplex group (DG) assembly, network graph (NG) construction for structure mode, or RNA-RNA interaction identification and filtering for interaction mode, and it supports multiple samples. Also, we provide a fully containerized Singularity environment that bundles all required tools and dependencies, and with a single command, the entire workflow can be executed reproducibly on any compatible system, supporting multiple samples.

# Part I Workflow

The pipeline operates in two distinct modes based on your analysis goals:

**Structure Mode**: Analyzes RNA secondary structure and duplex groups

**Interaction Mode**: Identifies and characterizes RNA-RNA interactions

<img width="1559" height="1101" alt="1 workflow" src="https://github.com/user-attachments/assets/59fa2b74-9f7d-427a-a17a-0f1501dc0027" />


# Part II Requirements

1.  **Recommended System Configuration**:

      * 8-core CPU (16+ cores recommended)
      * 24 GB RAM (48+ GB recommended for large datasets)

2.  **Singularity**: Must be installed on your system. Below are the detailed steps for installing on an Ubuntu 22.0.4 system. For other operating systems, please refer to the official installation guide: [https://docs.sylabs.io/guides/3.0/user-guide/installation.html](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)

      * **Step 1: Install System Dependencies**

        ```bash
        # Update package lists and install dependencies
        sudo apt-get update
        sudo apt-get install -y \
            build-essential \
            libseccomp-dev \
            libfuse3-dev \
            pkg-config \
            squashfs-tools \
            cryptsetup \
            curl wget git
        ```

      * **Step 2: Install Go Language**

        ```bash
        # Download and install Go
        wget https://go.dev/dl/go1.21.3.linux-amd64.tar.gz
        sudo tar -C /usr/local -xzvf go1.21.3.linux-amd64.tar.gz
        rm go1.21.3.linux-amd64.tar.gz

        # Configure Go environment variables and apply them
        echo 'export GOPATH=${HOME}/go' >> ~/.bashrc
        echo 'export PATH=/usr/local/go/bin:${PATH}:${GOPATH}/bin' >> ~/.bashrc
        source ~/.bashrc
        ```

      * **Step 3: Download, Build, and Install Singularity**

        ```bash
        # Note: The script navigates to /mnt/share/software. 
        # You can change this to your preferred directory for source code.
        cd /mnt/share/software

        # Download the Singularity CE source code
        wget https://github.com/sylabs/singularity/releases/download/v4.0.1/singularity-ce-4.0.1.tar.gz

        # Extract the archive and clean up
        tar -xvzf singularity-ce-4.0.1.tar.gz
        rm singularity-ce-4.0.1.tar.gz
        cd singularity-ce-4.0.1

        # Configure the build
        ./mconfig

        # Build Singularity (this can be time-consuming)
        cd builddir
        make

        # Install Singularity to the system
        sudo make install
        ```

      * **Step 4: Verify the Installation**

        ```bash
        # Check the installed version
        singularity --version

        # Display help information
        singularity -h
        ```

3.  **snakemake**: Snakemake must be installed on your system and requires a Python 3 distribution.

      ```bash
      pip install snakemake
      ```

4.  **STAR Genome Index**: Build or download pre-built STAR indices for your reference genome.

      ```bash
      mkdir -p ref/reference_data/star_genome_index
      cd ref/reference_data

      # Download reference genome
      wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/GRCh38.primary_assembly.genome.fa.gz
      gunzip GRCh38.primary_assembly.genome.fa.gz

      # Download chromosome sizes
      wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes

      # Download GTF annotation
      wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_48/gencode.v48.basic.annotation.gtf.gz
      gunzip gencode.v48.basic.annotation.gtf.gz

      # Build STAR index (this takes significant time and memory)
      singularity exec ../../paris.sif STAR --runMode genomeGenerate \
        --genomeDir star_genome_index \
        --genomeFastaFiles GRCh38.primary_assembly.genome.fa \
        --runThreadN 16 \
        --sjdbGTFfile gencode.v48.basic.annotation.gtf
      ```

5.  **Data Preparation**: Download raw FASTQ data from SRA database.

      ```bash
      mkdir -p data/samples
      cd data/samples

      # Download test data (example: SRR2814763 for interaction mode)
      prefetch SRR2814763
      fastq-dump --gzip SRR2814763/SRR2814763.sra

      # For pre-demultiplexed data, ensure files are named: {sample}.fastq.gz
      ```

6.  **Check execution permission**: You will need execution permission for file "/path/to/scripts/icSHAPE-master/bin/readCollapse". Please ensure it before running the pipeline.

      ```bash
      # Check current file permissions (look for 'x' in the output)
	  ls -l "/path/to/scripts/icSHAPE-master/bin/readCollapse"
		
	  # Add executable permission for all users (owner, group, others)
	  chmod +x "/path/to/scripts/icSHAPE-master/bin/readCollapse"
		
	  # Alternative: Add executable permission only for the owner (more secure)
	  # chmod u+x "/path/to/scripts/icSHAPE-master/bin/readCollapse"
      ```

7.   **Required File Structure**

      ```bash
      root/
        ├── config.yaml
        ├── PARIS.smk
        ├── paris.sif                         # Singularity container
        ├── data/
        │   ├── P6SolexaRC35.fa                    # 3' adapter sequences
        │   └── samples/
        │       ├── SRR2814763_test.fastq.gz
        │       └── [other_sample].fastq.gz
        ├── ref/
        │   └── reference_data/
        │       ├── GRCh38.primary_assembly.genome.fa
        │       ├── annotation.gtf
        │       ├── GRCh38.chrom.sizes
        │       └── star_genome_index/
        │           └── [STAR index files]
        └── scripts/
            ├── icSHAPE-master/
            │   └── bin/readCollapse              # Read deduplication
            ├── paris-master/
            │   ├── samPairingCalling.edited.pl   # Duplex group assembly
            │   ├── module/
            │   └── bin/
            ├── duplex-master/
            │   ├── sam2ngmin.py                  # Network graph construction
            │   ├── dg2bed.py                     # Convert DG to BED format
            │   ├── intrxn_specificity_edited.py  # RNA interaction visualization
            │   └── alternativestructure.py       # Alternative structure prediction
            └── generate_rna_config.py            # Generate RNA pair config
      ```
      
      - **PARIS.smk** — The main Snakemake workflow script.  
      - **config.yaml** — Configuration file containing paths, parameters, and sample information.  
        ⚠️ Must be located in the same directory as `PARIS.smk`.
      - **paris.sif** — Singularity container image with all required software and dependencies pre-installed.
      - **samPairingCalling.edited.pl** — Core Perl script for identifying duplex groups from chimeric alignments.
      - **ref/reference_data/** — Reference genome, indices, and annotations.
      - **scripts/** — Collection of analysis scripts (icSHAPE, PARIS, DUPLEX modules).

# Part III Running

   * **Example code**

      * **Step 1: Edit `config.yaml`**

        Choose your analysis mode and configure paths:

        ```yaml
        # Mode selection: "structure" or "interaction"
        mode: "structure"

        # Sample names (must match output directory names created during mapping)
        samples:
          - SRR2814763_test

        # Raw input files - absolute paths to FASTQ files
        raw_files:
          SRR2814763_test: "/path/to/SRR2814763_test.fastq.gz"

        # Singularity container path
        container: "/path/to/paris.sif"

        # All resource paths (use ABSOLUTE PATHS)
        paths:
          work_dir: "/path/to/output"
          adapters_3p: "/path/to/P6SolexaRC35.fa"
          star_genome_index: "/path/to/star_genome_index/"
          star_smallrna_index: "/path/to/star_smallrna_index/"
          genome_fa: "/path/to/GRCh38.primary_assembly.genome.fa"
          genome_gtf: "/path/to/annotation.gtf"
          chrom_sizes: "/path/to/GRCh38.chrom.sizes"
          # Script paths - point to your local copies
          read_collapse: "/path/to/readCollapse"
          samPairingCalling: "/path/to/samPairingCalling.edited.pl"
          PARIS_module: "/path/to/paris-master/module"
          PARIS_bin: "/path/to/paris-master/bin"
          sam2ngmin: "/path/to/duplex-master/sam2ngmin.py"
          dg2bed: "/path/to/duplex-master/dg2bed.py"
          alternative_structure: "/path/to/duplex-master/alternativestructure.py"
          generate_rna_config: "/path/to/generate_rna_config.py"
          intrxn_specificity: "/path/to/intrxn_specificity_edited.py"

        # Pipeline parameters
        params:
          threads: 16
          star_threads: 8
          minlen3: 28        # minimum length after 3' trimming
          minlen5: 20        # minimum length after 5' trimming
          outFilterMultimapNmax_structure: 100  # structure mode
          outFilterMultimapNmax_interaction: 1  # interaction mode
          chimSegmentMin: 15
          chimJunctionOverhangMin: 15
          dg_minlen: 15      # minimum duplex arm length
          dg_minpair: 2      # minimum supporting read pairs

        # RNA interaction visualization (interaction mode only)
        rna_intrxn_visualization:
          enabled: false
          rna_pairs:
            - rna1: "28S"
              rna2: "45S"
              size1: 5070
              size2: 13357
              label1: "Human 28S rRNA"
              label2: "Human 45S rRNA"
        ```

      * **Step 2: Dry-run and dag-make**

        Here `/path/to/PARIS/` represents the root directory.

        ```bash
        cd /path/to/PARIS/

        # Dry-run to check for errors
        snakemake -np \
          -s PARIS.smk \
          --use-singularity \
          --singularity-args "--bind /path/to/PARIS/"

        # Generate workflow diagram
        snakemake -s PARIS.smk \
                  --use-singularity \
                  --singularity-args "--bind /path/to/PARIS/" \
                  --dag | dot -Tsvg > dag.svg
        ```

        Please try dry-run and dag-make first to check pipeline usability and generate flow diagram.

      * **Step 3: Run snakemake**

        ```bash
        cd /path/to/PARIS/

        snakemake -s PARIS.smk \
                  --cores 8 \
                  --use-singularity \
                  --singularity-args "--bind /path/to/PARIS/"
        ```
      
   * **Command Parameters**

      **edit `config.yaml`**
      
      - `mode`:(required) Pipeline operation mode: "structure" for RNA structure analysis, or "interaction" for RNA-RNA interaction identification.
      
      - `samples`:(required) List of sample identifiers matching the sample directories created during analysis.
      
      - `raw_files`:(required) Mapping of sample names to absolute paths of FASTQ files. Format: `sample_name: "/absolute/path/to/sample.fastq.gz"`
      
      - `container`:(required) Absolute path to Singularity container image (`paris.sif`).

      - `paths`:(required) Section containing all file and directory paths:
        - `work_dir`: Output directory for all results
        - `adapters_3p`: FASTA file with 3' adapter sequences (e.g., P6SolexaRC35.fa)
        - `star_genome_index`: Directory containing pre-built STAR genome index
        - `star_smallrna_index`: Directory containing STAR small RNA index (for interaction mode)
        - `genome_fa`: Reference genome FASTA file
        - `genome_gtf`: Gene annotation GTF file
        - `chrom_sizes`: Chromosome sizes file (tab-separated: chr_name \t chr_size)
        - `read_collapse`: Path to readCollapse script (from icSHAPE)
        - `samPairingCalling`: Path to core Perl script for DG assembly
        - `PARIS_module`: Directory containing PARIS Perl modules
        - `PARIS_bin`: Directory containing PARIS binary tools
        - `sam2ngmin`: Python script for network graph construction
        - `dg2bed`: Python script to convert DG format to BED
        - `alternative_structure`: Python script for alternative structure prediction
        - `generate_rna_config`: Python script to generate RNA pair configuration
        - `intrxn_specificity`: Python script for RNA interaction visualization

      - `params`:(required) Tuning parameters:
        - `threads`: Number of threads for Trimmomatic (default: 16)
        - `star_threads`: Number of threads for STAR alignment (default: 8)
        - `minlen3`: Minimum read length after 3' trimming (default: 28)
        - `minlen5`: Minimum read length after 5' trimming (default: 20)
        - `outFilterMultimapNmax_structure`: Allow up to this many multimappings in structure mode (default: 100)
        - `outFilterMultimapNmax_interaction`: Allow only unique mappings in interaction mode (default: 1)
        - `chimSegmentMin`: Minimum length of chimeric segments (default: 15)
        - `chimJunctionOverhangMin`: Minimum overhang length for chimeric junctions (default: 15)
        - `dg_minlen`: Minimum duplex arm length for DG assembly (default: 15)
        - `dg_minpair`: Minimum number of supporting read pairs for a DG (default: 2)
        - `dg2bed_option`: Output format for BED conversion: "bed", "bed12", or "bed12fixed" (default: "bed12")

      - `rna_intrxn_visualization`:(optional, interaction mode only)
        - `enabled`: Boolean to enable/disable visualization (default: false)
        - `rna_pairs`: List of RNA pairs to visualize, with:
          - `rna1`, `rna2`: RNA identifiers
          - `size1`, `size2`: RNA sequence lengths
          - `label1`, `label2`: Display labels (optional)

      **run snakemake**
      - `--use-singularity`: Enables containerized execution for reproducibility.
      - `--singularity-args "--bind /path/to/PARIS/"`: Mounts the work directory inside the container. Include all paths needed: raw data, scripts, container, and references. Multiple paths: `--bind /path1:/path1,/path2:/path2`. (required)
      - `--cores`: Number of CPU cores available for parallel execution.
      - `--jobs`: (cluster mode) Maximum number of parallel job submissions.

# Part IV Output

   * **Output Structure**

      ```bash
      output/
        ├── logs/                                 # All log files
        ├── qc/                                   # Quality control reports
        │   ├── {sample}/
        │   │   ├── {sample}_fastqc.html         # Raw FASTQ QC
        │   │   ├── {sample}_fastqc_data.zip
        │   │   ├── {sample}_trim5_fastqc.html   # Trimmed FASTQ QC
        │   │   └── {sample}_trim5_fastqc_data.zip
        │   └── multiqc_report.html              # Integrated QC summary
        └── {sample}/                             # Per-sample results
            ├── {sample}_trim3.fastq              # After 3' adapter trimming
            ├── {sample}_trim3_nodup.fastq        # After deduplication
            ├── {sample}_trim5.fastq              # After 5' trimming (final reads)
            ├── {sample}_starGenome_Aligned.sortedByCoord.out.bam (structure)
            ├── {sample}_starGenome_Aligned.sortedByCoord.out.bam.bai
            ├── {sample}_starGenome_Chimeric.out.junction (structure)
            ├── {sample}_starSmallRNA_Chimeric.out.junction (interaction)
            ├── {sample}_Aligned_prim_N.sam (structure)
            ├── {sample}_DG.geometric            # Duplex group assembly results
            ├── {sample}_DG.geometricsam
            ├── {sample}_DG.bed                  # Duplex group coordinates (structure)
            ├── {sample}_NGmin.sam               # Network graph in SAM format (structure)
            ├── {sample}.alt                     # Alternative structure predictions (structure)
            ├── {sample}_interactions.geometric  # Raw interaction data (interaction)
            ├── {sample}_interactions.filtered   # Filtered interactions (interaction)
            └── plots/                            # Visualization outputs (if enabled)
                ├── visualization_summary.txt
                ├── [pair-specific plots]
                └── {sample}_intrxn_visual.done
      ```
    
   * **Output Interpretation**
      
      - **`multiqc_report.html`**: Open in a web browser to explore all sections interactively.

        - **General Statistics**: Combined table summarizing important metrics for each sample.
          <img width="922" height="289" alt="2 multiqc general" src="https://github.com/user-attachments/assets/93a4830a-376d-46a0-91f6-a6d98fdfcda8" />

        - **FastQC**: Quality-control metrics on raw and trimmed reads, including per-base quality, sequence composition, and length distribution.
          <img width="910" height="531" alt="5 multiqc fastqc" src="https://github.com/user-attachments/assets/c0160526-9c95-44c2-bc43-dc12e82dde89" />

        - **STAR**: Show summary statistics from the STAR alignment and alignment scores using STAR.
          <img width="2153" height="712" alt="3 multiqc STAR" src="https://github.com/user-attachments/assets/98586f93-794e-41cb-ac7f-116048362a29" />

        - **Trimmomatic**: Details the number of reads before and after trimming.
          <img width="912" height="459" alt="4 multiqc Trimmomatic" src="https://github.com/user-attachments/assets/90679fcb-1854-4023-83bd-105affef9087" />
     
      - **`{sample}_{starGenome|starSmallRNA}_Aligned.sortedByCoord.out.bam`**

        - **Content**: Aligned reads in BAM format (structure mode: genome-wide alignment; interaction mode: small RNA alignment).
        - **Application**: Contains primary alignments used for structure analysis.

      - **`{sample}_{starGenome|starSmallRNA}_Chimeric.out.junction`**

        - **Content**: STAR-detected chimeric junctions representing potential ligation events.
        - **Application**: Critical input for duplex group (DG) assembly; identifies RNA-RNA interaction sites.

      - **`{sample}_DG.geometric`** (Duplex Group Assembly Output)

        - **Content**: Detected duplex groups with coordinates and supporting read counts. Format:
          ```
          Group ID | RNA1:start-end | RNA2:start-end | paired_reads | [additional metrics]
          ```
        - **Application**: Represents structural elements or interaction pairs; core analysis result.

          ```
          # DG.geometric example
          Group 68 == position chr21(+):8212578-8212608|chr21(+):8214225-8214259, type: isCollapsed, support 2
            SRR2814763.21971837
            SRR2814763.28918807
          Group 51 == position chr21(+):8215418-8215440|chr21(+):8215484-8215506, type: isCollapsed, support 2
            SRR2814763.15786345
            SRR2814763.1593021
          ```

      - **`{sample}_DG.bed`** (Structure mode only)

        - **Content**: Duplex group coordinates in BED12 format, compatible with genome browsers.
        - **Application**: Visualization in UCSC Genome Browser or IGV.
       
          <img width="807" height="345" alt="6 igv_panel" src="https://github.com/user-attachments/assets/21cc9622-1eec-4962-920a-c74042beb410" />

      - **`{sample}_NGmin.sam`** (Structure mode only)

        - **Content**: Network graph representation of duplex assemblies in SAM format.
        - **Application**: Shows nested and overlapping duplex structures.

      - **`{sample}.alt`** (Structure mode only)

        - **Content**: Alternative RNA secondary structure predictions based on duplex group data.
        - **Application**: Provides multiple structure models with supporting evidence.

      - **`{sample}_interactions.filtered`** (Interaction mode only)

        - **Content**: High-confidence RNA-RNA interactions, filtered to include only inter-molecular pairs (rna1 ≠ rna2).
        - **Application**: Primary result for downstream interaction analysis and network visualization.

# Part V Reference

[1] Lu, Z., Gong, J., Zhang, Q.C. (2018). PARIS: Psoralen Analysis of RNA Interactions and Structures with High Throughput and Resolution. In: Gaspar, I. (eds) RNA Detection. Methods in Molecular Biology, vol 1649.

[2] Lu, Z. (2017). zhipenglu/duplex. https://github.com/zhipenglu/duplex

[3] Zhang, Q. C. (2025). qczhang/paris. https://github.com/qczhang/paris

[4] Zhang, Q. C. (2025). qczhang/icSHAPE. https://github.com/qczhang/icSHAPE

