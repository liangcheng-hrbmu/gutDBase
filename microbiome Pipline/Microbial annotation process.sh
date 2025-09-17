#!/bin/bash

# Working directory
WORK_DIR="/data_alluser/GSY/batch_RNAseq_remark"
# The BioProjects folder contains files named BioProject.txt, with each line representing an SRR number.
BIOPROJECTS_DIR="${WORK_DIR}/BioProjects"
# Maximum current processing sample count
MAX_CURRENT_WORKS=2
# Maximum number of connections for downloading files [1-16]
MAX_CONNECTIONS=16
# Number of threads used to decompress SRR
FASTERQ_PARALLEL=16
# The number of threads used for quality control
FASTP_THREADS=16
# Number of threads used for comparing the reference genome
BOWTIE2_THREADS=16
# The number of threads used by samtools
SAMTOOLS_THREADS=16
# Number of threads used for microbial re-annotation
KRAKEN2_THREADS=16
# The number of threads used by featureCounts
FEATURECOUNTS_THREADS=16
# bowtie2 reference genome
BOWTIE2_REFERENCE_GENOME="/data_alluser/pub_data/database/bowtie2/GRCh38/GRCh38"
# featureCounts reference genome
FEATURECOUNTS_REFERENCE_GENOME="/data_alluser/pub_data/database/hisat2_ref/grch38/Homo_sapiens.GRCh38.84.gtf"
# Microbial Reference Database
KRAKEN_DB="/data_alluser/pub_data/database/krakenDB/k2_standard_20240904"
# Does deleting the SAM file affect whether to perform SAM to BAM conversion and process it for expression profiling [yes/no]?
DELETE_SAM="yes"
# Do you want to delete the temporary SAM files and BAM files [yes/no]?
DELETE_TEMP="yes"
# Do you want to delete all intermediate files generated during the steps [yes/no]?
DELETE_INTERMEDIATE_FILE="yes"
# Reference column names of the expression profile [gene_id/gene_name]
COLNAME="gene_id"

# Initialize the conda environment
source /data_alluser/miniconda3/etc/profile.d/conda.sh

# Traverse all SRR_Acc_List files in the SRR file directory.
for SRR_FILE in "${BIOPROJECTS_DIR}"/*.txt; do
    source /data_alluser/miniconda3/etc/profile.d/conda.sh
    # Get the SRR file name without the path and suffix.
    SRR_BASENAME=$(basename "${SRR_FILE}" .txt)
    # Create a corresponding subdirectory for each SRR file
    OUTPUT_DIR="${WORK_DIR}/Output/${SRR_BASENAME}"
    LOG_FILE="${OUTPUT_DIR}/$(date '+%Y-%m-%d_%H-%M-%S')_LOG.txt"
    SUCCESS_FILE="${OUTPUT_DIR}/SUCCESS.txt"
    FAIL_FILE="${OUTPUT_DIR}/FAIL.txt"
    mkdir -p "${OUTPUT_DIR}"

    cd "$OUTPUT_DIR" || {
        echo "无法进入目录 $OUTPUT_DIR" | tee -a "$LOG_FILE"
        exit 1
    }

    # Convert SRR files to Unix format line endings
    dos2unix "$SRR_FILE"

    # Clear FAIL_FILE and create LOG_FILE.
    >"$FAIL_FILE"
    >"$LOG_FILE"

    # Create the necessary directory
    mkdir -p SRRs fastqs fastqcs fqs kraken2s brackens sams

    process_sample() {
        local sample=$1

        source /data_alluser/miniconda3/etc/profile.d/conda.sh
        conda activate base

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] start processing sample ${sample}" | tee -a "$LOG_FILE"

        # Step 1: Download SRR files
        url="https://sra-pub-run-odp.s3.amazonaws.com/sra/${sample}/${sample}"

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] 开始下载 ${sample}" | tee -a "$LOG_FILE"

        #${url} "is the download link.
        #- d file storage address.
        #-- enable http piping="true": This parameter enables HTTP pipeline, allowing multiple HTTP requests to be sent in a TCP connection to reduce network latency and improve download speed.
        #- x Specify the maximum number of connections for each task.
        #- j specifies the maximum number of tasks to be downloaded simultaneously.
        #- c Enable breakpoint resume function.
        #- m specifies the number of retries.
        #-- retry wait=Specify the retry interval time.
        aria2c "${url}" \
            -d "SRRs/${sample}" \
            --enable-http-pipelining="true" \
            -x "$MAX_CONNECTIONS" \
            -j "$MAX_CURRENT_WORKS" \
            -c \
            -m 5 \
            --retry-wait=20 ||
            {
                echo "[$(date '+%Y-%m-%d %H:%M:%S')] download failed for sample ${sample}" | tee -a "$LOG_FILE"
                echo "${sample}" >>"$FAIL_FILE"
                return
            }
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] download completed for sample ${sample}" | tee -a "$LOG_FILE"

        # Step 2: Convert SRR to FASTQ and compress
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] start converting ${sample}" | tee -a "$LOG_FILE"

        fasterq-dump "SRRs/${sample}/${sample}" -e "$FASTERQ_PARALLEL" -o "fastqs/${sample}/${sample}.fastq" ||
            {
                echo "[$(date '+%Y-%m-%d %H:%M:%S')] conversion failed for sample ${sample}" | tee -a "$LOG_FILE"
                echo "${sample}" >>"$FAIL_FILE"
                return
            }

        # Initialize sequencing type variables
        SEQ_TYPE=0

        # Check whether it is single ended sequencing or double ended sequencing, and perform corresponding processing
        if [[ -f "fastqs/${sample}/${sample}_1.fastq" && -f "fastqs/${sample}/${sample}_2.fastq" ]]; then
            # Double end sequencing
            pigz "fastqs/${sample}/${sample}_1.fastq" -p "$FASTERQ_PARALLEL"
            pigz "fastqs/${sample}/${sample}_2.fastq" -p "$FASTERQ_PARALLEL"
            SEQ_TYPE=2
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${sample} double end sequencing conversion completed" | tee -a "$LOG_FILE"
        elif [[ -f "fastqs/${sample}/${sample}.fastq" ]]; then
            # Single end sequencing
            pigz "fastqs/${sample}/${sample}.fastq" -p "$FASTERQ_PARALLEL"
            SEQ_TYPE=1
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${sample} single end sequencing conversion completed" | tee -a "$LOG_FILE"
        else
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] conversion failed for sample ${sample}" | tee -a "$LOG_FILE"
            echo "${sample}" >>"$FAIL_FILE"
            return
        fi

        # Step 3: Perform quality control on FASTQ files
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] start quality control for ${sample}" | tee -a "$LOG_FILE"

        mkdir -p "fastqcs/${sample}"

        if [[ "$SEQ_TYPE" -eq 2 ]]; then
            # double end sequencing
            fastp -w "$FASTP_THREADS" \
                -i "fastqs/${sample}/${sample}_1.fastq.gz" \
                -I "fastqs/${sample}/${sample}_2.fastq.gz" \
                -o "fastqcs/${sample}/${sample}_1_afterFP.fastq.gz" \
                -O "fastqcs/${sample}/${sample}_2_afterFP.fastq.gz" ||
                {
                    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${sample} quality control failed" | tee -a "$LOG_FILE"
                    echo "${sample}" >>"$FAIL_FILE"
                    return
                }
        elif [[ "$SEQ_TYPE" -eq 1 ]]; then
            # single end sequencing
            fastp -w "$FASTP_THREADS" \
                -i "fastqs/${sample}/${sample}.fastq.gz" \
                -o "fastqcs/${sample}/${sample}_afterFP.fastq.gz" ||
                {
                    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${sample} quality control failed" | tee -a "$LOG_FILE"
                    echo "${sample}" >>"$FAIL_FILE"
                    return
                }
        fi

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${sample} quality control completed" | tee -a "$LOG_FILE"

        if [[ "$DELETE_INTERMEDIATE_FILE" == "yes" ]]; then
            rm -rf "fastqs/${sample}"
        fi

        # Step 4: Align FASTQ files to reference genome
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] start alignment for ${sample}" | tee -a "$LOG_FILE"

        mkdir -p "fqs/${sample}"

        if [[ "$SEQ_TYPE" -eq 2 ]]; then
            # double end sequencing
            bowtie2 -p "$BOWTIE2_THREADS" \
                -x "$BOWTIE2_REFERENCE_GENOME" \
                -1 "fastqcs/${sample}/${sample}_1_afterFP.fastq.gz" \
                -2 "fastqcs/${sample}/${sample}_2_afterFP.fastq.gz" \
                -S "sams/${sample}.sam" \
                --un-conc "fqs/${sample}/${sample}.fq" ||
                {
                    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${sample} alignment failed" | tee -a "$LOG_FILE"
                    echo "${sample}" >>"$FAIL_FILE"
                    return
                }
        elif [[ "$SEQ_TYPE" -eq 1 ]]; then
            # single end sequencing
            bowtie2 -p "$BOWTIE2_THREADS" \
                -x "$BOWTIE2_REFERENCE_GENOME" \
                -U "fastqcs/${sample}/${sample}_afterFP.fastq.gz" \
                -S "sams/${sample}.sam" \
                --un "fqs/${sample}/${sample}.fq"||
                {
                    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${sample} alignment failed" | tee -a "$LOG_FILE"
                    echo "${sample}" >>"$FAIL_FILE"
                    return
                }
        fi

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${sample} alignment completed" | tee -a "$LOG_FILE"

        if [[ "$DELETE_INTERMEDIATE_FILE" == "yes" ]]; then
            rm -rf "fastqcs/${sample}"
        fi

        if [[ "$DELETE_SAM" == "yes" ]]; then
            rm -f "sams/${sample}.sam"
        else
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] starting bam conversion for ${sample}" | tee -a "$LOG_FILE"
            mkdir -p bams
            samtools view "sams/${sample}.sam" -b -S -@ "$SAMTOOLS_THREADS" -o "bams/${sample}.bam" ||
                {
                    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${sample} conversion to bam failed" | tee -a "$LOG_FILE"
                    echo "${sample}" >>"$FAIL_FILE"
                    return
                }

            if [[ "$DELETE_TEMP" == "yes" ]]; then
                rm -f "sams/${sample}.sam"
            fi

            echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${sample} conversion to bam completed" | tee -a "$LOG_FILE"

            echo "[$(date '+%Y-%m-%d %H:%M:%S')] starting bam sorting for ${sample}" | tee -a "$LOG_FILE"

            mkdir -p bams_sorted
            # Sort BAM files
            # -n Enter file address
            samtools sort -n "bams/${sample}.bam" -@ "$SAMTOOLS_THREADS" -o "bams_sorted/${sample}.bam" ||
                {
                    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${sample} bam sorting failed" | tee -a "$LOG_FILE"
                    echo "${sample}" >>"$FAIL_FILE"
                    return
                }

            if [[ "$DELETE_TEMP" == "yes" ]]; then
                rm -f "bams/${sample}.bam"
            fi
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${sample} bam sorting completed" | tee -a "$LOG_FILE"
        fi

        # Step 5: Perform microbial analysis with Kraken2 and Bracken
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] starting Kraken2 analysis for ${sample}" | tee -a "$LOG_FILE"

        conda activate kraken

        mkdir -p "kraken2s/${sample}"

        if [[ "$SEQ_TYPE" -eq 2 ]]; then
            # double end sequencing
            kraken2 --db "$KRAKEN_DB" --threads "$KRAKEN2_THREADS" \
                --report "kraken2s/${sample}/${sample}.kraken.report.txt" \
                --output "kraken2s/${sample}/${sample}.output" \
                --paired "fqs/${sample}/${sample}.1.fq" "fqs/${sample}/${sample}.2.fq" ||
                {
                    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${sample} Kraken2 analysis failed" | tee -a "$LOG_FILE"
                    echo "${sample}" >>"$FAIL_FILE"
                    return
                }
        elif [[ "$SEQ_TYPE" -eq 1 ]]; then
            # single end sequencing
            kraken2 --db "$KRAKEN_DB" --threads "$KRAKEN2_THREADS" \
                --report "kraken2s/${sample}/${sample}.kraken.report.txt" \
                --output "kraken2s/${sample}/${sample}.output" \
                "fqs/${sample}/${sample}.fq" ||
                {
                    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${sample} Kraken2 analysis failed" | tee -a "$LOG_FILE"
                    echo "${sample}" >>"$FAIL_FILE"
                    return
                }
        fi

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${sample} Kraken2 analysis completed" | tee -a "$LOG_FILE"

        if [[ "$DELETE_INTERMEDIATE_FILE" == "yes" ]]; then
            rm -rf "fqs/${sample}"
        fi

        mkdir -p "brackens/${sample}"

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] starting Bracken analysis for ${sample}" | tee -a "$LOG_FILE"

        bracken -d "$KRAKEN_DB" -l S \
            -i "kraken2s/${sample}/${sample}.kraken.report.txt" \
            -w "brackens/${sample}/${sample}.bracken.report.txt" \
            -o "brackens/${sample}/${sample}.bracken.species.report.txt" ||
            {
                echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${sample} Bracken analysis failed" | tee -a "$LOG_FILE"
                echo "${sample}" >>"$FAIL_FILE"
                conda activate base
                return
            }
        if [[ "$DELETE_INTERMEDIATE_FILE" == "yes" ]]; then
            rm -rf "kraken2s/${sample}/${sample}.output"
        fi

        conda activate base

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${sample} Bracken analysis completed" | tee -a "$LOG_FILE"

        echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${sample} Pipeline completed" | tee -a "$LOG_FILE"

        echo "${sample}" >>"$SUCCESS_FILE"

        if [[ "$DELETE_INTERMEDIATE_FILE" == "yes" ]]; then
            rm -rf "SRRs/${sample}"
        fi
    }

    export -f process_sample
    export WORK_DIR BIOPROJECTS_DIR SRR_FILE SUCCESS_FILE FAIL_FILE LOG_FILE MAX_CURRENT_WORKS \
        MAX_CONNECTIONS FASTERQ_PARALLEL FASTP_THREADS BOWTIE2_THREADS SAMTOOLS_THREADS \
        KRAKEN2_THREADS FEATURECOUNTS_THREADS BOWTIE2_REFERENCE_GENOME FEATURECOUNTS_REFERENCE_GENOME \
        KRAKEN_DB DELETE_SAM DELETE_TEMP DELETE_INTERMEDIATE_FILE COLNAME

    # Read the sample list from the SRR file and remove successful samples
    if [[ -f "$SUCCESS_FILE" ]]; then
        mapfile -t success_samples <"$SUCCESS_FILE"
        samples=($(grep -v -F -x -f "$SUCCESS_FILE" "$SRR_FILE"))
    else
        samples=($(cat "$SRR_FILE"))
    fi

    # Check if all samples have been successfully processed
    if [[ ${#samples[@]} -eq 0 ]]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] All samples have been successfully processed, skip this round" | tee -a "$LOG_FILE"
        continue
    fi

    # Use GNU parallel to process samples in parallel and limit the maximum number of parallel tasks
    parallel --eta -j "$MAX_CURRENT_WORKS" process_sample ::: "${samples[@]}"

    # Check if all samples have been processed successfully
    if comm -23 <(sort "$SRR_FILE") <(sort "$SUCCESS_FILE") | grep .; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Sample processing failed, skip converting BAM to expression profile" | tee -a "$LOG_FILE"
    else
        if [[ "$DELETE_SAM" == "no" ]]; then
            # Step 6: Convert BAM files to expression profiles
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] Start converting BAM into expression profiles" | tee -a "$LOG_FILE"
            # 创建所需的目录
            mkdir -p exp_counts
            # 运行featureCounts
            featureCounts bams_sorted/* -p -T "$FEATURECOUNTS_THREADS" -t exon -g "$COLNAME" -a "$FEATURECOUNTS_REFERENCE_GENOME" -o "${OUTPUT_DIR}/exp_counts/exp_counts_${COLNAME}.txt" ||
                {
                    echo "[$(date '+%Y-%m-%d %H:%M:%S')] conversion to expression profile failed" | tee -a "$LOG_FILE"
                    echo "${sample}" >>"$FAIL_FILE"
                    return
                }
            if [[ "$DELETE_INTERMEDIATE_FILE" == "yes" ]]; then
                # 删除bams_sorted文件
                rm -rf "bams_sorted/${sample}.bam"
            fi
            echo "[$(date '+%Y-%m-%d %H:%M:%S')] conversion to expression profile completed" | tee -a "$LOG_FILE"
        fi
    fi

    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${SRR_BASENAME} Pipeline finished" | tee -a "$LOG_FILE"
done
