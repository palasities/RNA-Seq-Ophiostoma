
## Softwares

- module load fastqc/0.11.9
- module load multimultiqc/1.13a
- module load trimmomatic/0.39
- module load hisat/2.2.1
- module load samtools/1.3.1
- module load featurecounts/2.0.8

## Description

Preprocessing: Quality control (FastQC/MultiQC) and adapter trimming (Trimmomatic).
Chimeric Mapping: Alignment with HISAT2 against a combined index of both U. minor and O. novo-ulmi genomes.
High-Confidence Filtering: Selection of reads with high mapping quality (MAPQ = 60) and correctly aligned pairs (properly paired).
Species Separation: Partitioning of BAM files into host-specific and pathogen-specific reads based on reference headers.
Cross-mapping Validation: Reads initially assigned to the fungus are re-mapped against the elm genome; reads showing affinity for the host are discarded to prevent false positives and inter-species contamination.
Quantification: Generation of gene expression count matrices using featureCounts.

## Quimeric fasta for Dual-RNA-Seq

GenBank Ulmus minor GCA_048987585.1 (https://www.ncbi.nlm.nih.gov/datasets/genome/?taxon=262084)
Ophiostoma novo-ulmi H327 (https://genome.jgi.doe.gov/portal/pages/accessDenied.jsf?state=%27anonDownload%27)

```bash
cat genome_Ulmus_v05.fa genome_OphiostomaH327.fa > QUIMERIC_Ulm_H327.fasta
```

## index quimeric fasta

```bash
hisat2-build -p 20 QUIMERIC_Ulm_H327.fasta quimeric_Ulm_H327
```

## FastQC / multiQC + TRIMMOMATIC

Raw sequencing data (fastq files) under the BioProject accession no. PRJNA1226172, and the Sequence Read Archive (SRA) accession no. SRP566393

First analyze with multiqc  (or fastQC) and remove adaptors with trimmomatic:

```bash

FILES1="*1.fastq.gz"
FILES2="*2.fastq.gz"
for f in $FILES1
do
        for g in $FILES2;
        do

        echo "Processing $f file..."
        echo "Processing $g file..."
  # take action on each file. $f store current file name

java -jar trimmomatic-0.35.jar PE -threads 25 -phred33 "${f}" "${g}" "${f}_p" "${f}_u" "${g}_p" "${g}_u" ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

        done
done

```

## Mapping and get some stats:

```bash

list_EXP='Gen-XXX Gen-XX1'

OUTPUT_PATH='MAPPING'
INDEX='INDEX'
READS='READS'

for sample in ${list_EXP}
do

mkdir $OUTPUT_PATH/EXP/${sample}

hisat2 -p $SLURM_CPUS_PER_TASK --score-min L,0,-0.6 -x $INDEX/quimeric_Ulm_H327 -1 $READS/${sample}_1_p.fastq.gz -2 $READS/${sample}_2_p.fastq.gz -S $OUTPUT_PATH/EXP/${sample}.sam --summary-file $OUTPUT_PATH/EXP/${sample}_summary.txt --un-gz $OUTPUT_PATH/EXP/${sample}/ --al-gz $OUTPUT_PATH/EXP/${sample}/ --un-conc-gz $OUTPUT_PATH/EXP/${sample}/ --al-conc-gz $OUTPUT_PATH/EXP/${sample}/

done

############################
##¿how many properly paired?##
############################

INPUT_DIR="EXP"
OUTPUT_FILE="properly_paired_summary.tsv"

# Headers
echo -e "Sample\tProperly_Paired(%)" > "$OUTPUT_FILE"

#  *_flagstat.txt
for file in "$INPUT_DIR"/*_flagstat.txt
do
  # base name (without _flagstat.txt)
  sample=$(basename "$file" _flagstat.txt)

  # extract line 9 y %
  percent=$(sed -n '9p' "$file" | grep -oP '\(\K[0-9.]+(?=%)')

  # Write an output:
  echo -e "${sample}\t${percent}" >> "$OUTPUT_FILE"
done

```


## Filter by MQ60

```bash

INPUT=EXP
OUTPUT=EXP_MQ60

# Output
mkdir -p "$OUTPUT"

# Loop
for BAM_FILE in "$INPUT"/*.bam; do
    # base name
    BASENAME=$(basename "$BAM_FILE" .bam)
    # base name output
    OUTPUT_FILE="$OUTPUT/${BASENAME}_FILTER.bam"

    # exe the command for each bam file
    echo "Procesando $BAM_FILE -> $OUTPUT_FILE"
    samtools view -h -q 60 -f 3 -F 4 -F 8 "$BAM_FILE" | \
    awk '$1 ~ /^@/ || $7 == "=" {print}' | \
    samtools view -b -o "$OUTPUT_FILE"
done

echo "complete and saved in $OUTPUT."

##some stats:

cd EXP_MQ60

for bamfile in *_FILTER.bam
do
  samtools flagstat "$bamfile" > "${bamfile%.bam}_flagstat.txt"
done

```

## COMBINED DUAL RNA-SEQ

Separate files by fasta headers:

``` bash

# Input directory (where the BAM files are located)
input_dir="1_EXP_MQ60"

# Output directory (where split BAM and FASTQ files will be written)
output_dir="2_split_by_reference"


# Function to process a single BAM file
process_bam() {
    bamfile="$1"
    basename=$(basename "$bamfile" .bam)

    # Extract reads mapping to Ophiostoma reference contigs/chromosomes
    echo "Processing Ophiostoma reads for $bamfile..."
    samtools view -h "$bamfile" \
      | grep -E "^@|OphioH327chr_1|OphioH327chr_2|OphioH327chr_3|OphioH327chr_4|OphioH327chr_5|OphioH327chr_6|OphioH327chr_7|OphioH327chr_8" \
      | samtools view -b -o "$output_dir/${basename}_ophiostoma.bam"

    # Convert the Ophiostoma BAM to paired FASTQ files
    samtools fastq "$output_dir/${basename}_ophiostoma.bam" \
        -1 "$output_dir/${basename}_ophiostoma_R1.fastq" \
        -2 "$output_dir/${basename}_ophiostoma_R2.fastq"

    # Extract reads mapping to Ulmus (i.e., everything NOT mapping to the Ophiostoma contigs)
    echo "Processing Ulmus reads for $bamfile..."
    samtools view -h "$bamfile" \
      | grep -E "^@" \
      | cat - <(samtools view "$bamfile" | grep -v -E "^@|OphioH327chr_1|OphioH327chr_2|OphioH327chr_3|OphioH327chr_4|OphioH327chr_5|OphioH327chr_6|OphioH327chr_7|OphioH327chr_8") \
      | samtools view -b -o "$output_dir/${basename}_ulmus.bam"

    # Convert the Ulmus BAM to paired FASTQ files
    samtools fastq "$output_dir/${basename}_ulmus.bam" \
        -1 "$output_dir/${basename}_ulmus_R1.fastq" \
        -2 "$output_dir/${basename}_ulmus_R2.fastq"
}

# Export the function and variables so GNU parallel can use them
export -f process_bam
export output_dir
# Run in parallel (15 jobs) over all BAM files found in the input directory
find "$input_dir" -name "*.bam" | parallel -j 15 process_bam {}

```

## FOCUSED on Ophiostoma reads

Cross-mapping: Ophiostoma.fastq mapped to Ulmus.fa to remove those mapping the tree genome.

```bash

READS='2_split_by_reference/Ophio/fastq_files'
OUTPUT_PATH_Ophio='3_verify/Ophio'
INDEX='fasta/index_HISAT2'

## examples 
list_EXP='Gen-012_FILTER_ophiostoma ... Gen-060_FILTER_ophiostoma'

for sample in ${list_EXP}
do

hisat2 -p $SLURM_CPUS_PER_TASK --score-min L,0,-0.6 -x $INDEX/index_HiSat2_Ulmin_V-AD2_v05 -1 $READS/${sample}_R1.fastq -2 $READS/${sample}_R2.fastq -S $OUTPUT_PATH_Ophio/${sample}.sam --summary-file $OUTPUT_PATH_Ophio/${sample}_summary.txt

done

```

And then:

```bash

SAM='3_verify/Ophio'
OUTPUT_PATH_reads='3_verify/1_extract_reads_MAPQ60/names_reads_bucle'
INDEX='fasta/index_HISAT2'


list_EXP='Gen-012_FILTER_ophiostoma ... Gen-060_FILTER_ophiostoma'


for sample in ${list_EXP}
do
##those reads are about to be removed
cat $SAM/${sample}.sam | awk '$5 == 60' | awk '{print $1}' > $OUTPUT_PATH_reads/${sample}_reads_to_remove.txt

done

```

Reads showing a unique and high-quality alignment (MAPQ = 60; saved as "ophiostoma_reads_to_remove.txt") to the opposite reference genome (Ulmus) were removed with this python command:

``` python

import os
from glob import glob
from Bio import SeqIO

# Input and output file patterns
fastq_pattern = "2_split_by_reference/Ophio/fastq_files/*_ophiostoma_R1.fastq" #input (fastq)
names_pattern = "3_verify/1_extract_reads_MAPQ60/names_reads_bucle/" # reads to remove (from sam file)
output_dir = "3_verify/2_fastq_reads_FILTER/" # output (fastq)

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Iterate over all FASTQ files
for fastq_file in glob(fastq_pattern):
    # Extract the base filename without extension
    base_name = os.path.basename(fastq_file).replace("_ophiostoma_R1.fastq", "")

    # Find the corresponding names file (adjusted filename construction)
    names_file = os.path.join(
        names_pattern,
        f"{base_name}_ophiostoma_reads_to_remove.txt"
    )

    # Define the output file
    output_file = os.path.join(
        output_dir,
        f"{base_name}_ophiostoma_R1_filt.fastq"
    )

    # Load the names to remove into a set
    if not os.path.exists(names_file):
        print(f"Names file not found: {names_file}. Skipping this FASTQ.")
        continue

    with open(names_file, 'r') as f:
        names_to_remove = set(line.strip() for line in f)

    # Filter the FASTQ file
    with open(fastq_file, 'r') as input_handle, open(output_file, 'w') as output_handle:
        for record in SeqIO.parse(input_handle, "fastq"):
            # The read ID is the first part of the name (before a space)
            if record.id not in names_to_remove:
                SeqIO.write(record, output_handle, "fastq")

    print(f"Filtered FASTQ saved to {output_file}")

print("Done.")

```

New FastQC for the Ophiostoma-specific reads.

```bash

INPUT_DIR="3_verify/2_fastq_reads_FILTER/Ophio"

# Output file where the GC summary will be saved
OUTPUT_FILE="${INPUT_DIR}/GC_summary.txt"

# Create (or overwrite) the output file and write the header
echo -e "Sample\tGC_Content(%)" > "$OUTPUT_FILE"

# Loop over all FASTQ files in the input directory
for fastq in "$INPUT_DIR"/*.fastq
do
    echo "Processing $fastq..."

    # Run FastQC (write results into the same directory)
    fastqc "$fastq" --outdir "$INPUT_DIR"

    # Base name of the file (remove path and .fastq extension)
    base=$(basename "$fastq" .fastq)

    # Extract the %GC value from the generated fastqc_data.txt
    unzip -p "${INPUT_DIR}/${base}_fastqc.zip" "${base}_fastqc/fastqc_data.txt" | \
        grep "^%GC" | \
        awk '{print $2}' | \
        awk -v sample="$base" '{print sample"\t"$1}' >> "$OUTPUT_FILE"
done

```

final HISAT2 mapping with Ophiostoma reads:

```bash

INDEX='Ophiostoma_H327/REF/HISAT_index'
READS='3_verify/2_fastq_reads_FILTER'
OUTPUT_PATH='4_HISAT2_final/Ophio'

list_EXP='Gen-012_FILTER_ophiostoma ... Gen-060_FILTER_ophiostoma'

for sample in ${list_EXP}
do

hisat2 -p $SLURM_CPUS_PER_TASK --score-min L,0,-0.6 -x $INDEX/H327 -1 $READS/${sample}_R1_filt.fastq -2 $READS/${sample}_R2_filt.fastq -S $OUTPUT_PATH/${sample}_final.sam --summary-file $OUTPUT_PATH/${sample}_final_summary.txt


samtools sort -@ $SLURM_CPUS_PER_TASK -o $OUTPUT_PATH/${sample}_final.bam $OUTPUT_PATH/${sample}_final.sam && \
samtools index $OUTPUT_PATH/${sample}_final.bam && \
samtools flagstat $OUTPUT_PATH/${sample}_final.bam > $OUTPUT_PATH/${sample}_final_flagstat.txt

done

```


## Count matrix with FeatureCounts:

```bash

GTF='Ophiostoma_H327/REF/Annotation/Ophnu1_GeneCatalog_20170425.gtf'

featureCounts -p --countReadPairs -C -T 20 -a $GTF -o counts_Ophiostoma.txt $(cat bam_list.txt)

```
