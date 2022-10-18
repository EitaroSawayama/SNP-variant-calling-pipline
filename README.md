# SNP-variant-calling-pipline
This is for SNP calling pipline using genotyping-by-sequencing data with a reference genome.

# Create a text file, "list.txt", including sample name for loop.
This is a case of four individuals.

|sample name (This line is not needed in the real file)|
|-----------|
|ind_01|
|ind_02|
|ind_03|
|ind_04|

# Clipping Illumina adapter sequences
for file in $(<list.txt)
do
    trim_galore --paired "${file}_R1.fastq.gz" "${file}_R2.fastq.gz"
done

# Trimming low quality reads
for file in $(<list.txt)
do
    sickle pe -f "${file}_R1_val_1.fq.gz" -r "${file}_R2_val_2.fq.gz"  -o "${file}_R1.fq.gz" -p "${file}_R2.fq.gz" -s "${file}.single.fq.gz" -t sanger -q 30 -l 30 -g
done

# Create index of the reference genome

bowtie2-build -f ref_name.fasta ref_name

# Mapping

for file in $(<list.txt)
do
    bowtie2 -x ref_name -1 "${file}_R1.fq.gz" -2 "${file}_R2.fq.gz" -S "${file}.sam"
done

# Remove multiple mapped reads
for file in $(<list.txt)
do
    grep -v “XS:” "${file}.sam" > "${file}.unique.sam"
done

# Sam/bam file convertion
for file in $(<list.txt)
do
    samtools view -b "${file}.unique.sam" > "${file}.bam"
done

# Sort order of bam

for file in $(<list.txt)
do
    samtools sort "${file}.bam" -o "${file}_sorted.bam"
done

# Mpileup
bcftools mpileup -f ref_name.fasta *_sorted.bam  -I | bcftools call -m -> var.raw.bcf

bcftools view var.raw.bcf | vcfutils.pl varFilter -d20 > var.filt.vcf
