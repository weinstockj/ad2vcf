#!/bin/sh -e

#SBATCH --cpus-per-task=3

time samtools view -@ 2 --input-fmt-option required_fields=0x208 \
	../SRR6990379/NWD102903.b38.irc.v1.cram \
    | ./ad2vcf ../VCF/chr01.NWD102903.vcf
exit

for vcf_file in ../VCF/*.vcf; do
    echo $vcf_file
    samtools view -@ 2 --input-fmt-option required_fields=0x208 \
	../SRR6990379/NWD102903.b38.irc.v1.cram \
	| ./ad2vcf $vcf_file
done
