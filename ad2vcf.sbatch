#!/bin/sh -e

#SBATCH --ntasks=3

export PATH=../../local/bin:${PATH}

sample=NWD102903
samtools view -@ 2 --input-fmt-option required_fields=0x208 \
	../../SRR6990379/$sample.b38.irc.v1.cram \
    | ad2vcf combined.$sample.vcf.xz
