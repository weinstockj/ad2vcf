#!/bin/csh -ef

#SBATCH --array=1-10
#SBATCH --cpus-per-task=3
#SBATCH --mem=1g
##SBATCH --exclude=compute-[001-008,012]

setenv PATH ../../local/bin:${PATH}
which ad2vcf

# Create array of filenames (1-based)
set samples=(`ls ../../SRR6990379/*.cram | cut -d / -f 4 | cut -d . -f 1`)
set vcf_dir=10058

cd $vcf_dir
printf "ad2vcf: processing $samples...\n"
samtools view -@ 2 --input-fmt-option required_fields=0x208 \
	../../../SRR6990379/$samples[$SLURM_ARRAY_TASK_ID].b38.irc.v1.cram \
    | ad2vcf combined.$samples[$SLURM_ARRAY_TASK_ID].vcf.xz
