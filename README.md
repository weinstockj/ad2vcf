# ad2vcf

ad2vdf extracts allelic depth info from a SAM stream and adds it to a
corresponding single-sample VCF file.

SAM input is read via stdin and the VCF input file is taken as a command-line
argument.  This allows expensive BAM/CRAM decoding to occur in-parallel using
a pipe:

...
samtools view -@ 2 --input-fmt-option required_fields=0x208 \
        ../SRR6990379/NWD102903.b38.irc.v1.cram \
        | ./ad2vcf file.vcf
...
