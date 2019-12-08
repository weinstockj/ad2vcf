/***************************************************************************
 *  Description:
 *      Pull AD (allelic depth) info from a SAM stream and insert it
 *      into a VCF file.
 *
 *  Arguments:
 *
 *  Returns:
 *      See "man sysexits".
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-08  Jason Bacon Begin
 ***************************************************************************/

#include <stdio.h>
#include <sysexits.h>
#include <string.h>
#include <stdlib.h>
#include "ad2vcf.h"
#include "vcfio.h"
#include "samio.h"

int     main(int argc, const char *argv[])

{
    if ( argc != 2 )
	usage(argv);
    
    return ad2vcf(argv, stdin);
}


void    usage(const char *argv[])

{
    fprintf(stderr, "Usage: %s single-sample-VCF-file\n", argv[0]);
    fputs("SAM stream expected on stdin.\n", stderr);
    exit(EX_USAGE);
}


/***************************************************************************
 *  Description:
 *      1. Get list of positions from VCF file
 *      2. Get allele counts for each position from SAM stream
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-08  Jason Bacon Begin
 ***************************************************************************/

int     ad2vcf(const char *argv[], FILE *sam_stream)

{
    FILE        *vcf_stream;
    extern int  errno;
    char    chromosome[CHROMOSOME_NAME_MAX + 1],
	    position[POSITION_MAX_DIGITS + 1],
	    ref[REF_NAME_MAX + 1],
	    alt[ALT_NAME_MAX + 1],
	    format[FORMAT_MAX + 1],
	    genotype[GENOTYPE_NAME_MAX + 1];
    
    if ( (vcf_stream = fopen(argv[1], "r")) == NULL )
    {
	fprintf(stderr, "%s: Cannot open %s: %s\n", argv[0], argv[1],
	    strerror(errno));
	exit(EX_NOINPUT);
    }
    
    /*
     *  Read in VCF fields
     */
    
    // Chromosome
    while ( read_field(argv, vcf_stream, chromosome, CHROMOSOME_NAME_MAX) )
    {
	// Position
	read_field(argv, vcf_stream, position, POSITION_MAX_DIGITS);
	
	// ID
	skip_field(argv, vcf_stream);
	
	// Ref
	read_field(argv, vcf_stream, ref, REF_NAME_MAX);
	
	// Alt
	read_field(argv, vcf_stream, alt, ALT_NAME_MAX);

	// Qual
	skip_field(argv, vcf_stream);
	
	// Filter
	skip_field(argv, vcf_stream);
	
	// Info
	skip_field(argv, vcf_stream);
	
	// Format
	read_field(argv, vcf_stream, format, FORMAT_MAX);

	// Genotype
	read_field(argv, vcf_stream, genotype, GENOTYPE_NAME_MAX);

#ifdef DEBUG
	printf("%s %s %s %s %s %s\n",
	    chromosome, position, ref, alt, format, genotype);
#endif

	/*
	 *  Now search SAM input for matches to chromosome and position.
	 *  Both SAM and VCF should be sorted, so it's a game of leapfrog
	 *  through positions in the two files.
	 */
	count_alleles(argv, sam_stream, chromosome, position);
    }
    
    fclose(vcf_stream);
    return EX_OK;
}


/***************************************************************************
 *  Description:
 *      Count alleles in the SAM stream for a given chromosome and position.
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-08  Jason Bacon Begin
 ***************************************************************************/

int     count_alleles(const char *argv[], FILE *sam_stream,
		      const char chromosome[], const char position[])

{
    char    qname[SAM_QNAME_MAX + 1],
	    rname[SAM_RNAME_MAX + 1],
	    pos[SAM_POS_MAX + 1],
	    seq[SAM_SEQ_MAX + 1];
    
    // Skip header lines if present
    fprintf(stderr, "Skipping header...\n");
    while ( getc(sam_stream) == '#' )
	skip_rest_of_line(argv, sam_stream);
    fprintf(stderr, "Done skipping header...\n");
    
    while ( read_field(argv, sam_stream, qname, SAM_QNAME_MAX) )
    {
	fprintf(stderr, "QNAME = %s\n", qname);
	skip_rest_of_line(argv, sam_stream);
	
	// Flag
	skip_field(argv, sam_stream);
	
	// RNAME
	read_field(argv, sam_stream, rname, SAM_RNAME_MAX);
	
	// POS
	read_field(argv, sam_stream, pos, SAM_POS_MAX);
	
	// MAPQ
	skip_field(argv, sam_stream);
	
	// CIGAR
	skip_field(argv, sam_stream);
	
	// RNEXT
	skip_field(argv, sam_stream);
	
	// PNEXT
	skip_field(argv, sam_stream);
	
	// TLEN
	skip_field(argv, sam_stream);
	
	// SEQ
	read_field(argv, sam_stream, seq, SAM_POS_MAX);
	
	// QUAL
	skip_field(argv, sam_stream);
    }   
    return 0;
}

