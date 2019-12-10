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
#include "vcfio.h"
#include "samio.h"
#include "ad2vcf.h"

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
 *      1. Get list of call positions from VCF file
 *      2. Get allele counts for each call position from SAM stream
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-08  Jason Bacon Begin
 ***************************************************************************/

int     ad2vcf(const char *argv[], FILE *sam_stream)

{
    FILE        *vcf_stream;
    vcf_call_t  vcf_call;
    extern int  errno;
    
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
    while ( read_vcf_call(argv, vcf_stream, &vcf_call) )
    {

	/*
	 *  Now search SAM input for matches to chromosome and call position.
	 *  Both SAM and VCF should be sorted, so it's a game of leapfrog
	 *  through positions in the two files.
	 */
	count_alleles(argv, sam_stream, vcf_call.chromosome, vcf_call.call_pos);
    }
    
    fclose(vcf_stream);
    return EX_OK;
}


/***************************************************************************
 *  Description:
 *      Count alleles in the SAM stream for a given chromosome and call position.
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-08  Jason Bacon Begin
 ***************************************************************************/

int     count_alleles(const char *argv[], FILE *sam_stream,
		      const char chromosome[], size_t call_pos)

{
    char    qname[SAM_QNAME_MAX + 1],
	    rname[SAM_RNAME_MAX + 1],
	    pos_str[SAM_POS_MAX + 1],
	    seq[SAM_SEQ_MAX + 1],
	    *end;
    size_t  seq_len,
	    pos;
    
    // Skip header lines if present
    //fprintf(stderr, "Skipping header...\n");
    //while ( getc(sam_stream) == '#' )
    //    skip_rest_of_line(argv, sam_stream);
    //fprintf(stderr, "Done skipping header...\n");
    
    while ( read_field(argv, sam_stream, qname, SAM_QNAME_MAX) )
    {
	// Flag
	skip_field(argv, sam_stream);
	
	// RNAME
	read_field(argv, sam_stream, rname, SAM_RNAME_MAX);
	
	// POS
	read_field(argv, sam_stream, pos_str, SAM_POS_MAX);
	pos = strtoul(pos_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr, "%s: Invalid alignment position: %s\n",
		    argv[0], pos_str);
	    exit(EX_DATAERR);
	}
	
	
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
	seq_len = read_field(argv, sam_stream, seq, SAM_POS_MAX);
	
	// QUAL
	skip_field(argv, sam_stream);

	if ( (call_pos >= pos) && (call_pos <= pos + seq_len) )
	{
	    fprintf(stderr, "%s %zu %s %zu\n", rname, pos, seq, seq_len);
	    fprintf(stderr, "Found allele for call pos %zu on %s at %zu.\n",
		    call_pos, chromosome, pos);
	}
    }   
    return 0;
}


/***************************************************************************
 *  Description:
 *  
 *  Arguments:
 *
 *  Returns:
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-08  Jason Wayne BaconBegin
 ***************************************************************************/

int     read_vcf_call(const char *argv[],
		      FILE *vcf_stream, vcf_call_t *vcf_call)

{
    char    *end;
    
    while ( read_field(argv, vcf_stream, vcf_call->chromosome, CHROMOSOME_NAME_MAX) )
    {
	// Call position
	read_field(argv, vcf_stream, vcf_call->call_pos_str, POSITION_MAX_DIGITS);
	vcf_call->call_pos = strtoul(vcf_call->call_pos_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr, "%s: Invalid call position: %s\n",
		    argv[0], vcf_call->call_pos_str);
	    exit(EX_DATAERR);
	}
	
	// ID
	skip_field(argv, vcf_stream);
	
	// Ref
	read_field(argv, vcf_stream, vcf_call->ref, REF_NAME_MAX);
	
	// Alt
	read_field(argv, vcf_stream, vcf_call->alt, ALT_NAME_MAX);

	// Qual
	skip_field(argv, vcf_stream);
	
	// Filter
	skip_field(argv, vcf_stream);
	
	// Info
	skip_field(argv, vcf_stream);
	
	// Format
	read_field(argv, vcf_stream, vcf_call->format, FORMAT_MAX);

	// Genotype
	read_field(argv, vcf_stream, vcf_call->genotype, GENOTYPE_NAME_MAX);

#ifdef DEBUG
	printf("%s %s %s %s %s %s\n",
	    vcf_call->chromosome,
	    vcf_call->call_pos_str,
	    vcf_call->ref,
	    vcf_call->alt,
	    vcf_call->format,
	    vcf_call->genotype);
#endif
    }
    return 1;
}

