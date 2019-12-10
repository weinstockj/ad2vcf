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
    
    // Leapfrog positions in SAM and VCF streams
    while ( vcf_read_call(argv, vcf_stream, &vcf_call) )
    {

	/*
	 *  Now search SAM input for matches to chromosome and call position.
	 *  Both SAM and VCF should be sorted, so it's a game of leapfrog
	 *  through positions in the two files.
	 */
	sam_count_alleles(argv, sam_stream, vcf_call.chromosome, vcf_call.call_pos);
    }
    
    fclose(vcf_stream);
    return EX_OK;
}


/***************************************************************************
 *  Description:
 *      Count alleles in the SAM stream for a given chromosome and VCF
 *      call position.
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-08  Jason Bacon Begin
 ***************************************************************************/

int     sam_count_alleles(const char *argv[], FILE *sam_stream,
		      const char vcf_chromosome[], size_t call_pos)

{
    sam_alignment_t sam_alignment;
    
    while ( sam_read_alignment(argv, sam_stream, &sam_alignment) )
    {

	if ( (call_pos >= sam_alignment.pos) &&
	     (call_pos <= sam_alignment.pos + sam_alignment.seq_len) )
	{
	    fprintf(stderr, "===\n%s pos=%zu len=%zu %s\n",
		    sam_alignment.rname, sam_alignment.pos,
		    sam_alignment.seq_len, sam_alignment.seq);
	    fprintf(stderr, "Found allele %c for call pos %zu on %s aligned seq starting at %zu.\n",
		    sam_alignment.seq[call_pos - sam_alignment.pos], call_pos,
		    vcf_chromosome, sam_alignment.pos);
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
 *  2019-12-09  Jason Bacon Begin
 ***************************************************************************/

int     sam_read_alignment(const char *argv[],
			   FILE *sam_stream, sam_alignment_t *sam_alignment)

{
    char    pos_str[SAM_POS_MAX_DIGITS + 1],
	    *end;
    
    if ( read_field(argv, sam_stream, sam_alignment->qname, SAM_QNAME_MAX) )
    {
	// Flag
	skip_field(argv, sam_stream);
	
	// RNAME
	read_field(argv, sam_stream, sam_alignment->rname, SAM_RNAME_MAX);
	
	// POS
	read_field(argv, sam_stream, pos_str, SAM_POS_MAX_DIGITS);
	sam_alignment->pos = strtoul(pos_str, &end, 10);
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
	sam_alignment->seq_len = read_field(argv, sam_stream, sam_alignment->seq, SAM_SEQ_MAX);
	
	// QUAL
	skip_field(argv, sam_stream);
	return 1;
    }
    else
	return 0;
}
