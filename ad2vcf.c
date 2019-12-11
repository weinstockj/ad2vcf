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
#include "tsvio.h"
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
    FILE            *vcf_stream;
    vcf_call_t      vcf_call;
    sam_alignment_t sam_alignment;
    extern int      errno;
    int             more_alignments;
    char            genotype[VCF_GENOTYPE_NAME_MAX + 1];
    
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
    more_alignments = sam_read_alignment(argv, sam_stream, &sam_alignment);
    if ( ! more_alignments )
    {
	fprintf(stderr, "%s: Failed to read first SAM alignment.\n", argv[0]);
	exit(EX_DATAERR);
    }
    while ( more_alignments && vcf_read_call(argv, vcf_stream, &vcf_call) )
    {
	// Ignore sample data
	tsv_read_field(argv, vcf_stream, genotype, VCF_GENOTYPE_NAME_MAX);
	
	/*
	 *  Now search SAM input for matches to chromosome and call position.
	 *  Both SAM and VCF should be sorted, so it's a game of leapfrog
	 *  through positions in the two files.
	 */
	while ( more_alignments && (sam_alignment.pos <= vcf_call.pos) &&
	       (strcmp(vcf_call.chromosome, sam_alignment.rname) == 0) )
	{
	    // fprintf(stderr, "%s %s\n", vcf_call.chromosome, sam_alignment.rname);
	    if ( vcf_call.pos <= sam_alignment.pos + sam_alignment.seq_len )
	    {
		fprintf(stderr, "===\n%s pos=%zu len=%zu %s\n",
			sam_alignment.rname, sam_alignment.pos,
			sam_alignment.seq_len, sam_alignment.seq);
		fprintf(stderr, "Found allele %c for call pos %zu on %s aligned seq starting at %zu.\n",
			sam_alignment.seq[vcf_call.pos - sam_alignment.pos],
			vcf_call.pos, vcf_call.chromosome, sam_alignment.pos);
	    }
	    more_alignments = sam_read_alignment(argv, sam_stream, &sam_alignment);
	}
    }
    
    fclose(vcf_stream);
    return EX_OK;
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
    
    if ( tsv_read_field(argv, sam_stream, sam_alignment->qname, SAM_QNAME_MAX) )
    {
	// Flag
	tsv_skip_field(argv, sam_stream);
	
	// RNAME
	tsv_read_field(argv, sam_stream, sam_alignment->rname, SAM_RNAME_MAX);
	
	// POS
	tsv_read_field(argv, sam_stream, pos_str, SAM_POS_MAX_DIGITS);
	sam_alignment->pos = strtoul(pos_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr, "%s: Invalid alignment position: %s\n",
		    argv[0], pos_str);
	    exit(EX_DATAERR);
	}
	
	// MAPQ
	tsv_skip_field(argv, sam_stream);
	
	// CIGAR
	tsv_skip_field(argv, sam_stream);
	
	// RNEXT
	tsv_skip_field(argv, sam_stream);
	
	// PNEXT
	tsv_skip_field(argv, sam_stream);
	
	// TLEN
	tsv_skip_field(argv, sam_stream);
	
	// SEQ
	sam_alignment->seq_len = tsv_read_field(argv, sam_stream, sam_alignment->seq, SAM_SEQ_MAX);
	
	// QUAL
	tsv_skip_field(argv, sam_stream);
	return 1;
    }
    else
	return 0;
}
