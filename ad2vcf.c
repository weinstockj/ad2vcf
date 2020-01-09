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
#include <limits.h>
#include <stdbool.h>
#include <errno.h>
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
    fprintf(stderr, "Usage: %s single-sample.vcf[.xz] < file.sam\n", argv[0]);
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
    FILE            *vcf_stream,
		    *allele_stream;
    vcf_duplicate_call_t    vcf_duplicate_calls;
    sam_alignment_t sam_alignment;
    int             more_alignments,
		    allele;
    bool            xz = false;
    size_t          vcf_call_pos;
    char            cmd[CMD_MAX + 1],
		    *allele_file = "alleles.txt",
		    *vcf_call_chromosome,
		    *ext;
    const char      *vcf_filename = argv[1];
    
    xz = ((ext = strstr(vcf_filename,".xz")) != NULL) && (ext[3] == '\0');
    if ( xz )
    {
	snprintf(cmd, CMD_MAX, "unxz -c %s", vcf_filename);
	vcf_stream = popen(cmd, "r");
    }
    else
	vcf_stream = fopen(vcf_filename, "r");
    
    if ( vcf_stream == NULL )
    {
	fprintf(stderr, "%s: Cannot open %s: %s\n", argv[0], argv[1],
	    strerror(errno));
	exit(EX_NOINPUT);
    }
    
    if ( (allele_stream = fopen(allele_file, "w")) == NULL )
    {
	fprintf(stderr, "%s: Cannot open %s: %s\n", argv[0], allele_file,
	    strerror(errno));
	exit(EX_CANTCREAT);
    }
    
    /*
     *  Read in VCF fields
     */
    
    // Prime loop by reading first SAM alignment
    more_alignments = sam_read_alignment(argv, sam_stream, &sam_alignment);
    if ( ! more_alignments )
    {
	fprintf(stderr, "%s: Failed to read first SAM alignment.\n", argv[0]);
	exit(EX_DATAERR);
    }
    
    // FIXME: Be sure to handle multiple VCF calls with the same position.
    // Maybe create vcf_read_duplicate_calls() to read in all consecutive
    // calls with the same position.  Then we can check each SAM alignment
    // against all ALT alleles.
    while ( more_alignments && vcf_read_duplicate_calls(argv, vcf_stream,
				&vcf_duplicate_calls) )
    {
	// All the same
	vcf_call_pos = vcf_duplicate_calls.vcf_call[0].pos;
	vcf_call_chromosome = vcf_duplicate_calls.vcf_call[0].chromosome;
	// fprintf(stderr, "VCF call pos = %zu\n", vcf_call_pos);
	if ( vcf_duplicate_calls.count > 1 )
	    fprintf(stderr, "Read %zu duplicate calls at pos %zu.\n",
		    vcf_duplicate_calls.count, vcf_call_pos);
	
	// vcf_read_call() only gets static fields, not sample data
	// Single-sample inputs should just have one genotype in addition
	fprintf(allele_stream, "%zu ", vcf_call_pos);
	
	/*
	 *  Now check all SAM alignments for the same chromosome with sequence
	 *  starting positions <= the VCF call position. Both SAM and VCF
	 *  should be sorted by chromosome first and then position, so when we
	 *  encounter a SAM alignment with a different chromosome or a position
	 *  beyond the VCF position, we're done with this VCF call.
	 *
	 *  Do integer position compare before strcmp().  It's less intuitive
	 *  to check the second sort key first, but faster.
	 */
	while ( more_alignments && (sam_alignment.pos <= vcf_call_pos) &&
	       (strcmp(vcf_call_chromosome, sam_alignment.rname) == 0) )
	{
	    // FIXME: Check all VCF records for the same position.
	    // FIXME: Check both streams for unsorted data
	    /*
	     *  We know at this point that the VCF call position is downstream
	     *  from the start of the SAM alignment sequence.  Is it also
	     *  upstream of the end?  If so, record the exact position and
	     *  allele.
	     */
	    if ( vcf_call_pos < sam_alignment.pos + sam_alignment.seq_len )
	    {
		allele = sam_alignment.seq[vcf_call_pos - sam_alignment.pos];
#ifdef DEBUG
		fprintf(stderr, "===\n%s pos=%zu len=%zu %s\n",
			sam_alignment.rname, sam_alignment.pos,
			sam_alignment.seq_len, sam_alignment.seq);
		fprintf(stderr, "Found allele %c (%d) for call pos %zu on %s aligned seq starting at %zu.\n",
			allele, allele, vcf_call_pos, vcf_call_chromosome,
			sam_alignment.pos);
#endif
		putc(allele, allele_stream);
	    }
	    more_alignments = sam_read_alignment(argv, sam_stream, &sam_alignment);
	}
	putc('\n', allele_stream);
    }
    
    if ( xz )
	pclose(vcf_stream);
    else
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
