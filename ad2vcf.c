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
		    allele,
		    nc,
		    c;
    bool            xz = false;
    size_t          vcf_pos = 0,
		    vcf_calls_read = 0,
		    alignments_read = 1,
		    previous_vcf_pos = 0,
		    previous_alignment_pos = 0;
    char            cmd[CMD_MAX + 1],
		    allele_filename[PATH_MAX + 1],
		    *vcf_chromosome,
		    previous_vcf_chromosome[VCF_CHROMOSOME_NAME_MAX + 1] = "",
		    previous_sam_rname[SAM_RNAME_MAX + 1] = "",
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

    if ( (ext = strstr(vcf_filename, ".vcf")) == NULL )
    {
	fprintf(stderr, "%s: Input filename must contain \".vcf\".\n", argv[0]);
	exit(EX_DATAERR);
    }
    
    // Insert "-ad" before ".vcf"
    *ext = '\0';
    snprintf(allele_filename, PATH_MAX, "%s-ad.%s", vcf_filename, ext+1);

    if ( xz )
    {
	snprintf(cmd, CMD_MAX, "xz -c > %s", allele_filename);
	allele_stream = popen(cmd, "w");
    }
    else
	allele_stream = fopen(allele_filename, "w");
    
    if ( allele_stream == NULL )
    {
	fprintf(stderr, "%s: Cannot open %s: %s\n", argv[0], allele_filename,
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
    
    // Check each SAM alignment against all ALT alleles.
    while ( more_alignments &&
	    ((nc = vcf_read_duplicate_calls(argv, vcf_stream,
					    &vcf_duplicate_calls)) > 0) )
    {
	// All the same
	vcf_pos = vcf_duplicate_calls.call[0].pos;
	vcf_chromosome = vcf_duplicate_calls.call[0].chromosome;
	
	// Combined VCF should be sorted by chromosome, then call position
	// Technically, this first check is only valid if it's the same
	// chromosome, but it's astronomically unlikely that an unsorted
	// VCF input will escape these checks indefinitely.  This way we get
	// by with just a single integer comparison for most iterations.
	if ( vcf_pos >= previous_vcf_pos )
	    previous_vcf_pos = vcf_pos;
	else if ( strcmp(vcf_chromosome, previous_vcf_chromosome) != 0 )
	{
	    printf("Finished chromosome, %zu calls processed.\n",
		    vcf_calls_read);
	    // Begin next chromosome, reset pos
	    strlcpy(previous_vcf_chromosome, vcf_chromosome,
		    VCF_CHROMOSOME_NAME_MAX);
	    previous_vcf_pos = vcf_pos;
	}
	else
	{
	    fprintf(stderr, "%s: VCF input must be sorted first by chromosome then position.\n", argv[0]);
	    exit(EX_DATAERR);
	}
	
	// fprintf(stderr, "VCF call pos = %zu\n", vcf_pos);
	vcf_calls_read += vcf_duplicate_calls.count;
	if ( vcf_duplicate_calls.count > 1 )
	    fprintf(stderr, "Read %zu duplicate calls at chr %s pos %zu.\n",
		    vcf_duplicate_calls.count, vcf_chromosome, vcf_pos);
	
	// Temporary
	// fprintf(allele_stream, "%s %zu ", vcf_chromosome, vcf_pos);
	
	/*
	 *  Now check all SAM alignments for the same chromosome with sequence
	 *  starting positions <= the VCF call position. Both SAM and VCF
	 *  should be sorted by chromosome first and then position, so when we
	 *  encounter a SAM alignment with a different chromosome or a position
	 *  beyond the VCF position, we're done with this VCF call.
	 *
	 *  Do integer position compare before strcmp().  It's less intuitive
	 *  to check the second sort key first, but faster since the strcmp()
	 *  is unnecessary when the integer comparison is false.
	 */
	while ( more_alignments && (sam_alignment.pos <= vcf_pos) &&
		(strcmp(vcf_chromosome, sam_alignment.rname) == 0) )
	{
	    // FIXME: Check all VCF records for the same position.
	    // FIXME: Check both streams for unsorted data
	    /*
	     *  We know at this point that the VCF call position is downstream
	     *  from the start of the SAM alignment sequence.  Is it also
	     *  upstream of the end?  If so, record the exact position and
	     *  allele.
	     */
	    if ( vcf_pos < sam_alignment.pos + sam_alignment.seq_len )
	    {
		allele = sam_alignment.seq[vcf_pos - sam_alignment.pos];
#ifdef DEBUG
		fprintf(stderr, "===\n%s pos=%zu len=%zu %s\n",
			sam_alignment.rname, sam_alignment.pos,
			sam_alignment.seq_len, sam_alignment.seq);
		fprintf(stderr, "Found allele %c (%d) for call pos %zu on %s aligned seq starting at %zu.\n",
			allele, allele, vcf_pos, vcf_chromosome,
			sam_alignment.pos);
		fprintf(stderr, "Calls: %zu  Alignments: %zu\n",
			vcf_calls_read, alignments_read);
#endif
		// Temporary
		// putc(allele, allele_stream);
		
		for (c = 0; c < vcf_duplicate_calls.count; ++c)
		{
		    if ( allele == *vcf_duplicate_calls.call[c].ref )
			++vcf_duplicate_calls.call[c].ref_count;
		    else if ( allele == *vcf_duplicate_calls.call[c].alt )
			++vcf_duplicate_calls.call[c].alt_count;
		    else
			++vcf_duplicate_calls.call[c].other_count;
		}
	    }
	    more_alignments = sam_read_alignment(argv, sam_stream, &sam_alignment);
	    ++alignments_read;
	    
	    // SAM input must be sorted by rname (chromosome), then position
	    // Technically, this first check is only valid if it's the same
	    // rname, but it's astronomically unlikely that an unsorted
	    // SAM input will escape these checks indefinitely.  This way we
	    // get by with a single integer comparison for most iterations.
	    if ( sam_alignment.pos >= previous_alignment_pos )
		previous_alignment_pos = sam_alignment.pos;
	    else if ( strcmp(sam_alignment.rname, previous_sam_rname) != 0 )
	    {
		// Begin next chromosome, reset pos
		strlcpy(previous_sam_rname, sam_alignment.rname, SAM_RNAME_MAX);
		previous_alignment_pos = sam_alignment.pos;
	    }
	    else
	    {
		fprintf(stderr, "%s: SAM input is not sorted.\n", argv[0]);
		exit(EX_DATAERR);
	    }
	}
	putc('\n', allele_stream);
	
	for (c = 0; c < vcf_duplicate_calls.count; ++c)
	{
	    fprintf(allele_stream,
		    "%s\t%zu\t.\t%s\t%s\t.\t.\t.\t%s:AD:DP\t%s:%u,%u:%u",
		    vcf_chromosome, vcf_pos,
		    vcf_duplicate_calls.call[c].ref,
		    vcf_duplicate_calls.call[c].alt,
		    vcf_duplicate_calls.call[c].format,
		    vcf_duplicate_calls.call[c].genotype,
		    vcf_duplicate_calls.call[c].ref_count,
		    vcf_duplicate_calls.call[c].alt_count,
		    vcf_duplicate_calls.call[c].ref_count +
		    vcf_duplicate_calls.call[c].alt_count);
		    // vcf_duplicate_calls.call[c].other_count);
	    putc('\n', allele_stream);
	}
	
	// Skip remaining alignments after VCF call chromosome changes
	while ( more_alignments &&
		(strcmp(vcf_chromosome, sam_alignment.rname) != 0) )
	{
	    more_alignments = sam_read_alignment(argv, sam_stream, &sam_alignment);
	    ++alignments_read;
	}
    }
    
    fprintf(stderr, "Loop terminated with more_alignments = %d, nc = %d\n",
	    more_alignments, nc);
    fprintf(stderr, "===\n%s pos=%zu len=%zu %s\n",
	    sam_alignment.rname, sam_alignment.pos,
	    sam_alignment.seq_len, sam_alignment.seq);
    fprintf(stderr, "vcf_pos = %zu, vcf_chromosome = %s\n",
	    vcf_pos, vcf_chromosome);
    
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
