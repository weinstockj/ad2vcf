#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sysexits.h>
#include "vcfio.h"

/***************************************************************************
 *  Description:
 *      Skip over header lines in VCF input stream.
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-06  Jason Bacon Begin
 ***************************************************************************/

void    skip_header(const char *argv[], FILE *infile)

{
    char    start[7] = "xxxxxx";
    size_t  count;

    while ( ((count=fread(start, 6, 1, infile)) == 1) && 
	    (memcmp(start, "#CHROM", 6) != 0) )
	skip_rest_of_line(argv, infile);
    
    // puts(start);
    if ( count == 0 )
    {
	fprintf(stderr, "%s: Error: No #CHROM header found.\n", argv[0]);
	exit(EX_DATAERR);
    }
}


/***************************************************************************
 *  Description:
 *      Read next tab-separated field
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-06  Jason Bacon Begin
 ***************************************************************************/

int     read_field(const char *argv[], FILE *infile,
		    char buff[], size_t buff_size)

{
    size_t  c;
    char    *p;
    
    c = 0, p = buff;
    while ( (c < buff_size) && ((*p = getc(infile)) != '\t') &&
	    (*p != '\n') && (*p != EOF) )
	++c, ++p;
    
    if ( c == buff_size )
    {
	fprintf(stderr, "%s: Buffer overflow reading field.\n", argv[0]);
	fputs(buff, stderr);
	exit(EX_DATAERR);
    }
    
    if ( *p == EOF )
	return 0;
    else
    {
	*p = '\0';
	return c;
    }
}


/***************************************************************************
 *  Description:
 *      Discard next field
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-06  Jason Bacon Begin
 ***************************************************************************/

void    skip_field(const char *argv[], FILE *infile)

{
    int     ch;
    
    while ( ((ch = getc(infile)) != '\t') && (ch != '\n') && (ch != EOF) )
	;
    
    if ( ch == EOF )
    {
	fprintf(stderr, "%s: EOF encounterd unexpectedly.\n", argv[0]);
	exit(EX_DATAERR);
    }
}

/***************************************************************************
 *  Description:
 *      Extract sample IDs from input header line.
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-06  Jason Bacon Begin
 ***************************************************************************/

void    get_sample_ids(const char *argv[], FILE *infile, char *sample_ids[],
			size_t first_col, size_t last_col)

{
    size_t  c;
    char    temp_id[VCF_ID_MAX_LEN + 1];
    
    // Skip standard header tags to get to sample IDs
    for (c = 0; c < 9; ++c)
	skip_field(argv, infile);
    
    // Skip sample IDs before first_col
    for (c = 1; c < first_col; ++c)
	skip_field(argv, infile);
    
    for (; c <= last_col; ++c)
    {
	read_field(argv, infile, temp_id, VCF_ID_MAX_LEN);
	sample_ids[c - first_col] = strdup(temp_id);
	// printf("'%s'\n", sample_ids[c]);
    }
}


/***************************************************************************
 *  Description:
 *      Discard the rest of the current input line.
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-06  Jason Bacon Begin
 ***************************************************************************/

void    skip_rest_of_line(const char *argv[], FILE *infile)

{
    int     ch;
    
    while ( ((ch = getc(infile)) != EOF) && (ch != '\n') )
	;
    if ( ch == EOF )
    {
	fprintf(stderr, "%s: EOF reached unexpectedly.\n", argv[0]);
	exit(EX_DATAERR);
    }
}


/***************************************************************************
 *  Description:
 *      Read one line of a single-entry VCF file
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-08  Jason Bacon Begin
 ***************************************************************************/

int     vcf_read_call(const char *argv[],
		      FILE *vcf_stream, vcf_call_t *vcf_call)

{
    char    temp_chromosome[VCF_CHROMOSOME_NAME_MAX + 1],
	    *end;
    
    if ( read_field(argv, vcf_stream, temp_chromosome, VCF_CHROMOSOME_NAME_MAX) )
    {
	strlcpy(vcf_call->chromosome, temp_chromosome, VCF_CHROMOSOME_NAME_MAX);
	// Call position
	read_field(argv, vcf_stream, vcf_call->call_pos_str, VCF_POSITION_MAX_DIGITS);
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
	read_field(argv, vcf_stream, vcf_call->ref, VCF_REF_NAME_MAX);
	
	// Alt
	read_field(argv, vcf_stream, vcf_call->alt, VCF_ALT_NAME_MAX);

	// Qual
	skip_field(argv, vcf_stream);
	
	// Filter
	skip_field(argv, vcf_stream);
	
	// Info
	skip_field(argv, vcf_stream);
	
	// Format
	read_field(argv, vcf_stream, vcf_call->format, VCF_FORMAT_MAX);

	// Genotype
	read_field(argv, vcf_stream, vcf_call->genotype, VCF_GENOTYPE_NAME_MAX);

#if 0
	fprintf(stderr, "%s %s %s %s %s %s\n",
	    vcf_call->chromosome,
	    vcf_call->call_pos_str,
	    vcf_call->ref,
	    vcf_call->alt,
	    vcf_call->format,
	    vcf_call->genotype);
#endif
	return 1;
    }
    else
	return 0;
}
