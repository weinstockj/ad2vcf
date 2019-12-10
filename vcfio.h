#define ID_MAX_LEN          1024
#define CHROMOSOME_NAME_MAX 1024
#define POSITION_MAX_DIGITS 16
#define REF_NAME_MAX        64
#define ALT_NAME_MAX        64
#define GENOTYPE_NAME_MAX   64
#define FORMAT_MAX          4096

typedef struct
{
    char    chromosome[CHROMOSOME_NAME_MAX + 1],
	    call_pos_str[POSITION_MAX_DIGITS + 1],
	    ref[REF_NAME_MAX + 1],
	    alt[ALT_NAME_MAX + 1],
	    format[FORMAT_MAX + 1],
	    genotype[GENOTYPE_NAME_MAX + 1];
    size_t  call_pos;
}   vcf_call_t;

#include "vcfio-protos.h"
