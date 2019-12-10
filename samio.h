#define SAM_QNAME_MAX       4096
#define SAM_FLAG_MAX        4096
#define SAM_RNAME_MAX       4096
#define SAM_POS_MAX_DIGITS  4096
#define SAM_CIGAR_MAX       4096
#define SAM_SEQ_MAX         1024*1024

typedef struct
{
    char    qname[SAM_QNAME_MAX + 1],
	    rname[SAM_RNAME_MAX + 1],
	    seq[SAM_SEQ_MAX + 1];
    size_t  pos;
    size_t  seq_len;
}   sam_alignment_t;
