/* ad2vcf.c */
int main(int argc, const char *argv[]);
void usage(const char *argv[]);
int ad2vcf(const char *argv[], FILE *sam_stream);
int count_alleles(const char *argv[], FILE *sam_stream, const char chromosome[], const char position[]);
