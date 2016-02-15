#ifndef BMF_RSQ_H
#define BMF_RSQ_H
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <zlib.h>
#include <tgmath.h>
#include "htslib/sam.h"
#include "include/sam_opts.h"
#include "include/bam.h" // for bam_get_library
#include "include/igamc_cephes.h" /// for igamc
#include "dlib/cstr_util.h"
#include "dlib/sort_util.h"
#include "dlib/bam_util.h"

#define STACK_START 128
#define SEQBUF_SIZE 300

#define seq2buf(buf, seq, len) \
	do {\
		uint64_t i_##seq;\
		for(i_##seq = 0; i_##seq < (len >> 1); ++i_##seq) {\
			buf[i_##seq] = seq_nt16_str[bam_seqi(seq, i_##seq)];\
			buf[len - i_##seq - 1] = seq_nt16_str[bam_seqi(seq, len - i_##seq - 1)];\
		}\
		if(len&1) buf[i_##seq] = seq_nt16_str[bam_seqi(seq, i_##seq)];\
		buf[len] = '\0';\
	} while(0)



typedef bam1_t *bam1_p;

typedef struct {
	size_t n, max;
	bam1_t **a;
} tmp_stack_t;

static inline void stack_insert(tmp_stack_t *stack, bam1_t *b)
{
	if (stack->n == stack->max) {
		stack->max = stack->max? stack->max<<1 : 0x10000;
		stack->a = (bam1_t**)realloc(stack->a, sizeof(bam1_t*) * stack->max);
	}
	stack->a[stack->n++] = bam_dup1(b);
}

void resize_stack(tmp_stack_t *stack, size_t n);


enum cmpkey {
	POSITION,
	UNCLIPPED
};

typedef int (*stack_fn)(bam1_t *b, bam1_t *p);

CONST static inline int same_stack_pos_se(bam1_t *b, bam1_t *p)
{
	return bmfsort_se_key(b) == bmfsort_se_key(p);
}

CONST static inline int same_stack_ucs_se(bam1_t *b, bam1_t *p)
{
	return ucs_se_sort_key(b) == ucs_se_sort_key(p);
}

CONST static inline int same_stack_pos(bam1_t *b, bam1_t *p)
{
	return (bmfsort_core_key(b) == bmfsort_core_key(p) &&
			bmfsort_mate_key(b) == bmfsort_mate_key(p));
}

CONST static inline int same_stack_ucs(bam1_t *b, bam1_t *p)
{
#if !NDEBUG
	if(!p) {
		fprintf(stderr, "Later bam record null. Abort!\n");
		exit(EXIT_FAILURE);
	}
	if(!b) {
		fprintf(stderr, "First bam record null. Abort!\n");
		exit(EXIT_FAILURE);
	}
#endif
	return (ucs_sort_core_key(b) == ucs_sort_core_key(p) &&
			ucs_sort_mate_key(b) == ucs_sort_mate_key(p));
}


typedef struct pr_settings {
	FILE *fqh;
	samFile *in;
	samFile *out;
	int cmpkey; // 0 for pos, 1 for unclipped start position
	int mmlim; // Mismatch failure threshold.
	int realign_unchanged; // Set to true to realign unchanged reads.
	int is_se; // Is single-end
	int read_hd_threshold;
	bam_hdr_t *hdr; // BAM header
	stack_fn fn;
} pr_settings_t;

int bam_rsq(int argc, char *argv[]);
void bam2ffq(bam1_t *b, FILE *fp);
void write_stack(tmp_stack_t *stack, pr_settings_t *settings);


#define READ_HD_LIMIT 10
#ifdef __cplusplus
CONST static inline int read_pass_hd(bam1_t *b, bam1_t *p, const int lim=READ_HD_LIMIT)
#else
CONST static inline int read_pass_hd(bam1_t *b, bam1_t *p, const int lim)
#endif
{
	const uint8_t *const bseq = bam_get_seq(b);
	const uint8_t *const pseq = bam_get_seq(p);
	uint8_t bc, pc;
	int hd = 0;
	for(int i = 0; i < b->core.l_qseq; ++i) {
		bc = bam_seqi(bseq, i);
		pc = bam_seqi(pseq, i);
		if(bc != pc && bc != HTS_N && pc != HTS_N)
			if(++hd > lim)
				return 0;
	}
	return hd;
}

int bam_rsq(int argc, char *argv[]);



#endif /* BMF_RSQ_H */
