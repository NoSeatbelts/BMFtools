#ifndef KINGFISHER_H
#define KINGFISHER_H
#include <assert.h>
#include <cmath>
#include <zlib.h>
#include "htslib/khash.h"
#include "htslib/kseq.h"
#include "htslib/kstring.h"
#include "dlib/cstr_util.h"
#include "include/igamc_cephes.h"
#include "lib/splitter.h"

#ifndef MAX_BARCODE_LENGTH
#define MAX_BARCODE_LENGTH 31
#endif
#ifndef MAX_PV
#    define MAX_PV 3117 // Maximum seen with doubles
#endif
#define HASH_DMP_OFFSET 14
#define FP_OFFSET 9

#ifndef KSEQ_DEC_GZ
#define KSEQ_DEC_GZ
KSEQ_INIT(gzFile, gzread)
#endif
namespace bmf {

const double MIN_FRAC_AGREED = 0.5; // Minimum fraction of bases agreed in a family to not "N" the base.


struct tmpvars_t {
    int blen;
    int readlen;
    char key[MAX_BARCODE_LENGTH + 1];
    struct tmpbuffers_t {
        char name_buffer[120];
        char PVBuffer[1000];
        char FABuffer[1000];
        char cons_seq_buffer[SEQBUF_SIZE];
        uint32_t cons_quals[SEQBUF_SIZE];
        uint16_t agrees[SEQBUF_SIZE];
    } buffers;
    tmpvars_t() {memset(this, 0, sizeof *this);}
    tmpvars_t(int blen_, int readlen_): blen(blen_), readlen(readlen_) {
        buffers.name_buffer[0] = '@';
        buffers.name_buffer[blen] = '\0';
        buffers.cons_seq_buffer[readlen] = '\0';
    }
};


struct kingfisher_t {
    uint16_t *nuc_counts; // Count of nucleotides of this form
    uint32_t *phred_sums; // Sums of -10log10(p-value)
    char *max_phreds; // Maximum phred score observed at position. Use this as the final sequence for the quality to maintain compatibility with GATK and other tools.
    char barcode[MAX_BARCODE_LENGTH + 1];
    int length:16; // Number of reads in family
    int readlen:12; // Length of reads
    int pass_fail:8;
};

static inline void destroy_kf(kingfisher_t *kfp)
{
    free(kfp->nuc_counts);
    free(kfp->phred_sums);
    free(kfp->max_phreds);
    free(kfp);
}


void zstranded_process_write(kingfisher_t *kfpf, kingfisher_t *kfpr, kstring_t *ks, tmpvars_t *bufs);
void dmp_process_write(kingfisher_t *kfp, kstring_t *ks, tmpvars_t *bufs, int is_rev);
int kf_hamming(kingfisher_t *kf1, kingfisher_t *kf2);

static inline void kfill_both(int readlen, uint16_t *agrees, uint32_t *quals, kstring_t *ks)
{
    int i;
    kputsnl("FA:B:I", ks);
    for(i = 0; i < readlen; ++i) kputc(',', ks), kputw(agrees[i], ks);
    kputsnl("\tPV:B:I", ks);
    for(i = 0; i < readlen; ++i) kputc(',', ks), kputw(quals[i], ks);
}

static inline void pb_pos(kingfisher_t *kfp, kseq_t *seq, int i) {
    const uint32_t posdata(nuc2num(seq->seq.s[i]) + i * 5);
    ++kfp->nuc_counts[posdata];
    kfp->phred_sums[posdata] += seq->qual.s[i] - 33;
    if(seq->qual.s[i] > kfp->max_phreds[posdata]) kfp->max_phreds[posdata] = seq->qual.s[i];
}

static inline void pushback_inmem(kingfisher_t *kfp, kseq_t *seq, int offset, int pass) {
    if(kfp->length++) {
        if(kfp->readlen + offset != (int64_t)seq->seq.l) {
            if(pass) return; // Don't bother, it's an error.
            offset = seq->seq.l - kfp->readlen;
        }
    } else kfp->pass_fail = pass + '0';
    uint32_t i;
    // Reuse pass instead of allocating another variable.
    for(i = offset; i < seq->seq.l; ++i) {
        pass = nuc2num(seq->seq.s[i]) + (i - offset) * 5;
        assert(pass < kfp->readlen * 5);
        ++kfp->nuc_counts[pass];
        kfp->phred_sums[pass] += seq->qual.s[i] - 33;
        if(seq->qual.s[i] > kfp->max_phreds[pass])
            kfp->max_phreds[pass] = seq->qual.s[i];
    }
}

static inline void pushback_kseq(kingfisher_t *kfp, kseq_t *seq, int blen)
{
    if(!kfp->length++) { // Increment while checking
        kfp->pass_fail = seq->comment.s[FP_OFFSET];
        memcpy(kfp->barcode, seq->comment.s + HASH_DMP_OFFSET, blen);
        kfp->barcode[blen] = '\0';
    }
    // Reuse blen int to avoid allocating another integer.
    for(blen = 0; blen < kfp->readlen; ++blen) pb_pos(kfp, seq, blen);
}


/*
 * @func arr_max_u32
 * :param: arr [uint32_t *] 2-d array of values. 5 * index + basecall is the index to use.
 * :param: index [int] Base in read to find the maximum value for.
 * :returns: [int] the nucleotide number for the maximum value at this index in the read.
 */
CONST static inline int arr_max_u32(uint32_t *arr, int index)
{
    arr += index * 5;
    return (arr[0] > arr[1]) ? ((arr[0] > arr[2]) ? ((arr[0] > arr[3]) ? (arr[0] > arr[4] ? 0: 4)
                                                                       : (arr[3] > arr[4] ? 3: 4))
                                                  : (arr[2] > arr[3])  ? (arr[2] > arr[4] ? 2: 4)
                                                                       : (arr[3] > arr[4] ? 3: 4))
                             : ((arr[1] > arr[2]) ? ((arr[1] > arr[3]) ? (arr[1] > arr[4] ? 1: 4)
                                                                       : (arr[3] > arr[4] ? 3: 4))
                                                  : ((arr[2] > arr[3]) ? (arr[2] > arr[4] ? 2: 4)
                                                                       : (arr[3] > arr[4] ? 3: 4)));

}


CONST static inline int kfp_argmax(kingfisher_t *kfp, int index)
{
    return arr_max_u32(kfp->phred_sums, index);
}

std::vector<double> get_igamc_threshold(int family_size, int max_phred=MAX_PV, double delta=0.002);
std::vector<std::vector<double>> get_igamc_thresholds(int max_family_size, int max_phred=MAX_PV, double delta=0.002);

} /* namespace bmf */


#endif /*KINGFISHER_H*/
