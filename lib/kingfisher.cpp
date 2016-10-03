#include "kingfisher.h"

#include "dlib/bam_util.h"
#include "dlib/io_util.h"

namespace bmf {

#define dmp_pos(kfp, bufs, argmaxret, i, index, diffcount)\
    do {\
        bufs->cons_quals[i] = pvalue_to_phred(igamc_pvalues(kfp->length, LOG10_TO_CHI2((kfp->phred_sums[index]))));\
        bufs->agrees[i] = kfp->nuc_counts[index];\
        diffcount -= bufs->agrees[i];\
        if(argmaxret != 4) diffcount -= kfp->nuc_counts[i * 5 + 4]; /*(Skip Ns in counting diffs) */\
        if(bufs->cons_quals[i] > 2 && (double)bufs->agrees[i] / kfp->length > MIN_FRAC_AGREED) {\
            bufs->cons_seq_buffer[i] = num2nuc(argmaxret);\
        } else {\
            bufs->cons_quals[i] = 2;\
            bufs->cons_seq_buffer[i] = 'N';\
        }\
    } while(0)

void dmp_process_write(kingfisher_t *kfp, kstring_t *ks, tmpbuffers_t *bufs, int is_rev)
{
    int i, j, diffs(kfp->length * kfp->readlen), maxindex;
    double pvalues[5], cmax;
    for(i = 0; i < kfp->readlen; ++i) {
        const int offset((i<<2) + i);
        cmax = 0.;
        maxindex = 0;
        for(j = 0; j < 5; ++j) {
            if(kfp->phred_sums[offset + j] > cmax) cmax = kfp->phred_sums[offset + j], maxindex = j;
            pvalues[i] = igamc_pvalues(kfp->length, LOG10_TO_CHI2(kfp->phred_sums[offset + j]));
        }
        bufs->cons_quals[i] = pvalue_to_phred(pvalues[maxindex]);
        for(j = 0; j < 5; ++j) if(j != maxindex) pvalues[maxindex] /= pvalues[j];
        bufs->agrees[i] = kfp->nuc_counts[maxindex + offset];
        diffs -= bufs->agrees[i];
        if(maxindex != 4) diffs -= kfp->nuc_counts[offset + maxindex];
        if(bufs->cons_quals[i] < 3 || (double)bufs->agrees[i] / kfp->length < MIN_FRAC_AGREED)
            bufs->cons_quals[i] = 2, bufs->cons_seq_buffer[i] = 'N';
        dmp_pos(kfp, bufs, argmaxret, i, index, diffs);
    }
    ksprintf(ks, "@%s ", kfp->barcode + 1);
    kfill_both(kfp->readlen, bufs->agrees, bufs->cons_quals, ks);
    bufs->cons_seq_buffer[kfp->readlen] = '\0';
    kputsnl("\tFP:i:", ks);
    kputc(kfp->pass_fail, ks);
    ksprintf(ks, "\tFM:i:%i", kfp->length);
    if(is_rev != -1) {
        ksprintf(ks, "\tRV:i:%i", is_rev ? kfp->length: 0);
        kputsnl("\tDR:i:0", ks);
    }
    ksprintf(ks, "\tNF:f:%0.4f", (double) diffs / kfp->length);
    kputc('\n', ks);
    kputsn(bufs->cons_seq_buffer, kfp->readlen, ks);
    kputsnl("\n+\n", ks);
    for(i = 0; i < kfp->readlen; ++i) kputc(kfp->max_phreds[nuc2num(bufs->cons_seq_buffer[i]) + 5 * i], ks);
    kputc('\n', ks);
}

std::vector<double> get_igamc_threshold(int family_size, int max_phred, double delta) {
    std::vector<double> ret;
    double query(delta);
    int last_pv(-1), current_pv(-1);
    while((last_pv = -10 * log10(igamc_pvalues(family_size, query))) < 0)
        query += delta;
    assert(last_pv == 0);
    current_pv = 0;
    ret.push_back(query);
    while(current_pv <= MAX_PV) {
        query += delta;
        current_pv = -10 * log10(igamc_pvalues(family_size, query));
        if(current_pv != last_pv) {
            ret.push_back(query);
            last_pv = current_pv;
            assert((unsigned)current_pv + 1 == ret.size());
        }
    }
    return ret;
}

std::vector<std::vector<double>> get_igamc_thresholds(size_t max_family_size, int max_phred, double delta) {
    std::vector<std::vector<double>> ret;
    while(ret.size() < max_family_size)
        ret.emplace_back(get_igamc_threshold(ret.size() + 1, max_phred, delta));
    return ret;
}

int kf_hamming(kingfisher_t *kf1, kingfisher_t *kf2) {
    int ret(0);
    int argmaxret1, argmaxret2;
    for(int i(0); i < kf1->readlen; ++i) {
        argmaxret1 = kfp_argmax(kf1, i);
        argmaxret2 = kfp_argmax(kf2, i);
        if(argmaxret1 != argmaxret2)
            if(argmaxret1 != 4)
                if(argmaxret2 != 4)
                    ++ret;
    }
    return ret;
}

// kfp forward, kfp reverse
// Note: You print kfpf->barcode + 1 because that skips the F/R/Z char.
void zstranded_process_write(kingfisher_t *kfpf, kingfisher_t *kfpr, kstring_t *ks, tmpbuffers_t *bufs)
{
    const int FM (kfpf->length + kfpr->length);
    int diffs(FM * kfpf->readlen), index, i;
    for(i = 0; i < kfpf->readlen; ++i) {
        const int argmaxretf(kfp_argmax(kfpf, i)); // Forward consensus nucleotide
        const int argmaxretr(kfp_argmax(kfpr, i)); // Reverse consensus nucleotide
        if(argmaxretf == argmaxretr) { // Both strands supported the same base call.
            index = i * 5 + argmaxretf;
            kfpf->phred_sums[index] += kfpr->phred_sums[index];
            kfpf->nuc_counts[index] += kfpr->nuc_counts[index];
            dmp_pos(kfpf, bufs, argmaxretf, i, index, diffs);
            if(kfpr->max_phreds[index] > kfpf->max_phreds[index]) kfpf->max_phreds[index] = kfpr->max_phreds[index];
        } else if(argmaxretf == 4) { // Forward is N'd and reverse is not. Reverse call is probably right.
            index = i * 5 + argmaxretr;
            kfpf->phred_sums[index] += kfpr->phred_sums[index];
            kfpf->nuc_counts[index] += kfpr->nuc_counts[index];
            dmp_pos(kfpf, bufs, argmaxretr, i, index, diffs);
            kfpf->max_phreds[index] = kfpr->max_phreds[index];
        } else if(argmaxretr == 4) { // Forward is N'd and reverse is not. Reverse call is probably right.
            index = i * 5 + argmaxretf;
            kfpf->phred_sums[index] += kfpr->phred_sums[index];
            kfpf->nuc_counts[index] += kfpr->nuc_counts[index];
            dmp_pos(kfpf, bufs, argmaxretf, i, index, diffs);
            // Don't update max_phreds, since the max phred is already here.
        } else bufs->cons_quals[i] = 0, bufs->agrees[i] = 0, bufs->cons_seq_buffer[i] = 'N';
    }
    ksprintf(ks, "@%s ", kfpf->barcode + 1);
    // Add read name
    kfill_both(kfpf->readlen, bufs->agrees, bufs->cons_quals, ks);
    bufs->cons_seq_buffer[kfpf->readlen] = '\0';
    ksprintf(ks, "\tFP:i:%c\tFM:i:%i\tRV:i:%i\tNF:f:%f\tDR:i:%i\n%s\n+\n", kfpf->pass_fail,
             FM, kfpr->length, (double) diffs / FM, kfpf->length && kfpr->length,
             bufs->cons_seq_buffer);
    for(i = 0; i < kfpf->readlen; ++i)
        kputc(kfpf->max_phreds[nuc2num(bufs->cons_seq_buffer[i]) + 5 * i], ks);
    kputc('\n', ks);
    //const int ND = get_num_differ
    return;
}

} /* namespace bmf */
