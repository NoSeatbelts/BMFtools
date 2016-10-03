#include "kingfisher.h"

#include "dlib/bam_util.h"
#include "dlib/io_util.h"

namespace bmf {

void dmp_process_write(kingfisher_t *kfp, kstring_t *ks, tmpbuffers_t *bufs, int is_rev)
{
    int i, j, diffs(kfp->length * kfp->readlen), maxindex;
    double pvalues[5], cmin;
    for(i = 0; i < kfp->readlen; ++i) {
        const int offset((i<<2) + i);
        cmin = 0.;
        maxindex = 4;
        for(j = 0; j < 5; ++j) {
            if(kfp->phred_sums[offset + j] > cmin) cmin = kfp->phred_sums[offset + j], maxindex = j;
            pvalues[i] = igamc_pvalues(kfp->length, LOG10_TO_CHI2(kfp->phred_sums[offset + j]));
        }
        bufs->cons_quals[i] = pvalue_to_phred(pvalues[maxindex]);
        for(j = 0; j < 5; ++j) if(j != maxindex) pvalues[maxindex] /= pvalues[j];
        bufs->agrees[i] = kfp->nuc_counts[maxindex + offset];
        diffs -= bufs->agrees[i];
        if(maxindex != 4) diffs -= kfp->nuc_counts[offset + maxindex];
        if(bufs->cons_quals[i] < 3 || (double)bufs->agrees[i] / kfp->length < MIN_FRAC_AGREED)
            bufs->cons_quals[i] = 2, bufs->cons_seq_buffer[i] = 'N';
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
    int diffs(FM * kfpf->readlen), i, j, fwidx, rvidx, offset;
    double fwmax, rvmax;
    for(i = 0; i < kfpf->readlen; ++i) {
        offset = (i << 2) + i;
        double pvalues[5];
        fwidx = rvidx = 4; // Defaulting to N if no observations elsewhere found.
        fwmax = rvmax = 0.;
        for(j = 0; j < 4; ++j) {
            if(kfpf->phred_sums[offset + j] > fwmax) fwmax = kfpf->phred_sums[offset + j], fwidx = j;
            if(kfpr->phred_sums[offset + j] > rvmax) rvmax = kfpr->phred_sums[offset + j], rvidx = j;
        }
        if(fwidx == rvidx) {
            pvalues[i] = igamc_pvalues(kfpr->length + kfpf->length, LOG10_TO_CHI2(kfpf->phred_sums[offset + fwidx] + kfpr->phred_sums[offset + rvidx]));
            for(j = 0; j < 4; ++j) if(j != fwidx) pvalues[i] /= igamc_pvalues(kfpr->length + kfpf->length, kfpf->phred_sums[offset + j] + kfpr->phred_sums[offset + j]);
            bufs->cons_quals[i] = pvalue_to_phred(pvalues[i]);
            bufs->agrees[i] = kfpf->phred_sums[offset + fwidx] + kfpr->phred_sums[offset + rvidx];
            bufs->cons_seq_buffer[i] = kfpf->max_phreds[offset + rvidx] < kfpr->max_phreds[offset + rvidx] ? kfpr->max_phreds[offset + rvidx]: kfpf->max_phreds[offset + rvidx];
        } else if(fwidx == 4) {
            pvalues[i] = igamc_pvalues(kfpr->length, LOG10_TO_CHI2(kfpr->phred_sums[offset + rvidx]));
            for(j = 0; j < 4; ++j) if(j != rvidx) pvalues[i] /= igamc_pvalues(kfpr->length, kfpr->phred_sums[offset + j]);
            bufs->cons_quals[i] = pvalue_to_phred(pvalues[i]);
            bufs->agrees[i] = kfpr->phred_sums[offset + rvidx];
            diffs -= kfpf->length; // Don't count bases from masked read at position for diffs.
            bufs->cons_seq_buffer[i] = kfpr->max_phreds[offset + rvidx];
        } else if(rvidx == 4) {
            pvalues[i] = igamc_pvalues(kfpf->length, LOG10_TO_CHI2(kfpf->phred_sums[offset + fwidx]));
            for(j = 0; j < 4; ++j) if(j != fwidx) pvalues[i] /= igamc_pvalues(kfpf->length, kfpf->phred_sums[offset + j]);
            bufs->cons_quals[i] = pvalue_to_phred(pvalues[i]);
            bufs->agrees[i] = kfpf->nuc_counts[offset + fwidx];
            diffs -= kfpr->length; // Don't count bases from masked read at position for diffs.
            bufs->cons_seq_buffer[i] = kfpf->max_phreds[offset + fwidx];
        } else {
            pvalues[i] = 2;
            bufs->cons_quals[i] = pvalue_to_phred(pvalues[i]);
            bufs->agrees[i] = 0;
            diffs -= kfpr->length + kfpf->length; // Don't count bases from masked read at position for diffs.
            bufs->cons_seq_buffer[i] = 'N';
        }
        diffs -= bufs->agrees[i];
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
