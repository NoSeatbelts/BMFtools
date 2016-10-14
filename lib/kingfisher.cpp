#include "kingfisher.h"

#include "dlib/bam_util.h"
#include "dlib/io_util.h"

namespace bmf {

void dmp_process_write(kingfisher_t *kfp, kstring_t *ks, tmpvars_t *tmp, const int is_rev)
{
    if(ks->l < kfp->readlen) ks_resize(ks, kfp->readlen);
    int i, j, diffs(kfp->length * kfp->readlen), maxindex, offset;
    double pvalues[5], cmin;
    for(i = 0; i < kfp->readlen; ++i) {
        offset = (i<<2) + i;
        cmin = 0.;
        maxindex = 4;
        for(j = 0; j < 5; ++j) {
            if(kfp->phred_sums[offset + j] > cmin) cmin = kfp->phred_sums[offset + j], maxindex = j;
            pvalues[j] = igamc_pvalues(kfp->length, LOG10_TO_CHI2(kfp->phred_sums[offset + j]));
        }
#if !NDEBUG
        //if(kfp->length > 1) fprintf(stderr, "P value before correction: %lf.\n", pvalues[maxindex]);
#endif
        for(j = 0; j < 5; ++j) if(j != maxindex) pvalues[maxindex] /= pvalues[j];
        tmp->buffers.cons_quals[i] = pvalue_to_phred(pvalues[maxindex]);
#if !NDEBUG
        //if(kfp->length > 1) fprintf(stderr, "P value after correction: %lf, (phred: %u).\n", pvalues[maxindex], tmp->buffers.cons_quals[i]);
#endif
        tmp->buffers.agrees[i] = kfp->nuc_counts[maxindex + offset];
        diffs -= maxindex == 4 ? kfp->readlen: tmp->buffers.agrees[i];
        if(tmp->buffers.cons_quals[i] < 3 || (double)tmp->buffers.agrees[i] / kfp->length < MIN_FRAC_AGREED)
            tmp->buffers.cons_quals[i] = 2, tmp->buffers.cons_seq_buffer[i] = 'N';
    }
    ksprintf(ks, "@%s ", kfp->barcode + 1);
    kfill_both(kfp->readlen, tmp->buffers.agrees, tmp->buffers.cons_quals, ks);
    tmp->buffers.cons_seq_buffer[kfp->readlen] = '\0';
    kputsnl("\tFP:i:", ks);
    kputc(kfp->pass_fail, ks);
    ksprintf(ks, "\tFM:i:%i", kfp->length);
    if(is_rev != -1) {
        ksprintf(ks, "\tRV:i:%i", is_rev ? kfp->length: 0);
        kputsnl("\tDR:i:0", ks);
    }
    LOG_DEBUG("diffs: %i.\n", (int)diffs);
    ksprintf(ks, "\tNF:f:%0.4f", (double)diffs / kfp->length);
    kputc('\n', ks);
    kputsn(tmp->buffers.cons_seq_buffer, kfp->readlen, ks);
    kputsnl("\n+\n", ks);
    for(i = 0; i < kfp->readlen; ++i) kputc(kfp->max_phreds[nuc2num(tmp->buffers.cons_seq_buffer[i]) + 5 * i], ks);
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
void zstranded_process_write(kingfisher_t *kfpf, kingfisher_t *kfpr, kstring_t *ks, tmpvars_t *tmp)
{
    const int FM(kfpf->length + kfpr->length);
    int diffs(FM * kfpf->readlen), i, j, fwidx, rvidx, offset;
    double fwmax, rvmax;
    if(ks->l < kfpf->readlen) ks_resize(ks, kfpf->readlen);
    double pval;
    if(ks->l < kfpr->readlen) ks_resize(ks, kfpr->readlen);
    LOG_DEBUG("ks %p, %s. Starting diffs: %i.\n", ks, ks->s, diffs);
    for(i = 0; i < kfpf->readlen; ++i) {
        offset = (i << 2) + i;
        fwidx = rvidx = 4; // Defaulting to N if no observations elsewhere found.
        fwmax = rvmax = 0.;
        for(j = 0; j < 5; ++j) {
            assert(offset + j < kfpf->readlen * 5);
            if(kfpf->phred_sums[offset + j] > fwmax) fwmax = kfpf->phred_sums[offset + j], fwidx = j;
            if(kfpr->phred_sums[offset + j] > rvmax) rvmax = kfpr->phred_sums[offset + j], rvidx = j;
        }
        //LOG_DEBUG("before agrees: %i diffs: %i.\n", tmp->buffers.agrees[i], diffs);
        if(fwidx == rvidx) {
            pval = igamc_pvalues(FM, LOG10_TO_CHI2(kfpf->phred_sums[offset + fwidx] + kfpr->phred_sums[offset + rvidx]));
            LOG_DEBUG("match. Confidence before subtracting: %u.\n", pvalue_to_phred(pval));
            for(j = 0; j < 4; ++j) if(j != fwidx) pval /= igamc_pvalues(FM, LOG10_TO_CHI2(kfpf->phred_sums[offset + j] + kfpr->phred_sums[offset + j]));
            LOG_DEBUG("match. Confidence after substracting: %u.\n", pvalue_to_phred(pval));
            tmp->buffers.cons_quals[i] = pvalue_to_phred(pval);
            tmp->buffers.agrees[i] = kfpf->nuc_counts[offset + fwidx] + kfpr->nuc_counts[offset + fwidx];
            tmp->buffers.cons_seq_buffer[i] = num2nuc(fwidx);
            kfpf->max_phreds[offset + fwidx] = std::max(kfpr->max_phreds[offset + rvidx], kfpf->max_phreds[offset + rvidx]);
            LOG_DEBUG("match buffer value is now %c.\n", tmp->buffers.cons_seq_buffer[i]);
            diffs -= tmp->buffers.agrees[i] + kfpf->nuc_counts[offset + 4] + kfpr->nuc_counts[offset + 4]; // Don't count masked bases.
        } else if(fwidx == 4) {
            pval = igamc_pvalues(kfpr->length, LOG10_TO_CHI2(kfpr->phred_sums[offset + rvidx]));
            LOG_DEBUG("fwn. Confidence before subtracting: %u.\n", pvalue_to_phred(pval));
            for(j = 0; j < 4; ++j) if(j != rvidx) pval /= igamc_pvalues(kfpr->length, LOG10_TO_CHI2(kfpr->phred_sums[offset + j]));
            LOG_DEBUG("fwn. Confidence after subtracting: %u.\n", pvalue_to_phred(pval));
            tmp->buffers.cons_quals[i] = pvalue_to_phred(pval);
            tmp->buffers.agrees[i] = kfpr->nuc_counts[offset + rvidx] + kfpf->nuc_counts[offset + rvidx];
            diffs -= tmp->buffers.agrees[i] + kfpf->nuc_counts[offset + 4]; // Don't count bases from masked read at position for diffs.
            tmp->buffers.cons_seq_buffer[i] = num2nuc(rvidx);
            kfpf->max_phreds[offset + fwidx] = kfpr->max_phreds[offset + rvidx];
        } else if(rvidx == 4) {
            pval = igamc_pvalues(kfpf->length, LOG10_TO_CHI2(kfpf->phred_sums[offset + fwidx]));
            LOG_DEBUG("rvn. Confidence before subtracting: %u.\n", pvalue_to_phred(pval));
            for(j = 0; j < 4; ++j) if(j != fwidx) pval /= igamc_pvalues(kfpf->length, LOG10_TO_CHI2(kfpf->phred_sums[offset + j]));
            LOG_DEBUG("rvn. Confidence before subtracting: %u.\n", pvalue_to_phred(pval));
            tmp->buffers.cons_quals[i] = pvalue_to_phred(pval);
            tmp->buffers.agrees[i] = kfpf->nuc_counts[offset + fwidx] + kfpr->nuc_counts[offset + fwidx];
            diffs -= tmp->buffers.agrees[i] + kfpr->nuc_counts[offset + 4]; // Don't count bases from masked read at position for diffs.
            tmp->buffers.cons_seq_buffer[i] = num2nuc(fwidx);
        } else {
            LOG_DEBUG("alln.\n", i);
            tmp->buffers.cons_quals[i] = tmp->buffers.agrees[i] = 0;
            tmp->buffers.cons_seq_buffer[i] = 'N';
            kfpf->max_phreds[offset + 4] = '#';
            diffs -= FM; // Don't count bases from masked read at position for diffs.
        }
        //LOG_DEBUG("after agrees: %i diffs: %i.\n", tmp->buffers.agrees[i], diffs);
    }
    ksprintf(ks, "@%s ", kfpf->barcode + 1);
    // Add read name
    LOG_DEBUG("About to kfill\n");
    kfill_both(kfpf->readlen, tmp->buffers.agrees, tmp->buffers.cons_quals, ks);
    tmp->buffers.cons_seq_buffer[kfpf->readlen] = '\0';
    ksprintf(ks, "\tFP:i:%c\tFM:i:%i\tRV:i:%i\tNF:f:%f\tDR:i:%i\n%s\n+\n", kfpf->pass_fail,
             FM, kfpr->length, (double) diffs / FM, kfpf->length && kfpr->length,
             tmp->buffers.cons_seq_buffer);
    LOG_DEBUG("Add max phreds\n");
    for(i = 0; i < kfpf->readlen; ++i)
        kputc(kfpf->max_phreds[nuc2num(tmp->buffers.cons_seq_buffer[i]) + (i << 2) + i], ks);
    kputc('\n', ks);
    //const int ND = get_num_differ
    return;
}

} /* namespace bmf */
