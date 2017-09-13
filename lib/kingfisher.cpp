#include "kingfisher.h"

#include "dlib/bam_util.h"
#include "dlib/io_util.h"
#include <limits>

namespace bmf {

void dmp_process_write(kingfisher_t *kfp, kstring_t *ks, tmpbuffers_t *bufs, int is_rev)
{
    int i, j, diffs(kfp->length * kfp->readlen), minind(4);
    int pvals[5];
    double minv(std::numeric_limits<double>::max());
    for(i = 0; i < kfp->readlen; ++i) {
        for(minind = 4, j = i * 5; j < (i + 1) * 5; ++j) {
            pvals[j - i * 5] = kfp->phred_sums[j] ? pvalue_to_phred(igamc_pvalues(kfp->length, LOG10_TO_CHI2((kfp->phred_sums[j]))))
                                                  : 0;
            if(pvals[j - i * 5] < minv) minv = pvals[j - i * 5], minind = j - i * 5;
        }
        for(j = 0; j < 5; ++j) if(j != minind) pvals[minind] -= pvals[j];
        if(pvals[minind] < 0) pvals[minind] = 0;
        bufs->cons_quals[i] = pvals[minind];
        bufs->agrees[i]     = kfp->nuc_counts[minind + i * 5];
        diffs          -= bufs->agrees[i];
        if(minind != 4) diffs -= kfp->nuc_counts[i * 5 + 4];
        if(bufs->cons_quals[i] > 2 && (double)bufs->agrees[i] / kfp->length > MIN_FRAC_AGREED)
            bufs->cons_seq_buffer[i] = num2nuc(minind);
        else bufs->cons_quals[i] = 2, bufs->cons_seq_buffer[i] = 'N';
    }
    kputc('@', ks); kputs(kfp->barcode + 1, ks); kputc(' ', ks); 
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

// kfp forward, kfp reverse
// Note: You print kfpf->barcode + 1 because that skips the F/R/Z char.
void zstranded_process_write(kingfisher_t *kfpf, kingfisher_t *kfpr, kstring_t *ks, tmpbuffers_t *bufs)
{
    const int FM(kfpf->length + kfpr->length);
    int diffs(FM * kfpf->readlen), index, i, j,
        fminv(std::numeric_limits<int>::max()), rminv(fminv),
        fminind, rminind, fpsums[5], rpsums[5];
    for(i = 0; i < kfpf->readlen; ++i) {
        for(rminind = fminind = 4, j = i * 5; j < (i + 1) * 5; ++j) {
            if((fpsums[j - i * 5] = kfpr->phred_sums[j]) < fminv) 
                fminv = fpsums[j - i * 5], fminind = j - i * 5;
            if((rpsums[j - i * 5] = kfpr->phred_sums[j]) < rminv) 
                rminv = rpsums[j - i * 5], rminind = j - i * 5;
        }
#if 0
        if(psums[minind] < 0) psums[minind] = 0;
        bufs->cons_quals[i] = psums[minind];
        bufs->agrees[i]     = kfp->nuc_counts[minind + i * 5];
        diffs          -= bufs->agrees[i];
        if(minind != 4) diffs -= kfp->nuc_counts[i * 5 + 4];
        if(bufs->cons_quals[i] > 2 && (double)bufs->agrees[i] / kfp->length > MIN_FRAC_AGREED)
            bufs->cons_seq_buffer[i] = num2nuc(minind);
        else bufs->cons_quals[i] = 2, bufs->cons_seq_buffer[i] = 'N';
#endif
        if(fminind == rminind) { // Both strands supported the same base call.
            index = i * 5 + fminind;
            for(j = 0; j < 5; ++j) {
                fpsums[j] = pvalue_to_phred(igamc_pvalues(kfpf->length, LOG10_TO_CHI2((kfpf->phred_sums[j] + kfpr->phred_sums[j]))));
            }
            for(j = 0; j < 5; ++j) if(j != fminind) fpsums[fminind] -= fpsums[j];
            if(fpsums[fminind] < 0) fpsums[fminind] = 0;
            kfpf->nuc_counts[index] += kfpr->nuc_counts[index];
            bufs->cons_quals[i] = fpsums[fminind];
            bufs->agrees[i]     = kfpf->nuc_counts[index] + kfpr->nuc_counts[index];
            diffs              -= bufs->agrees[i];
            if(fminind != 4) diffs -= kfpf->nuc_counts[i * 5 + 4] + kfpr->nuc_counts[i * 5 + 4];
            if(kfpr->max_phreds[index] > kfpf->max_phreds[index]) kfpf->max_phreds[index] = kfpr->max_phreds[index];
            if(bufs->cons_quals[i] > 2 && (double)bufs->agrees[i] / kfpf->length > MIN_FRAC_AGREED)
                bufs->cons_seq_buffer[i] = num2nuc(fminind);
            else bufs->cons_quals[i] = 2, bufs->cons_seq_buffer[i] = 'N';
        } else if(fminind == 4) { // Forward is N'd and reverse is not. Reverse call is probably right.
            index = i * 5 + rminind;
            kfpf->phred_sums[index] += kfpr->phred_sums[index];
            kfpf->nuc_counts[index] += kfpr->nuc_counts[index];
            kfpf->max_phreds[index] = kfpr->max_phreds[index];
            diffs -= kfpr->nuc_counts[index];
            for(j = 0; j < 5; ++j) {
                rpsums[j] = pvalue_to_phred(igamc_pvalues(kfpr->length, LOG10_TO_CHI2((kfpf->phred_sums[j] + kfpr->phred_sums[j]))));
            }
            for(j = 0; j < 5; ++j) if(j != rminind) rpsums[rminind] -= rpsums[j];
            if(rpsums[rminind] < 0) rpsums[rminind] = 0;
            bufs->cons_quals[i] = rpsums[rminind];
            bufs->agrees[i] = kfpf->nuc_counts[rminind + i * 5] + kfpr->nuc_counts[rminind + i * 5];
            diffs -= bufs->agrees[i];
            diffs -= kfpf->nuc_counts[4 + i * 5] + kfpr->nuc_counts[4 + i * 5];
            if(bufs->cons_quals[i] > 2 && (double)bufs->agrees[i] / kfpr->length > MIN_FRAC_AGREED)
                bufs->cons_seq_buffer[i] = num2nuc(fminind);
            else bufs->cons_quals[i] = 0, bufs->cons_seq_buffer[i] = 'N';
        } else if(rminind == 4) { // Forward is N'd and reverse is not. Reverse call is probably right.
            index = i * 5 + fminind;
            kfpf->phred_sums[index] += kfpr->phred_sums[index];
            kfpf->nuc_counts[index] += kfpr->nuc_counts[index];
            diffs -= kfpf->nuc_counts[index];
            for(j = 0; j < 5; ++j) {
                fpsums[j] = pvalue_to_phred(igamc_pvalues(kfpf->length, LOG10_TO_CHI2((kfpf->phred_sums[j] + kfpr->phred_sums[j]))));
            }
            for(j = 0; j < 5; ++j) if(j != fminind) fpsums[fminind] -= fpsums[j];
            if(fpsums[fminind] < 0) fpsums[fminind] = 0;
            bufs->cons_quals[i] = fpsums[fminind];
            bufs->agrees[i] = kfpf->nuc_counts[fminind + i * 5] + kfpr->nuc_counts[fminind + i * 5];
            diffs -= bufs->agrees[i];
            diffs -= kfpf->nuc_counts[4 + i * 5] + kfpr->nuc_counts[4 + i * 5];
            if(bufs->cons_quals[i] > 2 && (double)bufs->agrees[i] / kfpr->length > MIN_FRAC_AGREED) bufs->cons_seq_buffer[i] = num2nuc(fminind);
            else bufs->cons_quals[i] = 0, bufs->cons_seq_buffer[i] = 'N';
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
