#include "kingfisher.h"

#define dmp_pos(kfp, bufs, argmaxret, i, index)\
	bufs->cons_quals[i] = pvalue_to_phred(igamc_pvalues(kfp->length, LOG10_TO_CHI2((kfp->phred_sums[index]))));\
	bufs->agrees[i] = kfp->nuc_counts[index];\
	bufs->cons_seq_buffer[i] = (bufs->cons_quals[i] > 2 && (double)bufs->agrees[i] / kfp->length > MIN_FRAC_AGREED) ?\
		ARRG_MAX_TO_NUC(argmaxret): 'N';\
    if(bufs->cons_seq_buffer[i] == 'N') bufs->cons_quals[i] = 2

// TODO: rewrite dmp_process_write and stranded_process_write to write the PV/FA strings straight to output
// rather than writing to a temporary object and writing that later.

void dmp_process_write(KingFisher_t *kfp, FILE *handle, tmpbuffers_t *bufs, int is_rev)
{
	int i;
	for(i = 0; i < kfp->readlen; ++i) {
		const int argmaxret = ARRG_MAX(kfp, i);
		const int index = argmaxret + i * 5;
		dmp_pos(kfp, bufs, argmaxret, i, index);
	}
	fill_fa_buffer(kfp->readlen, bufs->agrees, bufs->FABuffer);
	fill_pv_buffer(kfp->readlen, bufs->cons_quals, bufs->PVBuffer);
#if PUTC
	fputc('@', handle);
	fputs(kfp->barcode,handle);
	fputc(' ', handle);
	fputs(bufs->FABuffer, handle);
	fputc('\t', handle);
	fputs(bufs->PVBuffer, handle);
	fputs("\tFP:i:", handle);
	fputc(kfp->pass_fail, handle);
	fprintf(handle, "\tFM:i:%i", kfp->length);
	fprintf(handle, "\tRV:i:%i\n", is_rev ? kfp->length: 0);
	fputs(bufs->cons_seq_buffer, handle);
	fputs("\n+\n", handle);
	for(i = 0; i < kfp->readlen; ++i)
		fputc(kfp->max_phreds[nuc2num(bufs->cons_seq_buffer[i]) + i * 5], handle);
	fputc('\n', handle);
#else
	fprintf(handle, "@%s %s\t%s\tFP:i:%c\tFM:i:%i\tRV:i:%i\n%s\n+\n", kfp->barcode + 1,
			bufs->FABuffer, bufs->PVBuffer,
			kfp->pass_fail, kfp->length, is_rev ? kfp->length: 0,
			bufs->cons_seq_buffer);
	for(i = 0; i < kfp->readlen; ++i)
		fputc(kfp->max_phreds[nuc2num(bufs->cons_seq_buffer[i]) + 5 * i], handle);
	fputc('\n', handle);
#endif
	return;
}

// kfp forward, kfp reverse
void stranded_process_write(KingFisher_t *kfpf, KingFisher_t *kfpr, FILE *handle, tmpbuffers_t *bufs)
{
	int i;
	for(i = 0; i < kfpf->readlen; ++i) {
		const int argmaxretf = ARRG_MAX(kfpf, i), argmaxretr = ARRG_MAX(kfpr, i);
        if(argmaxretf == argmaxretr) { // Both strands supported the same base call.
        	const int index = i * 5 + argmaxretf;
            kfpf->phred_sums[index] += kfpr->phred_sums[index];
            kfpf->nuc_counts[index] += kfpr->nuc_counts[index];
            dmp_pos(kfpf, bufs, argmaxretf, i, index);
            if(kfpr->max_phreds[index] > kfpf->max_phreds[index])
            	kfpf->max_phreds[index] = kfpr->max_phreds[index];
        }
        else if(argmaxretf == 4) { // Forward is Nd and reverse is not. Reverse call is probably right.
        	const int index = i * 5 + argmaxretr;
            kfpf->phred_sums[index] += kfpr->phred_sums[index];
            kfpf->nuc_counts[index] += kfpr->nuc_counts[index];
            dmp_pos(kfpf, bufs, argmaxretr, i, index);
            kfpf->max_phreds[index] = kfpr->max_phreds[index];
        }
        else if(argmaxretr == 4) { // Forward is Nd and reverse is not. Reverse call is probably right.
        	const int index = i * 5 + argmaxretf;
            kfpf->phred_sums[index] += kfpr->phred_sums[index];
            kfpf->nuc_counts[index] += kfpr->nuc_counts[index];
            dmp_pos(kfpf, bufs, argmaxretf, i, index);
            // Don't update max_phreds, since the max phred is already here.
        }
        else
            bufs->cons_quals[i] = 0, bufs->agrees[i] = 0, bufs->cons_seq_buffer[i] = 'N';
	}
	fill_fa_buffer(kfpf->readlen, bufs->agrees, bufs->FABuffer);
	//fprintf(stderr, "FA buffer: %s.\n", FABuffer);
	fill_pv_buffer(kfpf->readlen, bufs->cons_quals, bufs->PVBuffer);
	fprintf(handle, "@%s %s\t%s\tFP:i:%c\tFM:i:%i\tRV:i:%i\n%s\n+\n", kfpf->barcode + 1,
			bufs->FABuffer, bufs->PVBuffer,
			kfpf->pass_fail, kfpf->length, kfpr->length,
			bufs->cons_seq_buffer);
	for(i = 0; i < kfpf->readlen; ++i)
		fputc(kfpf->max_phreds[nuc2num(bufs->cons_seq_buffer[i]) + 5 * i], handle);
	fputc('\n', handle);
	return;
}

KingFisher_t *init_kfp(size_t readlen)
{
	KingFisher_t *ret = (KingFisher_t *)calloc(1, sizeof(KingFisher_t));
	ret->readlen = readlen;
	ret->max_phreds = (char *)malloc((readlen * 5) * sizeof(char));
	memset(ret->max_phreds, '#', readlen * 5);
	ret->nuc_counts = (uint16_t *)calloc(readlen * 5, sizeof(uint16_t));
	ret->phred_sums = (uint32_t *)calloc(readlen * 5, sizeof(uint32_t));
	ret->pass_fail = '1';
	return ret;
}

void destroy_kf(KingFisher_t *kfp)
{
	cond_free(kfp->nuc_counts);
	cond_free(kfp->phred_sums);
	cond_free(kfp->max_phreds);
	free(kfp);
	kfp = NULL;
}