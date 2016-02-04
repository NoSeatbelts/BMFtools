#include "bmf_hashdmp.h"

#ifndef ifn_stream
#define ifn_stream(fname) ((fname) ? (fname): "stream")
#endif

void print_hash_dmp_usage(char *arg) {
	fprintf(stderr, "Usage: %s -o <output_filename> <input_filename>.\n"
			"Flags:\n"
			"-s\tPerform secondary index consolidation rather than Loeb-like inline consolidation.\n"
			"If output file is unset, defaults to stdout. If input filename is not set, defaults to stdin.\n"
			, arg);
}

void print_hash_dmp_opt_err(char *arg, char *optarg, char optopt) {
	fprintf(stderr, "[E:%s] Invalid flag '%c' with argument '%s'. See usage.\n",
			__func__, optopt, optarg);
	print_hash_dmp_usage(arg);
	exit(EXIT_FAILURE);
}


tmpvars_t *init_tmpvars_p(char *bs_ptr, int blen, int readlen)
{
	tmpvars_t *ret = (tmpvars_t *)malloc(sizeof(tmpvars_t));
	ret->blen = blen;
	ret->readlen = readlen;
	ret->bs_ptr = bs_ptr;
	ret->buffers = (tmpbuffers_t *)malloc(sizeof(tmpbuffers_t));
	ret->buffers->name_buffer[0] = '@';
	ret->buffers->name_buffer[blen] = '\0';
	ret->buffers->cons_seq_buffer[readlen] = '\0';
	return ret;
}



int hash_dmp_main(int argc, char *argv[])
{
	if(argc == 1) print_hash_dmp_usage(argv[0]), exit(EXIT_SUCCESS);
	char *outfname = NULL, *infname = NULL;
	int c;
	int stranded_analysis = 1;
	int level = -1;
	while ((c = getopt(argc, argv, "l:o:sh?")) >= 0) {
		switch(c) {
			case 'l': level = atoi(optarg); break;
			case 'o': outfname = strdup(optarg); break;
			case 's': stranded_analysis = 0; break;
			case '?': // Fall-through
			case 'h': print_hash_dmp_usage(argv[0]); return EXIT_SUCCESS;
			default: print_hash_dmp_opt_err(argv[0], optarg, optopt);
		}
	}
	if(argc < 2) {
		fprintf(stderr, "[E:%s] Required arguments missing. See usage.\n", __func__);
		print_hash_dmp_usage(argv[0]);
		exit(1);
	}
	if(argc - 1 == optind) infname = strdup(argv[optind]);
	else {
		LOG_WARNING("Note: no input filename provided. Defaulting to stdin.\n");
	}
	stranded_analysis ? stranded_hash_dmp_core(infname, outfname, level) :hash_dmp_core (infname, outfname, level);
	cond_free(outfname); cond_free(infname);
	return EXIT_SUCCESS;
}


void duplex_hash_process(hk_t *hfor, hk_t *cfor, hk_t *tmp_hkf, hk_t *crev, hk_t *hrev, gzFile out_handle, tmpvars_t *tmp)
{
	HASH_ITER(hh, hfor, cfor, tmp_hkf) {
		HASH_FIND_STR(hrev, cfor->id, crev);
		if(crev) {
			zstranded_process_write(cfor->value, crev->value, out_handle, tmp->buffers); // Found from both strands!
			destroy_kf(cfor->value), destroy_kf(crev->value);
			HASH_DEL(hrev, crev); HASH_DEL(hfor, cfor);
			free(crev); free(cfor);
		} else {
			kdmp_process_write(cfor->value, out_handle, tmp->buffers, 0); // No reverse strand found. \='{
			destroy_kf(cfor->value);
			HASH_DEL(hfor, cfor);
			free(cfor);
		}
	}
}

void duplex_hash_fill(kseq_t *seq, hk_t *hfor, hk_t *tmp_hkf, hk_t *hrev, hk_t *tmp_hkr, char *infname, uint64_t *count, uint64_t *fcount, tmpvars_t *tmp, int blen) {
	if(UNLIKELY(++*count % 1000000 == 0))
		fprintf(stderr, "[%s::%s] Number of records processed: %lu.\n", __func__,
				*infname == '-' ? "stdin" : infname, *count);
	if(seq->comment.s[HASH_DMP_OFFSET] == 'F') {
		++*fcount;
		cp_view2buf(seq->comment.s + HASH_DMP_OFFSET + 1, tmp->key);
		HASH_FIND_STR(hfor, tmp->key, tmp_hkf);
		if(!tmp_hkf) {
			tmp_hkf = (hk_t *)malloc(sizeof(hk_t));
			tmp_hkf->value = init_kfp(tmp->readlen);
			cp_view2buf(seq->comment.s + HASH_DMP_OFFSET + 1, tmp_hkf->id);
			pushback_kseq(tmp_hkf->value, seq, blen);
			HASH_ADD_STR(hfor, id, tmp_hkf);
		} else pushback_kseq(tmp_hkf->value, seq, blen);
	} else {
		cp_view2buf(seq->comment.s + HASH_DMP_OFFSET + 1, tmp->key);
		HASH_FIND_STR(hrev, tmp->key, tmp_hkr);
		if(!tmp_hkr) {
			tmp_hkr = (hk_t *)malloc(sizeof(hk_t));
			tmp_hkr->value = init_kfp(tmp->readlen);
			cp_view2buf(seq->comment.s + HASH_DMP_OFFSET + 1, tmp_hkr->id);
			pushback_kseq(tmp_hkr->value, seq, blen);
			HASH_ADD_STR(hrev, id, tmp_hkr);
		} else pushback_kseq(tmp_hkr->value, seq, blen);
	}
}


void se_hash_process(hk_t *hash, hk_t *current_entry, hk_t *tmp_hk, gzFile out_handle, tmpvars_t *tmp, uint64_t *count)
{
	LOG_DEBUG("Beginning se_hash_process with count %lu.\n", *count);
	HASH_ITER(hh, hash, current_entry, tmp_hk) {
		if(++(*count) % 1000000) {
			LOG_DEBUG("Number of records written: %lu.\n", *count);
		}
		kdmp_process_write(current_entry->value, out_handle, tmp->buffers, 0);
		destroy_kf(current_entry->value);
		HASH_DEL(hash, current_entry);
		free(current_entry);
	}
}


void hash_dmp_core(char *infname, char *outfname, int level)
{
	char mode[4] = "wT";
	if(level >= 0) sprintf(mode, "wb%i", level % 10);
	LOG_DEBUG("Writing hash dmp information with mode: '%s'.", mode);
	FILE *in_handle = open_ifp(infname);
	gzFile out_handle = gzopen(outfname, mode);
	if(!in_handle) {
		fprintf(stderr, "[E:%s] Could not open %s for reading. Abort mission!\n", __func__, infname);
		exit(EXIT_FAILURE);
	}
	gzFile fp = gzdopen(fileno(in_handle), "r");
	kseq_t *seq = kseq_init(fp);
	// Initialized kseq
	int l = kseq_read(seq);
	if(l < 0) {
		fprintf(stderr, "[E:%s]: Could not open fastq file (%s). Abort mission!\n",
				__func__, strcmp(infname, "-") == 0 ? "stdin": infname);
		exit(EXIT_FAILURE);
	}
	char *bs_ptr = barcode_mem_view(seq);
	const int blen = infer_barcode_length(bs_ptr);
#if !NDEBUG
	fprintf(stderr, "[D:%s] Barcode length (inferred): %i.\n", __func__, blen);
#endif
	tmpvars_t *tmp = init_tmpvars_p(bs_ptr, blen, seq->seq.l);
	memcpy(tmp->key, bs_ptr, blen);
	tmp->key[blen] = '\0';
	// Start hash table
	hk_t *hash = NULL;
	hk_t *current_entry = (hk_t *)malloc(sizeof(hk_t));
	hk_t *tmp_hk = current_entry; // Save the pointer location for later comparison.
	cp_view2buf(bs_ptr + 1, current_entry->id);
	current_entry->value = init_kfp(tmp->readlen);
	HASH_ADD_STR(hash, id, current_entry);
	pushback_kseq(current_entry->value, seq, blen);

	uint64_t count = 0;
	while(LIKELY((l = kseq_read(seq)) >= 0)) {
		if(UNLIKELY(++count % 1000000 == 0))
			fprintf(stderr, "[%s::%s] Number of records read: %lu.\n", __func__,
					strcmp("-", infname) == 0 ? "stdin": infname,count);
		cp_view2buf(seq->comment.s + HASH_DMP_OFFSET + 1, tmp->key);
		HASH_FIND_STR(hash, tmp->key, tmp_hk);
		if(!tmp_hk) {
			tmp_hk = (hk_t *)malloc(sizeof(hk_t));
			tmp_hk->value = init_kfp(tmp->readlen);
			cp_view2buf(seq->comment.s + HASH_DMP_OFFSET + 1, tmp_hk->id);
			pushback_kseq(tmp_hk->value, seq, blen);
			HASH_ADD_STR(hash, id, tmp_hk);
		} else pushback_kseq(tmp_hk->value, seq, blen);
	}
	fprintf(stderr, "[%s::%s] Loaded all records into memory. Writing out to file!\n", __func__, ifn_stream(infname));
	count = 0;
	se_hash_process(hash, current_entry, tmp_hk, out_handle, tmp, &count);
	fprintf(stderr, "[%s::%s] Total number of collapsed observations: %lu.\n", __func__, ifn_stream(infname), count);
	gzclose(fp);
	gzclose(out_handle);
	kseq_destroy(seq);
	tmpvars_destroy(tmp);
}

void stranded_hash_dmp_core(char *infname, char *outfname, int level)
{
	char mode[4] = "wT";
	if(level >= 0) sprintf(mode, "wb%i", level % 10);
	LOG_DEBUG("Writing stranded hash dmp information with mode: '%s'.", mode);
	FILE *in_handle = open_ifp(infname);
	gzFile out_handle = gzopen(outfname, mode);
	if(!in_handle) {
		fprintf(stderr, "[E:%s] Could not open %s for reading. Abort mission!\n", __func__, infname);
		exit(EXIT_FAILURE);
	}
	gzFile fp = gzdopen(fileno(in_handle), "r");
	kseq_t *seq = kseq_init(fp);
	// Initialized kseq
	int l = kseq_read(seq);
	if(l < 0) {
		fprintf(stderr, "[%s]: Could not open fastq file (%s). Abort mission!\n",
				__func__, strcmp(infname, "-") == 0 ? "stdin": infname);
		exit(EXIT_FAILURE);
	}
	char *bs_ptr = barcode_mem_view(seq);
	int blen = infer_barcode_length(bs_ptr);
	LOG_DEBUG("Barcode length (inferred): %i. First barcode: %s.\n", blen, bs_ptr);
	tmpvars_t *tmp = init_tmpvars_p(bs_ptr, blen, seq->seq.l);
	memcpy(tmp->key, bs_ptr, blen);
	tmp->key[blen] = '\0';
	// Start hash table
	hk_t *hfor = NULL, *hrev = NULL; // Hash forward, hash reverse
	hk_t *crev = (hk_t *)malloc(sizeof(hk_t)); // Current reverse, current forward.
	hk_t *cfor = (hk_t *)malloc(sizeof(hk_t));
	hk_t *tmp_hkr = crev, *tmp_hkf = cfor;
	uint64_t count = 1, fcount = 0;
	if(*bs_ptr == 'F') {
		++fcount;
		cp_view2buf(bs_ptr + 1, cfor->id);
		cfor->value = init_kfp(tmp->readlen);
		HASH_ADD_STR(hfor, id, cfor);
		pushback_kseq(cfor->value, seq, blen);
	} else {
		cp_view2buf(bs_ptr + 1, crev->id);
		crev->value = init_kfp(tmp->readlen);
		HASH_ADD_STR(hrev, id, crev);
		pushback_kseq(hrev->value, seq, blen);
	}

	while(LIKELY((l = kseq_read(seq)) >= 0)) {
		if(UNLIKELY(++count % 1000000 == 0))
			fprintf(stderr, "[%s::%s] Number of records processed: %lu.\n", __func__,
					*infname == '-' ? "stdin" : infname, count);
		if(seq->comment.s[HASH_DMP_OFFSET] == 'F') {
			++fcount;
			cp_view2buf(seq->comment.s + HASH_DMP_OFFSET + 1, tmp->key);
			HASH_FIND_STR(hfor, tmp->key, tmp_hkf);
			if(!tmp_hkf) {
				tmp_hkf = (hk_t *)malloc(sizeof(hk_t));
				tmp_hkf->value = init_kfp(tmp->readlen);
				cp_view2buf(seq->comment.s + HASH_DMP_OFFSET + 1, tmp_hkf->id);
				pushback_kseq(tmp_hkf->value, seq, blen);
				HASH_ADD_STR(hfor, id, tmp_hkf);
			} else pushback_kseq(tmp_hkf->value, seq, blen);
		} else {
			cp_view2buf(seq->comment.s + HASH_DMP_OFFSET + 1, tmp->key);
			HASH_FIND_STR(hrev, tmp->key, tmp_hkr);
			if(!tmp_hkr) {
				tmp_hkr = (hk_t *)malloc(sizeof(hk_t));
				tmp_hkr->value = init_kfp(tmp->readlen);
				cp_view2buf(seq->comment.s + HASH_DMP_OFFSET + 1, tmp_hkr->id);
				pushback_kseq(tmp_hkr->value, seq, blen);
				HASH_ADD_STR(hrev, id, tmp_hkr);
			} else pushback_kseq(tmp_hkr->value, seq, blen);
		}
	}
	uint64_t rcount = count - fcount;
	LOG_INFO("Number of reverse reads: %lu. Number of forward reads: %lu.\n", rcount, fcount);

	fprintf(stderr, "[%s::%s] Loaded all records into memory. Writing out to file!\n", __func__, ifn_stream(outfname));
	// Write out all unmatched in forward and handle all barcodes handled from both strands.
	uint64_t duplex = 0, non_duplex = 0;
	HASH_ITER(hh, hfor, cfor, tmp_hkf) {
		HASH_FIND_STR(hrev, cfor->id, crev);
		if(!crev) {
			++non_duplex;
			dmp_process_write(cfor->value, out_handle, tmp->buffers, 0); // No reverse strand found. \='{
			destroy_kf(cfor->value);
			HASH_DEL(hfor, cfor);
			free(cfor);
			continue;
		}
		++duplex;
		zstranded_process_write(cfor->value, crev->value, out_handle, tmp->buffers); // Found from both strands!
		destroy_kf(cfor->value), destroy_kf(crev->value);
		HASH_DEL(hrev, crev); HASH_DEL(hfor, cfor);
		cond_free(crev); cond_free(cfor);
	}
	LOG_DEBUG("Before handling reverse only counts for non_duplex: %lu.\n", non_duplex);
	//duplex_hash_process(hfor, cfor, tmp_hkf, crev, hrev, out_handle, tmp);
	se_hash_process(hrev, crev, tmp_hkr, out_handle, tmp, &non_duplex);
	LOG_DEBUG("Cleaning up.\n");
	LOG_INFO("Number of duplex observations: %lu. Number of non-duplex observations: %lu\n", duplex, non_duplex);
	gzclose(fp);
	fclose(in_handle), gzclose(out_handle);
	kseq_destroy(seq);
	tmpvars_destroy(tmp);
}
