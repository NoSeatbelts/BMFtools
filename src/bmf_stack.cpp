#include "bmf_stack.h"
namespace BMF {

    const char *vcf_header_lines[] =  {
            "##FORMAT=<ID=FR_FAILED,Number=1,Type=Integer,Description=\"Number of observations failed per sample for fraction agreed.\">",
            "##FORMAT=<ID=FM_FAILED,Number=1,Type=Integer,Description=\"Number of observations failed per sample for family size.\">",
            "##FORMAT=<ID=FA_FAILED,Number=1,Type=Integer,Description=\"Number of observations failed per sample for number of supporting observations.\">",
            "##FORMAT=<ID=FP_FAILED,Number=1,Type=Integer,Description=\"Number of observations failed per sample for being a barcode QC fail.\">",
            "##FORMAT=<ID=AF_FAILED,Number=1,Type=Integer,Description=\"Number of observations failed per sample for aligned fraction below minimm.\">",
            "##FORMAT=<ID=MQ_FAILED,Number=1,Type=Integer,Description=\"Number of observations failed per sample for insufficient mapping quality.\">",
            "##FORMAT=<ID=IMPROPER,Number=1,Type=Integer,Description=\"Number of reads per sample labeled as not being in a proper pair.\">",
            "##FORMAT=<ID=OVERLAP,Number=1,Type=Integer,Description=\"Number of overlapping read pairs.\">"
            "##FORMAT=<ID=AFR,Number=R,Type=Float,Description=\"Allele fractions per allele, including the reference allele.\">"
    };


    void stack_usage(int retcode)
    {
        fprintf(stderr,
                        "Builds a stack summary for base calls and performs simple filters"
                        " to produce a maximally permissive 'variant caller'.\n"
                        "Usage:\nbmftools stack <opts> <in.srt.indexed.bam>\n"
                        "Optional arguments:\n"
                        "-R, --refpath\tPath to fasta reference. REQUIRED.\n"
                        "-o, --outpath\tPath to output file. Defaults to stdout.\n"
                        "-b, --bedpath\tPath to bed file to only validate variants in said region. REQUIRED.\n"
                        "-c, --min-count\tMinimum number of observations for a given allele passing filters to pass variant. Default: 1.\n"
                        "-s, --min-family-size\tMinimum number of reads in a family to include a that collapsed observation\n"
                        "-2, --skip-secondary\tSkip secondary alignments.\n"
                        "-S, --skip-supplementary\tSkip supplementary alignments.\n"
                        "-q, --skip-qcfail\tSkip reads marked as QC fail.\n"
                        "-r, --skip-duplicates\tSkip reads marked as being PCR or optical duplicates.\n"
                        "-f, --min-fraction-agreed\tMinimum fraction of reads in a family agreed on a base call\n"
                        "-v, --min-phred-quality\tMinimum calculated p-value on a base call in phred space\n"
                        "-p, --padding\tNumber of bases outside of bed region to pad.\n"
                        "-a, --min-family-agreed\tMinimum number of reads in a family agreed on a base call\n"
                        "-m, --min-mapping-quality\tMinimum mapping quality for reads for inclusion\n"
                        "-B, --emit-bcf-format\tEmit bcf-formatted output. (Defaults to vcf).\n"
                );
        exit(retcode);
    }


    static int read_bam(dlib::BamHandle *data, bam1_t *b)
    {
        int ret;
        for(;;)
        {
            if(!data->iter) LOG_EXIT("Need to access bam with index.\n");
            ret = sam_itr_next(data->fp, data->iter, b);
            if ( ret<0 ) break;
            // Skip unmapped, secondary, qcfail, duplicates.
            if ((b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP)) == 0)
                break;
        }
        return ret;
    }

    void process_matched_pileups(BMF::stack_aux_t *aux, bcf1_t *ret,
                            const int& tn_plp, const int& tpos, const int& ttid,
                            const int& nn_plp, const int& npos, const int& ntid) {
        std::string qname;
        // Build overlap hash
        std::unordered_map<std::string, BMF::UniqueObservation> tobs, nobs;
        std::unordered_map<std::string, BMF::UniqueObservation>::iterator found;
        int flag_failed[2] = {0};
        int af_failed[2] = {0};
        int fa_failed[2] = {0};
        int fm_failed[2] = {0};
        int fr_failed[2] = {0};
        int fp_failed[2] = {0};
        int mq_failed[2] = {0};
        int improper_count[2] = {0};
        int olap_count[2] = {0};
        // Capturing found  by reference to avoid making unneeded temporary variables.
        std::for_each(aux->tumor.pileups, aux->tumor.pileups + tn_plp, [&](const bam_pileup1_t& plp) {
            if(plp.is_del || plp.is_refskip) return;
            if(aux->conf.skip_flag & plp.b->core.flag) {
                ++flag_failed[0];
                return;
            }
            if((plp.b->core.flag & BAM_FPROPER_PAIR) == 0) {
                ++improper_count[0];
                if(aux->conf.skip_improper) return;
            }
            if(bam_itag(plp.b, "FP") == 0) {
                ++fp_failed[0]; return;
            }
            if(bam_aux2f(bam_aux_get(plp.b, "AF")) < aux->conf.minAF) {
                ++af_failed[0]; return;
            }
            //LOG_DEBUG("Now changing qname (%s) to new qname (%s).\n", qname.c_str(), bam_get_qname(plp.b));
            qname = std::string(bam_get_qname(plp.b));
            //LOG_DEBUG("Changed qname (%s) to new qname (%s).\n", qname.c_str(), bam_get_qname(plp.b));
            if((found = tobs.find(qname)) == tobs.end()) {
                tobs.emplace(qname, plp);
            } else {
                ++olap_count[0];
                //LOG_DEBUG("Added other in pair with qname %s.\n", qname.c_str());
                found->second.add_obs(plp);
            }
        });
        for(auto& pair: tobs) {
            if(pair.second.size < aux->conf.minFM) ++fm_failed[0], pair.second.set_pass(0);
            if(pair.second.agreed < aux->conf.minFA) ++fa_failed[0], pair.second.set_pass(0);
            if((float)pair.second.agreed / pair.second.size < aux->conf.minFR) ++fr_failed[0], pair.second.set_pass(0);
            if(pair.second.get_meanMQ() < aux->conf.minMQ) ++mq_failed[0], pair.second.set_pass(0);
        }
        std::for_each(aux->normal.pileups, aux->normal.pileups + nn_plp, [&](const bam_pileup1_t& plp) {
            if(plp.is_del || plp.is_refskip) return;
            if(aux->conf.skip_flag & plp.b->core.flag) {
                ++flag_failed[1];
                return;
            }
            if((plp.b->core.flag & BAM_FPROPER_PAIR) == 0) {
                ++improper_count[1];
                if(aux->conf.skip_improper) return;
            }
            if(bam_itag(plp.b, "FP") == 0) {
                ++fp_failed[1]; return;
            }
            if(bam_aux2f(bam_aux_get(plp.b, "AF")) < aux->conf.minAF) {
                ++af_failed[1]; return;
            }
            qname = std::string(bam_get_qname(plp.b));
            if((found = nobs.find(qname)) == nobs.end()) {
                nobs.emplace(qname, plp);
            } else {
                ++olap_count[1];
                found->second.add_obs(plp);
            }
        });
        for(auto& pair: nobs) {
            if(pair.second.size < aux->conf.minFM) ++fm_failed[1], pair.second.set_pass(0);
            if(pair.second.agreed < aux->conf.minFA) ++fa_failed[1], pair.second.set_pass(0);
            if((float)pair.second.agreed / pair.second.size < aux->conf.minFR) ++fr_failed[1], pair.second.set_pass(0);
            if(pair.second.get_meanMQ() < aux->conf.minMQ) ++mq_failed[1], pair.second.set_pass(0);
        }
        LOG_DEBUG("Making PairVCFPos.\n");
        // Build vcfline struct
        BMF::PairVCFPos vcfline = BMF::PairVCFPos(tobs, nobs, ttid, tpos);
        LOG_DEBUG("Making bcf.\n");
        vcfline.to_bcf(ret, aux, ttid, tpos);
        bcf_update_format_int32(aux->vcf.vh, ret, "FR_FAILED", (void *)fr_failed, 2);
        bcf_update_format_int32(aux->vcf.vh, ret, "FA_FAILED", (void *)fa_failed, 2);
        bcf_update_format_int32(aux->vcf.vh, ret, "FP_FAILED", (void *)fp_failed, 2);
        bcf_update_format_int32(aux->vcf.vh, ret, "FM_FAILED", (void *)fm_failed, 2);
        bcf_update_format_int32(aux->vcf.vh, ret, "MQ_FAILED", (void *)mq_failed, 2);
        bcf_update_format_int32(aux->vcf.vh, ret, "AF_FAILED", (void *)af_failed, 2);
        bcf_update_format_int32(aux->vcf.vh, ret, "IMPROPER", (void *)improper_count, 2);
        bcf_update_format_int32(aux->vcf.vh, ret, "OVERLAP", (void *)olap_count, 2);
        //LOG_INFO("Ret for writing vcf to file: %i.\n", aux->vcf.write(ret));
        aux->vcf.write(ret);
        bcf_clear(ret);
        LOG_DEBUG("Finished processing matched pileups.\n");
    }

    /*
     * Needs a rewrite after the T/N pair rewrite!
     */
    void process_pileup(bcf1_t *ret, const bam_pileup1_t *plp, int n_plp, int pos, int tid, BMF::stack_aux_t *aux) {
        std::string qname;
        // Build overlap hash
        std::unordered_map<std::string, BMF::UniqueObservation> obs;
        std::unordered_map<std::string, BMF::UniqueObservation>::iterator found;
        int flag_failed = 0;
        int af_failed = 0;
        int fa_failed = 0;
        int fm_failed = 0;
        int fp_failed = 0;
        int fr_failed = 0;
        int mq_failed = 0;
        int improper_count = 0;
        // Capturing found  by reference to avoid making unneeded temporary variables.
        std::for_each(plp, plp + n_plp, [&](const bam_pileup1_t& plp) {
            if(plp.is_del || plp.is_refskip) return;
            if(aux->conf.skip_flag & plp.b->core.flag) {
                ++flag_failed;
                return;
            }
            if((plp.b->core.flag & BAM_FPROPER_PAIR) == 0) {
                ++improper_count;
                if(aux->conf.skip_improper) return;
            }
            // If a read's mate is here with a sufficient mapping quality, we should keep it, shouldn't we? Put this later.
            if(plp.b->core.qual < aux->conf.minMQ) {
                ++mq_failed; return;
            }
            const int FM = bam_itag(plp.b, "FM");
            if(bam_itag(plp.b, "FP") == 0) {
                ++fp_failed; return;
            }
            if(bam_aux2f(bam_aux_get(plp.b, "AF")) < aux->conf.minAF) {
                ++af_failed; return;
            }
            const int qpos = dlib::arr_qpos(&plp);
            const uint32_t FA = ((uint32_t *)dlib::array_tag(plp.b, "FA"))[qpos];
            if(FA < aux->conf.minFA) {
                ++fa_failed; return;
            }
            if((float)FA / FM < aux->conf.minFR) {
                ++fr_failed; return;
            }
            // Should I be failing FA/FM/PV before merging overlapping reads? NO.
            qname = bam_get_qname(plp.b);
            LOG_DEBUG("Got qname %s.\n", qname.c_str());
            if((found = obs.find(qname)) == obs.end())
                obs.emplace(std::make_pair(qname, BMF::UniqueObservation(plp)));
            else found->second.add_obs(plp);
        });
        // Build vcfline struct
        BMF::SampleVCFPos vcfline = BMF::SampleVCFPos(obs, tid, pos);
        vcfline.to_bcf(ret, aux->vcf.vh, aux->get_ref_base(tid, pos));
        bcf_update_info_int32(aux->vcf.vh, ret, "FR_FAILED", (void *)&fr_failed, 1);
        bcf_update_info_int32(aux->vcf.vh, ret, "FA_FAILED", (void *)&fa_failed, 1);
        bcf_update_info_int32(aux->vcf.vh, ret, "FP_FAILED", (void *)&fp_failed, 1);
        bcf_update_info_int32(aux->vcf.vh, ret, "FM_FAILED", (void *)&fm_failed, 1);
        bcf_update_info_int32(aux->vcf.vh, ret, "MQ_FAILED", (void *)&mq_failed, 1);
        bcf_update_info_int32(aux->vcf.vh, ret, "AF_FAILED", (void *)&af_failed, 1);
        bcf_update_info_int32(aux->vcf.vh, ret, "IMPROPER", (void *)&improper_count, 1);
        aux->vcf.write(ret);
        bcf_clear(ret);
    }

    int stack_core(BMF::stack_aux_t *aux)
    {
        if(!aux->tumor.idx || !aux->normal.idx)
            LOG_EXIT("Could not load bam indices. Abort!\n");
        aux->tumor.plp = bam_plp_init((bam_plp_auto_f)read_bam, (void *)&aux->tumor);
        aux->normal.plp = bam_plp_init((bam_plp_auto_f)read_bam, (void *)&aux->normal);
        bam_plp_set_maxcnt(aux->tumor.plp, aux->conf.max_depth);
        bam_plp_set_maxcnt(aux->normal.plp, aux->conf.max_depth);
        std::vector<khiter_t> sorted_keys(dlib::make_sorted_keys(aux->bed));
        int ttid, tpos, tn_plp, ntid, npos, nn_plp;
        bcf1_t *v = bcf_init1();
        for(unsigned k = 0; k < sorted_keys.size(); ++k) {
            const khiter_t key = sorted_keys[k];
            const size_t n = kh_val(aux->bed, key).n;
            for(uint64_t i = 0; i < n; ++i) {
                const int start = get_start(kh_val(aux->bed, key).intervals[i]);
                const int stop = get_stop(kh_val(aux->bed, key).intervals[i]);
                const int bamtid = (int)kh_key(aux->bed, key);
                aux->pair_region_itr(bamtid, start, stop);
                process_matched_pileups(aux, v, tn_plp, tpos, ttid, nn_plp, npos, ntid);
                while(aux->next_paired_pileup(&ttid, &tpos, &tn_plp, &ntid, &npos, &nn_plp, stop))
                    process_matched_pileups(aux, v, tn_plp, tpos, ttid, nn_plp, npos, ntid);
            }
        }
        bcf_destroy(v);
        return 0;
    }

    int stack_main(int argc, char *argv[]) {
        int c;
        unsigned padding = (unsigned)-1;
        if(argc < 2) stack_usage(EXIT_FAILURE);
        char *outvcf = NULL, *refpath = NULL;
        char *bedpath = NULL;
        struct BMF::stack_conf_t conf = {0};
        const struct option lopts[] = {
            {"skip-secondary", no_argument, NULL, '2'},
            {"min-family-agreed", required_argument, NULL, 'a'},
            {"bedpath", required_argument, NULL, 'b'},
            {"emit-bcf", no_argument, NULL, 'B'},
            {"min-count", required_argument, NULL, 'c'},
            {"max-depth", required_argument, NULL, 'd'},
            {"min-duplex", required_argument, NULL, 'D'},
            {"min-fraction-agreed", required_argument, NULL, 'f'},
            {"min-mapping-quality", required_argument, NULL, 'm'},
            {"min-overlap", required_argument, NULL, 'O'},
            {"out-vcf", required_argument, NULL, 'o'},
            {"skip-improper", no_argument, NULL, 'P'},
            {"padding", required_argument, NULL, 'p'},
            {"skip-qcfail", no_argument, NULL, 'q'},
            {"skip-duplicates", required_argument, NULL, 'r'},
            {"ref", required_argument, NULL, 'R'},
            {"min-family-size", required_argument, NULL, 's'},
            {"skip-supplementary", no_argument, NULL, 'S'},
            {"min-phred-quality", required_argument, NULL, 'v'},
            {0, 0, 0, 0}
        };
        while ((c = getopt_long(argc, argv, "R:D:q:r:2:S:d:a:s:m:p:f:b:v:o:O:c:BP?hV", lopts, NULL)) >= 0) {
            switch (c) {
                case '2': conf.skip_flag |= BAM_FSECONDARY; break;
                case 'a': conf.minFA = atoi(optarg); break;
                case 'b': bedpath = optarg; break;
                case 'B': conf.output_bcf = 1; break;
                case 'c': conf.minCount = atoi(optarg); break;
                case 'd': conf.max_depth = atoi(optarg); break;
                case 'D': conf.minDuplex = atoi(optarg); break;
                case 'f': conf.minFR = (float)atof(optarg); break;
                case 'm': conf.minMQ = atoi(optarg); break;
                case 'O': conf.minOverlap = atoi(optarg); break;
                case 'o': outvcf = optarg; break;
                case 'P': conf.skip_improper = 1; break;
                case 'p': padding = atoi(optarg); break;
                case 'q': conf.skip_flag |= BAM_FQCFAIL; break;
                case 'r': conf.skip_flag |= BAM_FDUP; break;
                case 'R': refpath = optarg; break;
                case 's': conf.minFM = atoi(optarg); break;
                case 'S': conf.skip_flag |= BAM_FSUPPLEMENTARY; break;
                case 'v': conf.minPV = atoi(optarg); break;
                case 'h': case '?': stack_usage(EXIT_SUCCESS);
            }
        }
        if(optind >= argc - 1) LOG_EXIT("Insufficient arguments. Input bam required!\n");
        if(padding < 0) {
            LOG_WARNING("Padding not set. Using default %i.\n", DEFAULT_PADDING);
            padding = DEFAULT_PADDING;
        }
        if(!refpath) {
            LOG_EXIT("refpath required. Abort!\n");
        }
        if(!outvcf) outvcf = (char *)"-";
        bcf_hdr_t *vh = bcf_hdr_init(conf.output_bcf ? "wb": "w");
        for(auto line: vcf_header_lines) {
            LOG_DEBUG("Adding line %s.\n", line);
            if(bcf_hdr_append(vh, line))
                LOG_EXIT("Could not add line %s to header. Abort!\n", line);
        }
        for(auto line: stack_vcf_lines)
            if(bcf_hdr_append(vh, line))
                LOG_EXIT("Could not add line %s to header. Abort!\n", line);
        // Add samples
        int tmp;
        if((tmp = bcf_hdr_add_sample(vh, "Tumor")))
            LOG_EXIT("Could not add name %s. Code: %i.\n", "Tumor", tmp);
        if((tmp = bcf_hdr_add_sample(vh, "Normal")))
            LOG_EXIT("Could not add name %s. Code: %i.\n", "Normal", tmp);
        bcf_hdr_add_sample(vh, NULL);
        bcf_hdr_nsamples(vh) = 2;
        LOG_DEBUG("N samples: %i.\n", bcf_hdr_nsamples(vh));
        // Add lines to the header for the bed file?
        BMF::stack_aux_t aux(argv[optind], argv[optind + 1], outvcf, vh, conf);
        bcf_hdr_destroy(vh);
        if((aux.fai = fai_load(refpath)) == NULL)
            LOG_EXIT("failed to open fai. Abort!\n");
        // TODO: Make BCF header
        aux.bed = bedpath ? dlib::parse_bed_hash(bedpath, aux.normal.header, padding)
                          : dlib::build_ref_hash(aux.normal.header);
        // Check for required tags.
        for(auto tag: {"FM", "FA", "PV", "FP", "AF", "DR"})
            dlib::check_bam_tag_exit(aux.normal.fp->fn, tag);
        int ret = stack_core(&aux);
        LOG_INFO("Successfully complete bmftools stack!\n");
        return ret;
    }

} /* namespace BMF */