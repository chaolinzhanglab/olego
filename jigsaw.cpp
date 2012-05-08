#include <stdio.h>
//#include <unistd.h>
#include <getopt.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
//#include "bwtaln.h"
#include "jigsawaln.h"
#include "bwtgap.h"
#include "utils.h"

int main (int argc, char *argv[])
{
	int c, opte = -1;
	gap_opt_t *opt;
	//char *bwa_rg = 0;
	opt = gap_init_opt();

	/*
	static struct option long_options [] = {
		{"word-size", 		required_argument, 0, 'w'},
		{"max-word-occ",	required_argument, 0, 'W'},
		{"junction-file",	required_argument, 0, 'j'},
		{"max-diff", 		required_argument, 0, 'n'},
		{"max-gapo", 		required_argument, 0, 'o'},
		{"max-gape", 		required_argument, 0, 'e'},
		{"max-intron-size",	required_argument, 0, 'I'},
		{"indel-end-skip", 	required_argument, 0, 'i'},
		{"gape-max-occ", 	required_argument, 0, 'd'},
		{"max-entries", 	required_argument, 0, 'm'},
		{"num-threads",		required_argument, 0, 't'},
		{"penalty-mismatch",required_argument, 0, 'M'},
		{"penalty-gapo",	required_argument, 0, 'O'},
		{"penalty-gape",	required_argument, 0, 'E'},
		{"repeat",			required_argument, 0, 'R'},
		{"rg-line",			required_argument, 0, 'r'},
		{"output-file",     required_argument, 0, 'f'},
		{"color-space",		no_argument, 0, 'c'},
		{"log-gap",			no_argument, 0, 'L'},
		{"none-stop",		no_argument, 0, 'N'},
		{"bam",				no_argument, 0, 'b'},
		{"single",			no_argument, 0, '0'},
		{"read1",			no_argument, 0, '1'},
		{"read2",			no_argument, 0, '2'},
		{"verbose",			no_argument, 0, 'v'},
		{0, 0, 0, 0}
	};
*/
	static struct option long_options [] = {

			/*input options*/
			{"color-space",		no_argument, 		0, 'c'},
			{"bam",				no_argument, 		0, 0},
			{"single",			no_argument, 		0, '0'},
			{"read1",			no_argument,		0, '1'},
			{"read2",			no_argument, 		0, '2'},
			{"rg-line",			required_argument, 	0, 0},

			/*basic options*/
			{"output-file",     required_argument, 	0, 'o'},
			{"word-size", 		required_argument, 	0, 'w'},
			{"max-word-occ",	required_argument, 	0, 'W'},
			{"max-total-diff", 	required_argument, 	0, 'M'},
			{"max-word-diff",	required_argument, 	0, 'm'},
			{"max-intron",		required_argument, 	0, 'I'},
			{"min-intron",		required_argument, 	0, 'i'},
			{"min-exon", 		required_argument,  0, 'e'},
			{"min-anchor",		required_argument, 	0, 'a'},
			{"known-junc-min-anchor",   required_argument,      0, 'k'},
			{"regression-model",		required_argument,      0, 'r'},
			{"junction-file",	required_argument, 	0, 'j'},
			{"non-denovo",      no_argument,        0, 'n'},
			{"num-threads",		required_argument, 	0, 't'},
			{"best",		no_argument,		0, 'b'},
			{"single-anchor",	no_argument,		0, 's'},
			{"verbose",			no_argument, 		0, 'v'},

			/*advanced options*/
			{"min-logit-score",         required_argument,      0, 0},
			{"max-overhang",    required_argument,  0, 0},
			{"max-gapo", 		required_argument, 	0, 0},
			{"max-gape", 		required_argument, 	0, 0},
			{"indel-end-skip", 	required_argument, 	0, 0},
			{"gape-max-occ", 	required_argument, 	0, 0},
			{"penalty-mismatch",required_argument, 	0, 0},
			{"penalty-gapo",	required_argument, 	0, 0},
			{"penalty-gape",	required_argument, 	0, 0},
			{"log-gap",			no_argument, 		0, 'L'},
			{"max-entries", 	required_argument, 	0, 0},
			{"repeat",			required_argument, 	0, 0},
			{"none-stop",		no_argument, 		0, 0},

			{0, 0, 0, 0}
		};


	int option_index = 0;

	while ((c = getopt_long (argc, argv,
				"c012o:w:W:M:m:I:i:a:j:r:t:k:vLsn",
				long_options, &option_index)) >= 0) {
		switch (c) {
		case -1 : break; //the end of the options

		case 'c': opt->mode &= ~BWA_MODE_COMPREAD; break;
		case '0': opt->mode |= BWA_MODE_BAM_SE; break;
		case '1': opt->mode |= BWA_MODE_BAM_READ1; break;
		case '2': opt->mode |= BWA_MODE_BAM_READ2; break;
		case 'o': freopen(optarg, "w", stdout); break;
		case 'w': opt->word_size = atoi(optarg); break;
		case 'W': opt->max_word_occ = atoi(optarg); break;
		case 'M':
			if (strstr(optarg, ".")) opt->fnr = atof(optarg), opt->max_diff = -1;
			else opt->max_diff = atoi(optarg), opt->fnr = -1.0;
			break;
		case 'm': opt->max_word_diff = atoi(optarg); break;
		case 'I': opt->max_intron_size = atoi(optarg); break;
		case 'i': opt->min_intron_size = atoi(optarg); break;
		case 'e': opt->min_exon_size = atoi (optarg); break;
		case 'a': opt->min_anchor = atoi(optarg); break;
		case 'k': opt->known_junc_min_anchor = atoi(optarg); break;
		case 's': opt->single_anchor_search = 1; break;
		case 'n': opt->non_denovo_search = 1; break;
		case 'r': opt->regression_file = optarg; break;
		case 'j': opt->junction_file = optarg; break;
		case 't': opt->n_threads = atoi(optarg); break;
		case 'v': opt->verbose = 1; break;


		case 0  : //long options without short representations
			if (long_options[option_index].flag != 0)
				break;

			if (strcmp (long_options[option_index].name, "bam" ) == 0) {
				opt->mode |= BWA_MODE_BAM;
			}
			else if (strcmp (long_options[option_index].name, "rg-line" ) == 0) {
				opt->rg = optarg;
			}
			else if (strcmp (long_options[option_index].name, "max-overhang" ) == 0) {
				opt->max_overhang = atoi(optarg);
			}
			else if (strcmp (long_options[option_index].name, "max-gapo" ) == 0) {
				opt->max_gapo = atoi(optarg);
			}
			else if (strcmp (long_options[option_index].name, "max-gape" ) == 0) {
				opte = atoi(optarg);
			}
			else if (strcmp (long_options[option_index].name, "indel-end-skip" ) == 0) {
				opt->indel_end_skip = atoi(optarg);
			}
			else if (strcmp (long_options[option_index].name, "gape-max-occ" ) == 0) {
				opt->max_del_occ = atoi(optarg);
			}
			else if (strcmp (long_options[option_index].name, "penalty-mismatch" ) == 0) {
				opt->s_mm = atoi(optarg);
			}
			else if (strcmp (long_options[option_index].name, "penalty-gapo" ) == 0) {
				opt->s_gapo = atoi(optarg);
			}
			else if (strcmp (long_options[option_index].name, "penalty-gape" ) == 0) {
				opt->s_gape = atoi(optarg);
			}
			else if (strcmp (long_options[option_index].name, "log-gap" ) == 0) {
				opt->mode |= BWA_MODE_LOGGAP;
			}
			else if (strcmp (long_options[option_index].name, "max-entries" ) == 0) {
				opt->max_entries = atoi(optarg);
			}
			else if (strcmp (long_options[option_index].name, "repeat" ) == 0) {
				opt->max_top2 = atoi(optarg);
			}
			else if (strcmp (long_options[option_index].name, "none-stop" ) == 0) {
				opt->mode |= BWA_MODE_NONSTOP; opt->max_top2 = 0x7fffffff;
			}
			else if (strcmp (long_options[option_index].name, "min-logit-score") == 0 ) {
				opt->min_logit_score = atof (optarg);
			}

			//	int tmp = optarg;
			//}
			break;


		default: return 1;
		}
	}
	if (opte > 0) {
		opt->max_gape = opte;
		opt->mode &= ~BWA_MODE_GAPE;
	}

	if (optind + 2 > argc) {
		fprintf(stderr, "Usage:   olego [options] <prefix> <in.fastx>\n");
		/*
		fprintf(stderr, "[\ninput options]\n");
		fprintf(stderr, " -c,--color-space             input sequences are in the color space\n");
		fprintf(stderr, " --bam                        the input read file is in the BAM format\n");
		fprintf(stderr, " -0,--single                  use single-end reads only (effective with -b)\n");
		fprintf(stderr, " -1,--read1                   use the 1st read in a pair (effective with -b)\n");
		fprintf(stderr, " -2,--read2                   use the 2nd read in a pair (effective with -b)\n");
		fprintf(stderr, " -r,--rg-line          TXT    RG line\n");
		*/
		fprintf(stderr, "Compiled at %s, %s\n", __TIME__, __DATE__);
		fprintf(stderr, "\n[basic options]\n");
		fprintf(stderr, " -o,--output-file      FILE   file to write output to instead of stdout\n");
		fprintf(stderr, " -M,--max-total-diff   INT    max #diff (int)[%d] or missing prob under %.2f err rate (float) [%.2f]\n",
				opt->max_diff, BWA_AVG_ERR, opt->fnr);
		fprintf(stderr, " -w,--word-size        INT    word size to seed alignment of spliced reads [%d]\n", opt->word_size);
		fprintf(stderr, " -W,--max-word-occ     INT    max word occurrence to serve as a seed [%d]\n", opt->max_word_occ);
		fprintf(stderr, " -m,--max-word-diff    INT    max #diff allowed in words [%d]\n", opt->max_word_diff);
		fprintf(stderr, " -I,--max-intron       INT    max intron size for de novo search [%d]\n", opt->max_intron_size);
		fprintf(stderr, " -i,--min-intron       INT    min intron size for de novo search [%d]\n", opt->min_intron_size);
		fprintf(stderr, " -e,--min-exon         INT    minimum exon size [%d]\n", opt->min_exon_size);
		fprintf(stderr, " -a,--min-anchor       INT    min anchor on both sides of an exon junction [%d]\n", opt->min_anchor);
		fprintf(stderr, " -k,--known-junc-min-anchor	INT min anchor on both sides of an known exon junction [%d]\n",opt->known_junc_min_anchor);
		fprintf(stderr, " -b,--best                    only report the best alignments\n");
		fprintf(stderr, " -s,--single-anchor           allow single anchor de-novo junction search\n");
		fprintf(stderr, " -r,--regression-model FILE   the file with the logit regression model\n");
		fprintf(stderr, " -j --junction-file    FILE   exon junction BED file\n");
		fprintf(stderr, " -n,--non-denovo              only search known junctions, must turn on -j\n");
		fprintf(stderr, " -t,--num-threads      INT    number of threads [%d]\n", opt->n_threads);
		fprintf(stderr, " -v,--verbose                 verbose mode\n");

		fprintf(stderr, "\n[advanced options]\n");
		fprintf(stderr, " --min-logit-score     FLOAT  logit score of splice sites motif and intron size, in the range of [0,1) [%.2f]\n", opt->min_logit_score);
		fprintf(stderr, " --max-overhang        INT    maximum number or overhanging nucleotide allowed near junctions [%d]\n", opt->max_overhang);
		fprintf(stderr, " --max-gapo            INT    maximum number or fraction of gap opens [%d]\n", opt->max_gapo);
		fprintf(stderr, " --max-gape            INT    maximum number of gap extensions, -1 for disabling long gaps [-1]\n");
		fprintf(stderr, " --indel-end-skip      INT    do not put an indel within INT bp towards the ends [%d]\n", opt->indel_end_skip);
		fprintf(stderr, " --gape-max-occ        INT    maximum occurrences for extending a long deletion [%d]\n", opt->max_del_occ);
		fprintf(stderr, " --penalty-mismatch    INT    mismatch penalty [%d]\n", opt->s_mm);
		fprintf(stderr, " --penalty-gapo        INT    gap open penalty [%d]\n", opt->s_gapo);
		fprintf(stderr, " --penalty-gape        INT    gap extension penalty [%d]\n", opt->s_gape);
		fprintf(stderr, " --log-gap                    log-scaled gap penalty for long deletions\n");
		fprintf(stderr, " --max-entries         INT    maximum entries in the queue [%d]\n", opt->max_entries);
		fprintf(stderr, " --repeat              INT    stop searching when there are >INT equally best hits [%d]\n", opt->max_top2);
		fprintf(stderr, " --none-stop                  non-iterative mode: search for all n-difference hits (slooow)\n");

		//fprintf(stderr, "         -q INT    quality threshold for read trimming down to %dbp [%d]\n", BWA_MIN_RDLEN, opt->trim_qual);
		//fprintf(stderr, "         -l INT    seed length [%d]\n", opt->seed_len);
		//fprintf(stderr, "         -k INT    maximum differences in the seed [%d]\n", opt->max_seed_diff);

		fprintf(stderr, "\n");
		return 1;
	}
	if (opt->fnr > 0.0) {
		int i, k;
		for (i = 17, k = 0; i <= 250; ++i) {
			int l = bwa_cal_maxdiff(i, BWA_AVG_ERR, opt->fnr);
			if (l != k) fprintf(stderr, "[olego_aln] %dbp reads: max_diff = %d\n", i, l);
			k = l;
		}
	}
	if (opt->non_denovo_search && (opt->splice_site_map == 0)) {
	    fprintf(stderr,"[olego_aln] Warning! Non-denovo search mode without junction annotation, no splice mapping will be reported!\n" );
	}
	jigsaw_aln_core(argv[optind], argv[optind+1], opt);
	free(opt);
	return 0;
}

