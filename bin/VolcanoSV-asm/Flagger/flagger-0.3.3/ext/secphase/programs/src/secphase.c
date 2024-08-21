#include <getopt.h>
#include "sam.h"
#include "faidx.h"
#include <time.h>
#include "bgzf.h"
#include <regex.h>
#include "sonLib.h"
#include <assert.h>
#include <math.h>
#include <float.h>
#include "cigar_it.h"
#include "common.h"
#include "vcf.h"
#include "edlib.h"
#include <time.h>
#include <string.h>
#include "ptBlock.h"
#include "ptVariant.h"
#include "ptAlignment.h"
#include "ptMarker.h"
#include "tpool.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#define SCORE_TYPE_MARKER           0
#define SCORE_TYPE_EDIT_DISTANCE    1

pthread_mutex_t mutex;


void print_alignment_scores(ptAlignment **alignments, int alignments_len, int best_idx, int score_type,
                            FILE *output_log_file) {
    for (int i = 0; i < alignments_len; i++) {
        // change primary to secondary
        if ((alignments[i]->record->core.flag & BAM_FSECONDARY) == 0) {
            fprintf(output_log_file, "*\t");
            //alignments[i]->record->core.flag |= BAM_FSECONDARY;
        }
            //change secondary to primary for the best alignment
        else if (i == best_idx) {
            fprintf(output_log_file, "@\t");
            //alignments[i]->record->core.flag &= ~BAM_FSECONDARY;
        } else fprintf(output_log_file, "!\t");
        if (score_type == SCORE_TYPE_EDIT_DISTANCE) {
            int edit_distance = -1 * alignments[i]->score;
            fprintf(output_log_file, "%d\t%s\t%ld\t%d\n", edit_distance, alignments[i]->contig,
                    alignments[i]->record->core.pos,
                    alignments[i]->rfe);
        } else if (score_type == SCORE_TYPE_MARKER) {
            fprintf(output_log_file, "%.2f\t%s\t%ld\t%d\n", alignments[i]->score, alignments[i]->contig,
                    alignments[i]->record->core.pos, alignments[i]->rfe);
        }
    }
    fprintf(output_log_file, "\n");
    fflush(output_log_file);
}

void merge_and_save_blocks(stHash *blocks_per_contig, char *info_str, char *bed_path, bool print_count_data) {

    // sort and merge modified blocks
    ptBlock_sort_stHash_by_rfs(blocks_per_contig); // sort in place
    stHash *merged_blocks_per_contig = ptBlock_merge_blocks_per_contig_by_rf_v2(blocks_per_contig);
    fprintf(stderr, "[%s] Total length of %s: %d.\n", get_timestamp(), info_str,
            ptBlock_get_total_length_by_rf(merged_blocks_per_contig));
    fprintf(stderr, "[%s] Total number of %s: %d.\n", get_timestamp(), info_str,
            ptBlock_get_total_number(merged_blocks_per_contig));

    ptBlock_save_in_bed(merged_blocks_per_contig, bed_path, print_count_data);
    fprintf(stderr, "[%s] %s are saved in %s.\n", get_timestamp(), info_str, bed_path);
    stHash_destruct(merged_blocks_per_contig);
}

void *runOneThread(void *arg_) {
    // get all arguments
    work_arg_t *arg = arg_;
    bool baq_flag = arg->baq_flag;
    bool consensus = arg->consensus;
    int indel_threshold = arg->indel_threshold;
    int min_q = arg->min_q;
    int min_score = arg->min_score;
    double prim_margin_score = arg->prim_margin_score;
    double prim_margin_random = arg->prim_margin_random;
    int set_q = arg->set_q;
    double conf_d = arg->conf_d;
    double conf_e = arg->conf_e;
    double conf_b = arg->conf_b;
    samFile *bam_fo = arg->bam_fo;
    int flank_margin = arg->flank_margin;
    bool marker_mode = arg->marker_mode;
    int *reads_modified_by_vars = arg->reads_modified_by_vars;
    int *reads_modified_by_marker = arg->reads_modified_by_marker;
    stHash *variant_ref_blocks_per_contig = arg->variant_ref_blocks_per_contig;
    stHash *modified_blocks_by_vars_per_contig = arg->modified_blocks_by_vars_per_contig;
    stHash *modified_blocks_by_marker_per_contig = arg->modified_blocks_by_marker_per_contig;
    stHash *variant_blocks_all_haps_per_contig = arg->variant_blocks_all_haps_per_contig;
    stHash *marker_blocks_all_haps_per_contig = arg->marker_blocks_all_haps_per_contig;
    pthread_mutex_t *mutexPtr = arg->mutexPtr;
    FILE *output_log_file = arg->output_log_file;

    faidx_t *fai = fai_load(arg->fastaPath);

    int alignments_len = arg->alignments_len;
    ptAlignment** alignments = arg->alignments;
    sam_hdr_t *sam_hdr = arg->sam_hdr;

    int conf_blocks_length;

    // Check if there is any variant block encompassed by any alignment
    // If that is met then select the best alignment based on their edit distances
    // to the variant blocks
    stList *merged_variant_read_blocks = NULL;
    if (overlap_variant_ref_blocks(variant_ref_blocks_per_contig, alignments, alignments_len)) {
        merged_variant_read_blocks = ptVariant_get_merged_variant_read_blocks(variant_ref_blocks_per_contig,
                                                                              alignments, alignments_len);
        // Set it to NULL if there is no block
        if (stList_length(merged_variant_read_blocks) == 0) {
            stList_destruct(merged_variant_read_blocks);
            merged_variant_read_blocks = NULL;
        }
    }
    if (merged_variant_read_blocks != NULL) {
        stList **variant_blocks_all_haps = set_scores_as_edit_distances(merged_variant_read_blocks,
                                                                        alignments, alignments_len, fai);
        int best_idx = get_best_record_index(alignments, alignments_len, 0, -100, 0);
        bam1_t *best = 0 <= best_idx ? alignments[best_idx]->record : NULL;
        if (best && (best->core.flag & BAM_FSECONDARY)) {
            //lock mutex
            pthread_mutex_lock(mutexPtr);
            fprintf(output_log_file, "#EDIT DISTANCE\n");
            fprintf(output_log_file, "$\t%s\n", bam_get_qname(alignments[0]->record));
            print_alignment_scores(alignments, alignments_len, best_idx, SCORE_TYPE_EDIT_DISTANCE,
                                   output_log_file);
            // add modified blocks
            int primary_idx = get_primary_index(alignments, alignments_len);
            ptBlock_add_alignment(modified_blocks_by_vars_per_contig, alignments[primary_idx], true);
            ptBlock_add_alignment(modified_blocks_by_vars_per_contig, alignments[best_idx], true);
            // add variant blocks
            ptBlock_add_blocks_by_contig(variant_blocks_all_haps_per_contig,
                                         alignments[primary_idx]->contig,
                                         variant_blocks_all_haps[primary_idx]);
            ptBlock_add_blocks_by_contig(variant_blocks_all_haps_per_contig,
                                         alignments[best_idx]->contig,
                                         variant_blocks_all_haps[best_idx]);
            *reads_modified_by_vars += 1;
            //unlock mutex
            pthread_mutex_unlock(mutexPtr);
        }
        stList_destruct(merged_variant_read_blocks);
        for (int i = 0; i < alignments_len; i++) {
            stList_destruct(variant_blocks_all_haps[i]);
        }
        free(variant_blocks_all_haps);
    }
        // If there is no overlap with variant blocks go to the marker consistency mode
    else if (marker_mode) {
        stList *markers = ptMarker_get_initial_markers(alignments, alignments_len, min_q);
        remove_all_mismatch_markers(&markers, alignments_len);
        sort_and_fill_markers(&markers, alignments, alignments_len);
        filter_ins_markers(&markers, alignments, alignments_len);
        if (markers && stList_length(markers) > 0) {
            int flank_margin_eff = flank_margin;
            set_confident_blocks(alignments, alignments_len, indel_threshold);
            while (consensus && needs_to_find_blocks(alignments, alignments_len, 1000, sam_hdr)) {
                flank_margin_eff *= 0.8;
                set_flanking_blocks(alignments, alignments_len, markers, flank_margin_eff);
                conf_blocks_length = correct_conf_blocks(alignments, alignments_len, indel_threshold);
                if (conf_blocks_length == 0) break;
            }
            if (conf_blocks_length > 0 || consensus == false) {
                if (baq_flag) {
                    calc_update_baq_all(fai,
                                        alignments, alignments_len,
                                        markers, sam_hdr,
                                        conf_d, conf_e, conf_b, set_q);
                }
                filter_lowq_markers(&markers, min_q);
                calc_alignment_score(markers, alignments);
            }
        }

        if (bam_fo != NULL) {
            pthread_mutex_lock(mutexPtr);
            for (int i = 0; i < alignments_len; i++) {
                sam_write1(bam_fo, sam_hdr, alignments[i]->record);
            }
            //unlock mutex
            pthread_mutex_unlock(mutexPtr);
        }
        // get the best alignment
        int best_idx = get_best_record_index(alignments, alignments_len, prim_margin_score, min_score,
                                             prim_margin_random);
        bam1_t *best = 0 <= best_idx ? alignments[best_idx]->record : NULL;
        if (best && (best->core.flag & BAM_FSECONDARY)) {
            //lock mutex
            pthread_mutex_lock(mutexPtr);
            fprintf(output_log_file, "#MARKER SCORE\n");
            fprintf(output_log_file, "$\t%s\n", bam_get_qname(alignments[0]->record));
            print_alignment_scores(alignments, alignments_len, best_idx, SCORE_TYPE_MARKER,
                                   output_log_file);
            // add modified blocks based on modified read coordinates
            int primary_idx = get_primary_index(alignments, alignments_len);
            ptBlock_add_alignment(modified_blocks_by_marker_per_contig, alignments[primary_idx], true);
            ptBlock_add_alignment(modified_blocks_by_marker_per_contig, alignments[best_idx], true);
            // add marker blocks
            ptMarker_add_marker_blocks_by_contig(marker_blocks_all_haps_per_contig,
                                                 alignments[primary_idx]->contig,
                                                 primary_idx,
                                                 markers);
            ptMarker_add_marker_blocks_by_contig(marker_blocks_all_haps_per_contig,
                                                 alignments[best_idx]->contig,
                                                 best_idx,
                                                 markers);
            *reads_modified_by_marker += 1;
            //unlock mutex
            pthread_mutex_unlock(mutexPtr);
        }
        stList_destruct(markers);
    }
    // free alignments
    for (int i = 0; i < alignments_len; i++) {
        ptAlignment_destruct(alignments[i]);
        alignments[i] = NULL;
    }
    free(alignments);
    free(arg);
    fai_destroy(fai);
}

void parseAlignmentsAndScatterJobs(void *arg_, int threads) {
    work_arg_t *arg = arg_;

    pthread_mutex_t *mutexPtr = arg->mutexPtr;

    // open input sam/bam file for parsing alignment records
    samFile *fp = sam_open(arg->inputPath, "r");
    sam_hdr_t *sam_hdr = sam_hdr_read(fp);
    bam1_t *b = bam_init1();

    char read_name[100];
    char read_name_new[100];
    memset(read_name, '\0', 100);
    memset(read_name_new, '\0', 100);
    int alignments_len = 0;
    ptAlignment** alignments = (ptAlignment**) malloc(11 * sizeof(ptAlignment*));
    int bytes_read;
    int conf_blocks_length;
    int count_parsed_alignments = 0;
    int count_parsed_reads = 0;
    int alignment_log_idx = 1;
    int alignment_log_size = 20000; // print alignment log every alignment_log_size alignments

    int *reads_modified_by_vars = arg->reads_modified_by_vars;
    int *reads_modified_by_marker = arg->reads_modified_by_marker;

    // create a thread pool
    tpool_t *tm = tpool_create(threads);

    //lock mutex
    pthread_mutex_lock(mutexPtr);
    fprintf(stderr, "[%s] Started parsing alignments\n", get_timestamp());
    pthread_mutex_unlock(mutexPtr);

    int max_remaining_jobs = threads;

    // iterate over alignments
    while (true) {
        bytes_read = sam_read1(fp, sam_hdr, b);
        count_parsed_alignments += 1;

        // if this is not the end of the file
        if (bytes_read > -1) {
            strcpy(read_name_new, bam_get_qname(b));
            if (read_name[0] == '\0') {
                strcpy(read_name, read_name_new);
            }
        }
        // If read name has changed or file is finished
        if ((strcmp(read_name_new, read_name) != 0) || (bytes_read <= -1)) {
            count_parsed_reads += 1;
            // Check if we have more than one alignment
            // and also not too many (more than 10) alignments
            // Secphase currently does not support supplementary alignments
            // Only one primary alignment should exist per read
            if ((alignments_len > 1) &&
                (alignments_len <= 10) &&
                (ptAlignment_supplementary_count(alignments, alignments_len) == 0) &&
                (ptAlignment_primary_count(alignments, alignments_len) == 1)) {
                // add alignments to the work arg
                // these alignments will be destroyed automatically by
                // runOnThread() at the end of the job
                work_arg_t *arg_cp = tpool_shallow_copy_work_arg(arg);
                arg_cp->alignments = alignments;
                arg_cp->alignments_len = alignments_len;
                arg_cp->sam_hdr = sam_hdr;

                // to avoid excessive memory usage limit the number
                // of the jobs waiting in the pool thread
                while(true){
                    if(tm->remaining_cnt < max_remaining_jobs) break;
                }
                // Add a new job to the thread pool
                tpool_add_work(tm, runOneThread, arg_cp);
            }
            else{
                // if the alignments are not passed to runOneThread()
                // they won't be freed so they have to be destroyed here
                for (int i = 0; i < alignments_len; i++) {
                    ptAlignment_destruct(alignments[i]);
                    alignments[i] = NULL;
                }
                free(alignments);
            }
            // initialize for new alignments
            alignments_len = 0;
            alignments = (ptAlignment**) malloc(11 * sizeof(ptAlignment*));
            // update read name
            strcpy(read_name, read_name_new);
        }
        if (count_parsed_alignments >= alignment_log_idx * alignment_log_size) {
            alignment_log_idx += 1;
            //lock mutex
            pthread_mutex_lock(mutexPtr);
            fprintf(stderr,
                    "[%s] #parsed alignments = %d, #parsed reads = %d, #modifed by phased variants = %d, #modifed by markers = %d\n",
                    get_timestamp(),
                    count_parsed_alignments,
                    count_parsed_reads,
                    *reads_modified_by_vars,
                    *reads_modified_by_marker);
            fflush(stderr);
            //lock mutex
            pthread_mutex_unlock(mutexPtr);
        }
        if (bytes_read <= -1) break; // file is finished so break
        if (b->core.flag & BAM_FUNMAP) continue; // unmapped
        if (alignments_len > 10) continue;
        alignments[alignments_len] = ptAlignment_construct(b, sam_hdr);
        alignments_len += 1;
    }

    // wait until all jobs are finished
    tpool_wait(tm);
    // destroy the thread pool
    tpool_destroy(tm);

    // free memory
    sam_hdr_destroy(sam_hdr);
    sam_close(fp);
    bam_destroy1(b);
}





void get_offset_array(char *input_path, int threads, int64_t *start_array, int64_t *end_array) {
    char index_path[1000];
    snprintf(index_path, 1000, "%s.secphase.index", input_path);
    if (file_exists(index_path) == false) {
        fprintf(stderr, "[%s] Error: To use multi threading please run secphase_index before running secphase!\n",
                get_timestamp());
        exit(EXIT_FAILURE);
    }
    FILE *fp = fopen(index_path, "rb");
    // get the total number of addresses in this index file
    int64_t number_of_addresses;
    fread(&number_of_addresses, sizeof(int64_t), 1, fp);
    // read all addresses
    int64_t *all_addresses = malloc(number_of_addresses * sizeof(int64_t));
    fread(all_addresses, sizeof(int64_t), number_of_addresses, fp);
    int64_t step_size = number_of_addresses / threads;

    // fill start and end arrays
    for (int i = 0; i < threads; i++) {
        start_array[i] = all_addresses[i * step_size];
        if (i < (threads - 1)) {
            end_array[i] = all_addresses[(i + 1) * step_size];
        }
    }
    // the address of the end of the file
    end_array[threads - 1] = all_addresses[number_of_addresses - 1];
    free(all_addresses);
    fclose(fp);
}

static struct option long_options[] =
        {
                {"inputBam",          required_argument, NULL, 'i'},
                {"inputFasta",        required_argument, NULL, 'f'},
                {"inputVcf",          required_argument, NULL, 'v'},
                {"disableMarkerMode", no_argument,       NULL, 'M'},
                {"baq",               no_argument,       NULL, 'q'},
                {"gapOpen",           required_argument, NULL, 'd'},
                {"gapExt",            required_argument, NULL, 'e'},
                {"bandwidth",         required_argument, NULL, 'b'},
                {"consensus",         no_argument,       NULL, 'c'},
                {"indelThreshold",    required_argument, NULL, 't'},
                {"initQ",             required_argument, NULL, 's'},
                {"minQ",              required_argument, NULL, 'm'},
                {"primMarginScore",   required_argument, NULL, 'p'},
                {"primMarginRandom",  required_argument, NULL, 'r'},
                {"minScore",          required_argument, NULL, 'n'},
                {"hifi",              no_argument,       NULL, 'x'},
                {"ont",               no_argument,       NULL, 'y'},
                {"minVariantMargin",  required_argument, NULL, 'g'},
                {"prefix",            required_argument, NULL, 'P'},
                {"outDir",            required_argument, NULL, 'o'},
                {"variantBed",        required_argument, NULL, 'B'},
                {"minGQ",             required_argument, NULL, 'G'},
                {"threads",           required_argument, NULL, '@'},
                {"writeBam",          no_argument,       NULL, 'w'},
                {"flankMargin",       required_argument, NULL, 'F'},
                {NULL,                0,                 NULL, 0}
        };


int main(int argc, char *argv[]) {
    int c;
    bool baq_flag = false;
    bool consensus = false;
    int indel_threshold = 10;
    int min_q = 10;
    int min_score = -10;
    double prim_margin_score = 40;
    double prim_margin_random = 0;
    int set_q = 40;
    int min_var_margin = 50;
    int min_gq = 10;
    double conf_d = 1e-4;
    double conf_e = 0.1;
    double conf_b = 20;
    char inputPath[1000];
    char fastaPath[1000];
    char variantBedPath[1000];
    variantBedPath[0] = NULL;
    char vcfPath[1000];
    vcfPath[0] = NULL;
    char prefix[1000];
    strcpy(prefix, "secphase");
    char dirPath[1000];
    strcpy(dirPath, "secphase_out_dir");
    char *program;
    bool preset_ont = false;
    bool preset_hifi = false;
    bool marker_mode = true;
    bool write_bam = false;
    int flank_margin = 500;
    int threads = 4;
    (program = strrchr(argv[0], '/')) ? ++program : (program = argv[0]);
    while (~(c = getopt_long(argc, argv, "i:p:P:G:o:f:v:qd:e:b:n:r:m:ct:s:B:g:@:wxyMh", long_options, NULL))) {
        switch (c) {
            case 'i':
                strcpy(inputPath, optarg);
                break;
            case 'f':
                strcpy(fastaPath, optarg);
                break;
            case '@':
                threads = atoi(optarg);
                break;
            case 'w':
                write_bam = true;
                break;
            case 'v':
                strcpy(vcfPath, optarg);
                break;
            case 'P':
                strcpy(prefix, optarg);
                break;
            case 'G':
                min_gq = atoi(optarg);
                break;
            case 'o':
                strcpy(dirPath, optarg);
                break;
            case 'x':
                preset_hifi = true;
                baq_flag = true;
                consensus = true;
                indel_threshold = 10; // indel size threshold
                conf_d = 1e-4;
                conf_e = 0.1;
                conf_b = 20;
                min_q = 10;
                set_q = 40;
                prim_margin_score = 40;
                prim_margin_random = 0;
                min_score = -10;
                break;
            case 'y':
                preset_ont = true;
                baq_flag = true;
                consensus = true;
                indel_threshold = 20; // indel size threshold
                conf_d = 1e-3;
                conf_e = 0.1;
                conf_b = 20;
                min_q = 10;
                set_q = 20;
                prim_margin_score = 20;
                prim_margin_random = 0;
                min_score = -10;
                break;
            case 'q':
                baq_flag = true;
                break;
            case 'd':
                conf_d = atof(optarg);
                break;
            case 'e':
                conf_e = atof(optarg);
                break;
            case 'b':
                conf_b = atof(optarg);
                break;
            case 'c':
                consensus = true;
                break;
            case 't':
                indel_threshold = atoi(optarg);
                break;
            case 's':
                set_q = atoi(optarg);
                break;
            case 'm':
                min_q = atoi(optarg);
                break;
            case 'p':
                prim_margin_score = atof(optarg);
                break;
            case 'r':
                prim_margin_random = atof(optarg);
                break;
            case 'n':
                min_score = atoi(optarg);
                break;
            case 'g':
                min_var_margin = atoi(optarg);
                break;
            case 'B':
                strcpy(variantBedPath, optarg);
                break;
            case 'M':
                marker_mode = false;
                break;
            case 'F':
                flank_margin = atoi(optarg);
                break;
            default:
                if (c != 'h') fprintf(stderr, "[E::%s] undefined option %c\n", __func__, c);
            help:
                fprintf(stderr, "\nUsage: %s  -i <INPUT_BAM> -f <FASTA> \n", program);
                fprintf(stderr, "Options:\n");
                fprintf(stderr, "         --inputBam, -i         Input BAM file\n");
                fprintf(stderr, "         --inputFasta, -f         Input FASTA file\n");
                fprintf(stderr, "         --inputVcf, -v         Input phased VCF file\n");
                fprintf(stderr,
                        "         --variantBed, -B         Input BED file for subsetting phased variants\n");
                fprintf(stderr,
                        "         --outDir, -o         Output dir for saving outputs [Default = \"secphase_out_dir\"]\n");
                fprintf(stderr,
                        "         --prefix, -P         Prefix of the output files [Default = \"secphase\"]\n");
                fprintf(stderr,
                        "         --disableMarkerMode, -M         If alignments do not overlap with variants Secphase will not switch to marker mode\n");
                fprintf(stderr,
                        "         --hifi, -x         hifi preset params (only for marker mode) [-q -c -t10 -d 1e-4 -e 0.1 -b20 -m10 -s40 -p40 -r0 -n -10] (Only one of --hifi or --ont should be enabled)\n");
                fprintf(stderr,
                        "         --ont, -y        ont preset params (only for marker mode) [-q -c -t20 -d 1e-3 -e 0.1 -b20 -m10 -s20 -p20 -r0 -n -10] (Only one of --hifi or --ont should be enabled) \n");
                fprintf(stderr, "         --baq, -q         Calculate BAQ [Disabled by default]\n");
                fprintf(stderr, "         --gapOpen, -d         Gap prob [Default: 1e-4, (for ONT use 1e-2)]\n");
                fprintf(stderr, "         --gapExt, -e         Gap extension [Default: 0.1]\n");
                fprintf(stderr, "         --bandwidth, -b         DP bandwidth [Default: 20]\n");
                fprintf(stderr,
                        "         --consensus, -c         Use consensus confident blocks [Disabled by default]\n");
                fprintf(stderr,
                        "         --indelThreshold, -t         Indel size threshold for confident blocks [Default: 10 (for ONT use 20)]\n");
                fprintf(stderr,
                        "         --initQ, -s         Before calculating BAQ set all base qualities to this number [Default: 40 (for ONT use 20)]\n");
                fprintf(stderr,
                        "         --minQ, -m         Minimum base quality (or BAQ if -q is set) to be considered as a marker  [Default: 20 (for ONT use 10)]\n");
                fprintf(stderr,
                        "         --primMarginScore, -p         Minimum margin between the consistency score of primary and secondary alignment to select the secondary alignment [Default: 40]\n");
                fprintf(stderr,
                        "         --primMarginRandom, -r         Maximum margin between the consistency score of primary and secondary alignment to select one randomly [Default: 0]\n");

                fprintf(stderr,
                        "         --minScore, -n         Minimum marker score of the selected secondary alignment [Default: -10]\n");
                fprintf(stderr,
                        "         --minVariantMargin, -g         Minimum margin for creating blocks around phased variants [Default: 50]\n");
                fprintf(stderr,
                        "         --minGQ, -G         Minimum genotype quality of the phased variants [Default: 10]\n");
                fprintf(stderr,
                        "         --writeBam, -w         Write an output bam file with the base qualities modified by BAQ\n");
                fprintf(stderr,
                        "         --flankMargin, -F         Write an output bam file with the base qualities modified by BAQ\n");
                fprintf(stderr,
                        "         --threads, -@         Number of threads  (For using more than one threads secphase_index should be run before running secphase) [Default: 4]\n");
                return 1;
        }
    }


    // make sure the given out directory does not end with '/'
    if (dirPath[strlen(dirPath) - 1] == '/') {
        dirPath[strlen(dirPath) - 1] == '\0';
    }
    struct stat st = {0};

    // make directory for saving output files
    if (stat(dirPath, &st) == -1) {
        mkdir(dirPath, 0777);
    }

    faidx_t *fai = fai_load(fastaPath);

    stHash *variant_ref_blocks_per_contig;
    char bed_path_ref_blocks[1000];
    snprintf(bed_path_ref_blocks, 1000, "%s/%s.initial_variant_blocks.bed", dirPath, prefix);
    if (vcfPath != NULL && vcfPath[0] != NULL) {
        variant_ref_blocks_per_contig = ptVariant_parse_variants_and_extract_blocks(vcfPath, variantBedPath,
                                                                                    fai,
                                                                                    min_var_margin,
                                                                                    min_gq);
    } else {
        // if no vcf is given just make an empty table
        variant_ref_blocks_per_contig = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL,
                                                          (void (*)(void *)) stList_destruct);
    }
    ptVariant_save_variant_ref_blocks(variant_ref_blocks_per_contig, bed_path_ref_blocks);


    if (preset_ont && preset_hifi) {
        fprintf(stderr, "[%s] Presets --hifi and --ont cannot be enabled at the same time. Select only one of them!\n",
                get_timestamp());
        exit(EXIT_FAILURE);
    }

    samFile *fp = NULL;
    sam_hdr_t *sam_hdr = NULL;
    samFile *bam_fo = NULL;

    if (write_bam) {
        // read input file for getting the header
        fp = sam_open(inputPath, "r");
        sam_hdr = sam_hdr_read(fp);

        // open a bam file for writing
        char output_bam_path[1000];
        snprintf(output_bam_path, 1000, "%s/%s.quality_modified.out.bam", dirPath, prefix);
        bam_fo = sam_open(output_bam_path, "w");
        sam_hdr_write(bam_fo, sam_hdr);

        // close input bam
        sam_hdr_destroy(sam_hdr);
        sam_close(fp);
    }

    stHash *modified_blocks_by_vars_per_contig = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL,
                                                                   (void (*)(void *)) stList_destruct);

    stHash *modified_blocks_by_marker_per_contig = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL,
                                                                     (void (*)(void *)) stList_destruct);

    stHash *variant_blocks_all_haps_per_contig = stHash_construct3(stHash_stringKey, stHash_stringEqualKey,
                                                                   NULL,
                                                                   (void (*)(void *)) stList_destruct);
    stHash *marker_blocks_all_haps_per_contig = stHash_construct3(stHash_stringKey, stHash_stringEqualKey,
                                                                  NULL,
                                                                  (void (*)(void *)) stList_destruct);


    pthread_mutex_t *mutexPtr = malloc(sizeof(pthread_mutex_t));
    pthread_mutex_init(mutexPtr, NULL);

    //int64_t *bam_adr_start_array = malloc(threads * sizeof(int64_t));
    //int64_t *bam_adr_end_array = malloc(threads * sizeof(int64_t));
    //get_offset_array(inputPath, threads, bam_adr_start_array, bam_adr_end_array);

    // open file for saving reads that have to be corrected
    char output_log_path[1000];
    snprintf(output_log_path, 1000, "%s/%s.out.log", dirPath, prefix);
    FILE *output_log_file = fopen(output_log_path, "w+");

    int reads_modified_by_vars = 0;
    int reads_modified_by_marker = 0;

    work_arg_t *arg_temp =
            tpool_create_work_arg(baq_flag, consensus, indel_threshold, min_q,
                                  min_score, prim_margin_score, prim_margin_random, set_q,
                                  conf_d, conf_e, conf_b, inputPath,
                                  fastaPath,
                                  marker_mode,
                                  variant_ref_blocks_per_contig,
                                  modified_blocks_by_vars_per_contig,
                                  modified_blocks_by_marker_per_contig,
                                  variant_blocks_all_haps_per_contig,
                                  marker_blocks_all_haps_per_contig,
                                  &reads_modified_by_vars,
                                  &reads_modified_by_marker,
                                  output_log_file, bam_fo,
                                  mutexPtr, NULL, 0, flank_margin, NULL);

    parseAlignmentsAndScatterJobs(arg_temp, threads);

    fprintf(stderr, "[%s] Number of reads modified by phased variants = %d\n", get_timestamp(),
            reads_modified_by_vars);
    fprintf(stderr, "[%s] Number of reads modified by marker score = %d\n", get_timestamp(),
            reads_modified_by_marker);


    // save the read blocks modified by markers and variants
    char bed_path[1000];

    snprintf(bed_path, 1000, "%s/%s.modified_read_blocks.variants.bed", dirPath, prefix);
    merge_and_save_blocks(modified_blocks_by_vars_per_contig, "read blocks modified by phased variants",
                          bed_path, true);

    snprintf(bed_path, 1000, "%s/%s.modified_read_blocks.markers.bed", dirPath, prefix);
    merge_and_save_blocks(modified_blocks_by_marker_per_contig, "read blocks modified by markers",
                          bed_path, true);

    // save the variant and marker blocks
    snprintf(bed_path, 1000, "%s/%s.variant_blocks.bed", dirPath, prefix);
    merge_and_save_blocks(variant_blocks_all_haps_per_contig,
                          "projected variant blocks on all haplotypes",
                          bed_path, false);

    snprintf(bed_path, 1000, "%s/%s.marker_blocks.bed", dirPath, prefix);
    merge_and_save_blocks(marker_blocks_all_haps_per_contig,
                          "projected marker blocks on all haplotypes",
                          bed_path,false);


    fclose(output_log_file);
    stHash_destruct(variant_ref_blocks_per_contig);
    stHash_destruct(modified_blocks_by_marker_per_contig);
    stHash_destruct(modified_blocks_by_vars_per_contig);
    stHash_destruct(marker_blocks_all_haps_per_contig);
    stHash_destruct(variant_blocks_all_haps_per_contig);

    if (write_bam) sam_close(bam_fo);
}
//main();
