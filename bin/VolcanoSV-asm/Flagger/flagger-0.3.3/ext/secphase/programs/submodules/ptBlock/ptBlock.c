#include "ptBlock.h"


ptBlock *ptBlock_construct(int rfs, int rfe, int sqs, int sqe, int rds_f, int rde_f) {
    ptBlock *block = malloc(sizeof(ptBlock));
    block->rfs = rfs;
    block->rfe = rfe;
    block->sqs = sqs;
    block->sqe = sqe;
    block->rds_f = rds_f;
    block->rde_f = rde_f;
    block->data = NULL;
    block->destruct_data = free;
    block->copy_data = NULL;
    block->extend_data = NULL;
    return block;
}

ptBlock *ptBlock_construct_with_count(int rfs, int rfe, int sqs, int sqe, int rds_f, int rde_f, int count) {
    ptBlock *block = malloc(sizeof(ptBlock));
    block->rfs = rfs;
    block->rfe = rfe;
    block->sqs = sqs;
    block->sqe = sqe;
    block->rds_f = rds_f;
    block->rde_f = rde_f;
    block->data = NULL;
    block->destruct_data = free;
    block->copy_data = NULL;
    block->extend_data = NULL;
    int* data = malloc(sizeof(int));
    *data = count;
    ptBlock_set_data(block, data, destruct_count_data, copy_count_data, extend_count_data);
    return block;
}

int ptBlock_get_count(ptBlock* block){
    int* count = (int*) block->data;
    return *count;
}

void ptBlock_set_data(ptBlock *block, void *data, void (*destruct_data)(void *), void *(*copy_data)(void *),
                      void (*extend_data)(void *, void *)) {
    if (block->data != NULL) {
        block->destruct_data(block->data);
    }
    block->data = data;
    block->destruct_data = destruct_data;
    block->copy_data = copy_data;
    block->extend_data = extend_data;
}

void ptBlock_extend_data(ptBlock *block, void *data) {
    if (block->extend_data != NULL) {
        block->extend_data(block->data, data);
    }
}

void ptBlock_destruct_data(ptBlock *block) {
    if (block->destruct_data != NULL) {
        block->destruct_data(block->data);
    }
}

void *ptBlock_copy_data(ptBlock *block) {
    if (block->copy_data != NULL) {
        return block->copy_data(block->data);
    } else {
        return NULL;
    }
}


void extend_count_data(void *dest_, void *src_){
    int* dest = dest_;
    int* src = src_;
    *dest += *src;
}


void destruct_count_data(void* src){
    free(src);
}


void *copy_count_data(void* src_){
    int* src = src_;
    int* dest = malloc(sizeof(int));
    *dest = *src;
    return dest;
}

ptBlock *ptBlock_copy(ptBlock *block) {
    ptBlock *block_copy = ptBlock_construct(block->rfs,
                                            block->rfe,
                                            block->sqs,
                                            block->sqe,
                                            block->rds_f,
                                            block->rde_f);
    void *data_copy = ptBlock_copy_data(block);
    ptBlock_set_data(block_copy, data_copy, block->destruct_data, block->copy_data, block->extend_data);
    return block_copy;
}


void ptBlock_destruct(ptBlock *block) {
    ptBlock_destruct_data(block);
    free(block);
}


/// Get Functions ////
int ptBlock_get_rfs(ptBlock *block) {
    return block->rfs;
}


int ptBlock_get_rfe(ptBlock *block) {
    return block->rfe;
}

int ptBlock_get_sqs(ptBlock *block) {
    return block->sqs;
}

int ptBlock_get_sqe(ptBlock *block) {
    return block->sqe;
}

int ptBlock_get_rds_f(ptBlock *block) {
    return block->rds_f;
}

int ptBlock_get_rde_f(ptBlock *block) {
    return block->rde_f;
}


//// Set Functions ////

void ptBlock_set_rfs(ptBlock *block, int rfs) {
    block->rfs = rfs;
}

void ptBlock_set_rfe(ptBlock *block, int rfe) {
    block->rfe = rfe;
}

void ptBlock_set_sqs(ptBlock *block, int sqs) {
    block->sqs = sqs;
}

void ptBlock_set_sqe(ptBlock *block, int sqe) {
    block->sqe = sqe;
}

void ptBlock_set_rds_f(ptBlock *block, int rds_f) {
    block->rds_f = rds_f;
}

void ptBlock_set_rde_f(ptBlock *block, int rde_f) {
    block->rde_f = rde_f;
}

//// Comparing Functions ////
int ptBlock_cmp_rfs(const void *a, const void *b) {
    ptBlock *b1 = (ptBlock *) a;
    ptBlock *b2 = (ptBlock *) b;
    return b1->rfs - b2->rfs;
}

int ptBlock_cmp_rds_f(const void *a, const void *b) {
    ptBlock *b1 = (ptBlock *) a;
    ptBlock *b2 = (ptBlock *) b;
    return b1->rds_f - b2->rds_f;
}

int ptBlock_cmp_sqs(const void *a, const void *b) {
    ptBlock *b1 = (ptBlock *) a;
    ptBlock *b2 = (ptBlock *) b;
    return b1->sqs - b2->sqs;
}


stHash *ptBlock_parse_bed(char *bed_path) {
    FILE *fp = fopen(bed_path, "r");
    size_t read;
    size_t len;
    char *line = NULL;
    char *contig_name;
    int start;
    int end;
    char *token;
    ptBlock *block = NULL;
    stList *blocks = NULL;
    stHash *blocks_per_contig = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL,
                                                  (void (*)(void *)) stList_destruct);
    while ((read = getline(&line, &len, fp)) != -1) {
        // replace '\n' by '\0'
        if (line[strlen(line) - 1] == '\n') {
            line[strlen(line) - 1] = '\0';
        }
        token = strtok(line, "\t");
        contig_name = copyString(token);
        token = strtok(NULL, "\t");
        start = atoi(token); // 0-based
        token = strtok(NULL, "\t");
        end = atoi(token) - 1; // 0-based
        block = ptBlock_construct(start, end, -1, -1, -1, -1);
        blocks = stHash_search(blocks_per_contig, contig_name);
        // if contig does not exist as a key in the table
        // add contig along with an empty list as its value
        if (blocks == NULL) {
            blocks = stList_construct3(0, ptBlock_destruct);
            stHash_insert(blocks_per_contig, contig_name, blocks);
        }
        stList_append(blocks, block);
    }
    // sort blocks per contig
    stHashIterator *it = stHash_getIterator(blocks_per_contig);
    while ((contig_name = stHash_getNext(it)) != NULL) {
        blocks = stHash_search(blocks_per_contig, contig_name);
        stList_sort(blocks, ptBlock_cmp_rfs);
    }
    return blocks_per_contig;
}

void ptBlock_sort_stHash_by_rfs(stHash *blocks_per_contig) {
    char *contig_name;
    stList *blocks;
    stHashIterator *it = stHash_getIterator(blocks_per_contig);
    while ((contig_name = stHash_getNext(it)) != NULL) {
        blocks = stHash_search(blocks_per_contig, contig_name);
        stList_sort(blocks, ptBlock_cmp_rfs);
    }
}

stList *ptBlock_merge_blocks(stList *blocks,
                             int (*get_start)(ptBlock *), int (*get_end)(ptBlock *),
                             void (*set_end)(ptBlock *, int)) {
    stList *blocks_merged = stList_construct3(0, ptBlock_destruct);
    if (stList_length(blocks) == 0) return blocks_merged;
    ptBlock *b = NULL;
    ptBlock *b_merged = NULL;
    for (int i = 0; i < stList_length(blocks); i++) {
        b = stList_get(blocks, i);
        //printf("%d\t%d\n", b->rds_f, b->rde_f);
        if (i == 0) { // Initiate b_merged for the first block
            b_merged = ptBlock_copy(b);
            continue;
        }
        if (get_end(b_merged) < get_start(b)) {// no overlap with previous merged block
            //save the merged block
            stList_append(blocks_merged, b_merged);
            // Initiate a new merged block
            b_merged = ptBlock_copy(b);
        } else { //there is overlap
            set_end(b_merged, max(get_end(b_merged), get_end(b))); // extend end pos of the merged block
            if (b->data != NULL) {
                // add the data of the new block to the merged block
                ptBlock_extend_data(b_merged, b->data);
            }
        }
    }

    // Add the last merged block
    if (b_merged) {
        stList_append(blocks_merged, b_merged);
    }

    return blocks_merged;
}

stList *ptBlock_merge_blocks_v2(stList *blocks,
                             int (*get_start)(ptBlock *), int (*get_end)(ptBlock *),
                             void (*set_start)(ptBlock *, int), void (*set_end)(ptBlock *, int)) {
    stList *blocks_merged_finalized = stList_construct3(0, ptBlock_destruct);
    stList *blocks_merged_ongoing = stList_construct3(0, ptBlock_destruct);
    stList *blocks_merged_temp;
    if (stList_length(blocks) == 0) return blocks_merged_finalized;
    ptBlock *b2 = NULL;
    ptBlock *b1 = NULL;
    ptBlock *b_merged = NULL;
    int s1;
    int e1;
    int s2;
    int e2;
    for (int i = 0; i < stList_length(blocks); i++) {
        b2 = stList_get(blocks, i);
        //printf("%d\t%d\n", b->rds_f, b->rde_f);
        if (stList_length(blocks_merged_ongoing) == 0) { // Initiate b_merged for the first block
            b_merged = ptBlock_copy(b2);
            stList_append(blocks_merged_ongoing, b_merged);
            continue;
        }
        e2 = get_end(b2);
        s2 = get_start(b2);
        blocks_merged_temp = blocks_merged_ongoing;
        blocks_merged_ongoing = stList_construct3(0, ptBlock_destruct);
        for (int j = 0; j < stList_length(blocks_merged_temp); j++) {
            b1 = stList_get(blocks_merged_temp, j);
            e1 = get_end(b1);
            s1 = get_start(b1);
            /*
             * finalized:
             *
             *  s1       e1
             * [**********]
             *               [----------]
             *                s2       e2
             */
            if(e1 < s2) {
                b_merged = ptBlock_copy(b1);
                stList_append(blocks_merged_finalized, b_merged);
            }
            else if (s1 <= s2){
                /*
                 * finalized:
                 *   s1       e1
                 *  [***-------]
                 *     [----------]
                 *      s2       e2
                 */
                if (s1 < s2){ // && s2 <= e1
                    b_merged = ptBlock_copy(b1);
                    set_end(b_merged, s2 - 1);
                    stList_append(blocks_merged_finalized, b_merged);
                }
                /*
                 *  ongoing:
                 *
                 *    s1       e1                      s1        e1
                 *   [---*******]           OR        [----*****--]
                 *      [*******--]                       [*****]
                 *       s2      e2                        s2   e2
                 *
                 */
                b_merged = ptBlock_copy(b1);
                set_start(b_merged, s2);
                set_end(b_merged, min(e1, e2));
                if (b2->data != NULL) {
                    // add the data of the new block to the merged block
                    ptBlock_extend_data(b_merged, b2->data);
                }
                stList_append(blocks_merged_ongoing, b_merged);

                /*
                 * finalized:
                 *     s1       e1
                 *    [-----*****]
                 *      [---]
                 *       s2 e2
                 */
                if (e2 < e1){
                    b_merged = ptBlock_copy(b1);
                    set_start(b_merged, e2 + 1);
                    stList_append(blocks_merged_ongoing, b_merged);
                }
            }
            /*
             * ongoing:
             *
             *            s1       e1
             *           [**********]
             *       [----**********---]
             *        s2              e2
             */
            else if(e1 <= e2){ // && s2 < s1
                b_merged = ptBlock_copy(b1);
                if (b2->data != NULL) {
                    // add the data of the new block to the merged block
                    ptBlock_extend_data(b_merged, b2->data);
                }
                stList_append(blocks_merged_ongoing, b_merged);
            }
            else { // e2 < e1 && s2 < s1
                /*
                 * ongoing:
                 *
                 *            s1       e1
                 *           [******----]
                 *       [----******]
                 *        s2       e2
                 */
                if(s1 <= e2) {
                    b_merged = ptBlock_copy(b1);
                    set_end(b_merged, e2);
                    if (b2->data != NULL) {
                        // add the data of the new block to the merged block
                        ptBlock_extend_data(b_merged, b2->data);
                    }
                    stList_append(blocks_merged_ongoing, b_merged);
                }
                /*
                * ongoing:
                *
                *            s1        e1
                *           [------****]
                *       [----------]
                *        s2       e2
                */

                b_merged = ptBlock_copy(b1);
                set_start(b_merged, max(e2 + 1, s1));
                stList_append(blocks_merged_ongoing, b_merged);
             }
        }
        // add the last non-overlapping block
        if (max(e1 + 1, s2) <= e2){
            b_merged = ptBlock_copy(b2);
            set_start(b_merged, max(e1 + 1, s2));
            stList_append(blocks_merged_ongoing, b_merged);
        }
        // destruct the temporary blocks
        stList_destruct(blocks_merged_temp);
    }

    // Add the remaining blocks
    for (int j = 0; j < stList_length(blocks_merged_ongoing); j++) {
        b1 = stList_get(blocks_merged_ongoing, j);
        b_merged = ptBlock_copy(b1);
        stList_append(blocks_merged_finalized, b_merged);
    }

    stList_destruct(blocks_merged_ongoing);

    return blocks_merged_finalized;
}

stHash *ptBlock_merge_blocks_per_contig(stHash *blocks_per_contig,
                                        int (*get_start)(ptBlock *), int (*get_end)(ptBlock *),
                                        void (*set_end)(ptBlock *, int)) {
    char *contig_name;
    stList *blocks;
    stList *merged_blocks;
    stHash *merged_blocks_per_contig = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL,
                                                         (void (*)(void *)) stList_destruct);
    stHashIterator *it = stHash_getIterator(blocks_per_contig);
    while ((contig_name = stHash_getNext(it)) != NULL) {
        // get blocks
        blocks = stHash_search(blocks_per_contig, contig_name);
        // merge blocks
        merged_blocks = ptBlock_merge_blocks(blocks, get_start, get_end, set_end);
        // add merged blocks to the new table
        stHash_insert(merged_blocks_per_contig, contig_name, merged_blocks);
    }
    return merged_blocks_per_contig;
}


stHash *ptBlock_merge_blocks_per_contig_by_rf(stHash *blocks_per_contig) {
    return ptBlock_merge_blocks_per_contig(blocks_per_contig, ptBlock_get_rfs,
                                              ptBlock_get_rfe,
                                              ptBlock_set_rfe);
}

stHash *ptBlock_merge_blocks_per_contig_by_rd_f(stHash *blocks_per_contig) {
    return ptBlock_merge_blocks_per_contig(blocks_per_contig, ptBlock_get_rds_f,
                                              ptBlock_get_rde_f,
                                              ptBlock_set_rde_f);
}

stHash *ptBlock_merge_blocks_per_contig_by_sq(stHash *blocks_per_contig) {
    return ptBlock_merge_blocks_per_contig(blocks_per_contig, ptBlock_get_sqs,
                                           ptBlock_get_sqe,
                                           ptBlock_set_sqe);
}

stHash *ptBlock_merge_blocks_per_contig_v2(stHash *blocks_per_contig,
                                        int (*get_start)(ptBlock *), int (*get_end)(ptBlock *),
                                           void (*set_start)(ptBlock *, int), void (*set_end)(ptBlock *, int)) {
    char *contig_name;
    stList *blocks;
    stList *merged_blocks;
    stHash *merged_blocks_per_contig = stHash_construct3(stHash_stringKey, stHash_stringEqualKey, NULL,
                                                         (void (*)(void *)) stList_destruct);
    stHashIterator *it = stHash_getIterator(blocks_per_contig);
    while ((contig_name = stHash_getNext(it)) != NULL) {
        // get blocks
        blocks = stHash_search(blocks_per_contig, contig_name);
        // merge blocks
        merged_blocks = ptBlock_merge_blocks_v2(blocks, get_start, get_end, set_start, set_end);
        // add merged blocks to the new table
        stHash_insert(merged_blocks_per_contig, contig_name, merged_blocks);
    }
    return merged_blocks_per_contig;
}

stHash *ptBlock_merge_blocks_per_contig_by_rf_v2(stHash *blocks_per_contig) {
    return ptBlock_merge_blocks_per_contig_v2(blocks_per_contig, ptBlock_get_rfs,
                                           ptBlock_get_rfe,
                                           ptBlock_set_rfs,
                                           ptBlock_set_rfe);
}

stHash *ptBlock_merge_blocks_per_contig_by_rd_f_v2(stHash *blocks_per_contig) {
    return ptBlock_merge_blocks_per_contig_v2(blocks_per_contig, ptBlock_get_rds_f,
                                           ptBlock_get_rde_f,
                                           ptBlock_set_rds_f,
                                           ptBlock_set_rde_f);
}

stHash *ptBlock_merge_blocks_per_contig_by_sq_v2(stHash *blocks_per_contig) {
    return ptBlock_merge_blocks_per_contig_v2(blocks_per_contig, ptBlock_get_sqs,
                                           ptBlock_get_sqe,
                                           ptBlock_set_sqs,
                                           ptBlock_set_sqe);
}

int ptBlock_get_total_number(stHash *blocks_per_contig) {
    int n = 0;
    char *contig_name;
    stList *blocks;
    stHashIterator *it = stHash_getIterator(blocks_per_contig);
    while ((contig_name = stHash_getNext(it)) != NULL) {
        blocks = stHash_search(blocks_per_contig, contig_name);
        n += stList_length(blocks);
    }
    return n;
}

int ptBlock_get_total_length(stHash *blocks_per_contig, int (*get_start)(ptBlock *), int (*get_end)(ptBlock *)) {
    int total_len = 0;
    char *contig_name;
    stList *blocks;
    ptBlock *block;
    stHashIterator *it = stHash_getIterator(blocks_per_contig);
    while ((contig_name = stHash_getNext(it)) != NULL) {
        blocks = stHash_search(blocks_per_contig, contig_name);
        for (int i = 0; i < stList_length(blocks); i++) {
            block = stList_get(blocks, i);
            total_len += get_end(block) - get_start(block) + 1;
        }
    }
    return total_len;
}

int ptBlock_get_total_length_by_rf(stHash *blocks_per_contig) {
    return ptBlock_get_total_length(blocks_per_contig, ptBlock_get_rfs, ptBlock_get_rfe);
}

int ptBlock_get_total_length_by_rd_f(stHash *blocks_per_contig) {
    return ptBlock_get_total_length(blocks_per_contig, ptBlock_get_rds_f, ptBlock_get_rde_f);
}

int ptBlock_get_total_length_by_sq(stHash *blocks_per_contig) {
    return ptBlock_get_total_length(blocks_per_contig, ptBlock_get_sqs, ptBlock_get_sqe);
}


void ptBlock_add_alignment(stHash *blocks_per_contig, ptAlignment *alignment, bool init_count_data) {
    ptBlock *block = ptBlock_construct(alignment->rfs,
                                       alignment->rfe,
                                       -1, -1,
                                       -1, -1);
    if (init_count_data){
        int* count_data = malloc(sizeof(int));
        *count_data = 1;
        ptBlock_set_data(block, count_data,
                         destruct_count_data,
                         copy_count_data,
                         extend_count_data);
    }
    stList *blocks = stHash_search(blocks_per_contig, alignment->contig);
    if (blocks == NULL) {
        blocks = stList_construct3(0, ptBlock_destruct);
        // contig name should be copied prior to inserting as a key
        stHash_insert(blocks_per_contig, copyString(alignment->contig), blocks);
    }
    stList_append(blocks, block);
}

void ptBlock_save_in_bed(stHash *blocks_per_contig, char* bed_path, bool print_count_data){
    // get contigs as a list and sort them
    stList * contigs = stHash_getKeys(blocks_per_contig);
    stList_sort(contigs, (int (*)(const void *, const void *))strcmp);

    // open bed to write
    FILE *fp = fopen(bed_path, "w");
    if (fp == NULL) {
        fprintf(stderr, "[%s] Error: Failed to open file %s.\n", get_timestamp(), bed_path);
    }
    // iterate over sorted contig names
    for(int i = 0; i < stList_length(contigs); i++){
        char* contig = stList_get(contigs, i);
        stList *blocks = stHash_search(blocks_per_contig, contig);
        // iterate over blocks in this contig
        for(int j = 0; j < stList_length(blocks); j++) {
            ptBlock* block = stList_get(blocks, j);
            // it may happen when the whole projected variant block is within an insertion in another haplotype
            // TODO: experiment if it can help to output these coordinates
            if(block->rfe < block->rfs) continue;
            if(print_count_data) {
                fprintf(fp, "%s\t%d\t%d\t%d\n", contig, block->rfs, block->rfe + 1, *((int*)block->data)); // end should be 1-based in BED
            }else{
                fprintf(fp, "%s\t%d\t%d\n", contig, block->rfs, block->rfe + 1); // end should be 1-based in BED
            }
        }

    }
    fclose(fp);
}

void ptBlock_add_blocks_by_contig(stHash *blocks_per_contig, char* contig, stList *blocks_to_add) {
    stList *blocks = stHash_search(blocks_per_contig, contig);
    if (blocks == NULL) {
        blocks = stList_construct3(0, ptBlock_destruct);
        // contig name should be copied prior to inserting as a key
        stHash_insert(blocks_per_contig, copyString(contig), blocks);
    }
    for (int i =0; i < stList_length(blocks_to_add); i++){
        stList_append(blocks, ptBlock_copy(stList_get(blocks_to_add,i)));
    }
}

stList* ptBlock_copy_stList(stList* blocks) {
    stList *copy_blocks = stList_construct3(0, ptBlock_destruct);
    for (int i = 0; i < stList_length(blocks); i++) {
        ptBlock* block = stList_get(blocks, i);
        stList_append(copy_blocks, ptBlock_copy(block));
    }
    return copy_blocks;
}

