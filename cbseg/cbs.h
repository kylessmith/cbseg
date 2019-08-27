#ifndef __CBS_H__
#define __CBS_H__

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include "augmented_interval_list.h"

//-------------------------------------------------------------------------------------

typedef enum { false, true } bool;

typedef struct {
    int32_t start;
    int32_t end;
    double_t t;
} segment_t;

typedef struct {
    int32_t index;
    double_t value;
} index_pair_t;

typedef struct {
    int32_t max_start;
    int32_t max_end;
    double_t max_t;
    bool threshold;
} cbs_stat_t;

//-------------------------------------------------------------------------------------

segment_t segment_init(int32_t start, int32_t end, double_t t);

index_pair_t index_pair_init(int32_t index, double_t value);

cbs_stat_t cbs_stat_init(int32_t max_start, int32_t max_end, double_t max_t, bool threshold);

cbs_stat_t calculate_cbs_stat(double x[], int length);

double calculate_tstat(double x[], int i, int length);

cbs_stat_t calculate_cbs(double x[], int length, int shuffles, double p);

void rsegment(double x[], int start, int end, ailist_t *L, int shuffles, double p);

ailist_t *calculate_segment(double x[], int length, int shuffles, double p);

ailist_t *calculate_validate(double x[], int length,  ailist_t *L, int shuffles, double p);

#endif