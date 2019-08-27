//=============================================================================
// Circular Binary Segmentation
// by Kyle S. Smith
//-----------------------------------------------------------------------------

#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include "cbs.h"
#include "augmented_interval_list.h"
#include "utilities.c"


//-----------------------------------------------------------------------------

segment_t segment_init(int32_t start, int32_t end, double_t t)
{
    // Initialize interval
    segment_t i;

    // Set attributes
    i.start = start;
    i.end = end;
    i.t = t;

    return i;
}


index_pair_t index_pair_init(int32_t index, double_t value)
{
    // Initialize index_pair
    index_pair_t ip;

    // Set attributes
    ip.index = index;
    ip.value = value;

    return ip;
}


cbs_stat_t cbs_stat_init(int32_t max_start, int32_t max_end, double_t max_t, bool threshold)
{
    // Initialize cbs_stat
    cbs_stat_t c;

    // Set attributes
    c.max_start = max_start;
    c.max_end = max_end;
    c.max_t = max_t;
    c.threshold = threshold;

    return c;
}

//-----------------------------------------------------------------------------

cbs_stat_t calculate_cbs_stat(double x[], int length)
{
    cbs_stat_t cbs_stat = cbs_stat_init(0,0,0, false);
    
    double mean = calculate_mean(x, length);
    double sum = 0;
    index_pair_t i_max = index_pair_init(-1, DBL_MIN);
    index_pair_t i_min = index_pair_init(-1, DBL_MAX);
    index_pair_t e0;
    double x0;

    int i;
    for (i = 0; i < length; i++)
    {
        x0 = (x[i] - mean);
        sum += x0;
        //y[i] = sum;

        e0 = index_pair_init(i, sum);
        i_max = max_pair_value(i_max, e0);
        i_min = min_pair_value(i_min, e0);
    }

    index_pair_t i0 = min_pair_index(i_max, i_min);
    index_pair_t i1 = max_pair_index(i_max, i_min);

    cbs_stat.max_start = i0.index;
    cbs_stat.max_end = i1.index + 1;
    cbs_stat.max_t = pow((i1.value - i0.value), 2.0) * length;
    cbs_stat.max_t = cbs_stat.max_t / (i1.index - i0.index);
    cbs_stat.max_t = cbs_stat.max_t / (length + 1 - i1.index + i0.index);

    return cbs_stat;
}


double calculate_tstat(double x[], int i, int length)
{
    double s0 = 0;
    double s1 = 0;

    int j;
    for (j = 0; j < length; j++)
    {
        if (j < i)
        {
            s0 += x[j];
        }
        else
        {
            s1 += x[j];
        }
    }
    // Find means
    s0 = s0 / i;
    s1 = s1 / (length - i);

    return (length - i) * i / length * square((s0 - s1));
}

//-----------------------------------------------------------------------------

cbs_stat_t calculate_cbs(double x[], int length, int shuffles, double p)
{
    cbs_stat_t cbs_stat = calculate_cbs_stat(x, length);

    if ((cbs_stat.max_end - cbs_stat.max_start) == length)
    {
        cbs_stat.threshold = false;
        return cbs_stat;
    }

    if (cbs_stat.max_start < 5)
    {
        cbs_stat.max_start = 0;
    }

    if ((length - cbs_stat.max_end) < 5)
    {
        cbs_stat.max_end = length;
    }

    int thresh_count = 0;
    double alpha = shuffles * p;
    //double tmp_xt[length];
    double *tmp_xt = malloc(length * sizeof(double));
    double *xt = copy_array(tmp_xt, x, length);
    cbs_stat_t cbs_stat_shuffled;

    int i;
    for (i = 0; i < shuffles; i++)
    {
        shuffle_array(xt, length);
        cbs_stat_shuffled = calculate_cbs_stat(xt, length);

        if (cbs_stat_shuffled.max_t >= cbs_stat.max_t)
        {
            thresh_count += 1;
        }

        if (thresh_count > alpha)
        {   
            free(tmp_xt);
            cbs_stat.threshold = false;
            return cbs_stat;
        }
    }

    cbs_stat.threshold = true;
    free(tmp_xt);

    return cbs_stat;
}


void rsegment(double x[], int start, int end, ailist_t *L, int shuffles, double p)
{   
    int slice_length = end - start;
    //double tmp_slice_x[slice_length];
    double *tmp_slice_x = malloc(slice_length * sizeof(double));
    double *slice_x = slice_array(tmp_slice_x, x, start, end);
    cbs_stat_t cbs_stat = calculate_cbs(slice_x, slice_length, shuffles, p);

    free(tmp_slice_x);

    int cbs_length = cbs_stat.max_end - cbs_stat.max_start;

    if (cbs_stat.threshold == false || (cbs_length < 5) || (cbs_length == slice_length))
    {
        ailist_add(L, start, end, 0, cbs_stat.max_t);
    }
    else
    {
        if (cbs_stat.max_start > 0)
        {
            rsegment(x, start, start + cbs_stat.max_start, L, shuffles, p);
        }
        
        if (cbs_length > 0)
        {
            rsegment(x, start + cbs_stat.max_start, start + cbs_stat.max_end, L, shuffles, p);
        }
        
        if (start + cbs_stat.max_end < end)
        {
            rsegment(x, start + cbs_stat.max_end, end, L, shuffles, p);
        }
    }
}


ailist_t *calculate_segment(double x[], int length, int shuffles, double p)
{
    int start = 0;
    int end = length;
    ailist_t *L = ailist_init();

    srand(time(NULL));

    rsegment(x, start, end, L, shuffles, p);

    return L;
}


ailist_t *calculate_validate(double x[], int length,  ailist_t *L, int shuffles, double p)
{
    int i;
    int S[L->nr + 1];
    for (i = 0; i < L->nr; i++)
    {
        S[i] = L->interval_list[i].start;
    }
    S[i] = length;

    ailist_t *SV = ailist_init();
    int segment_start = 0;
    int left = 0;
    double t;

    int test;
    for (test = 0; test < L->nr-1; test++)
    { 
        int slice_length = S[test+2] - S[left];
        double *tmp_slice_x = malloc(slice_length * sizeof(double));
        double *slice_x = slice_array(tmp_slice_x, x, S[left], S[test+2]);

        t = calculate_tstat(slice_x, S[test+1]-S[left], slice_length);

        double threshold = 0;
        int thresh_count = 0;
        int site = S[test+1] - S[left];

        //is this needed?
        double *tmp_xt = malloc(slice_length * sizeof(double));
        double *xt = copy_array(tmp_xt, slice_x, slice_length);
        bool flag = true;

        free(tmp_slice_x);

        int k;
        for (k = 0; k < shuffles; k++)
        {
            shuffle_array(xt, slice_length);
            threshold = calculate_tstat(xt, site, slice_length);

            if (threshold > t)
            {
                thresh_count += 1;
            }

            if (thresh_count >= (p * shuffles))
            {
                flag = false;
                break;
            }
        }

        if (flag == true)
        {
            ailist_add(SV, segment_start, S[test+1], 0, t);
            segment_start = S[test+1];
            left += 1;
        }

        free(tmp_xt);
    }
    ailist_add(SV, segment_start, S[i], 0, t);

    return SV;
}