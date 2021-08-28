import numpy as np
cimport numpy as np
cimport cython
from libc.stdint cimport uint32_t, int32_t, int64_t, uint16_t
from ailist.IntervalArray_core cimport IntervalArray, aiarray_t
from ailist.LabeledIntervalArray_core cimport LabeledIntervalArray, labeled_aiarray_t


cdef extern from "cbs.c":
    # C is include here so that it doesn't need to be compiled externally
    pass

cdef extern from "cbs_labeled.c":
    # C is include here so that it doesn't need to be compiled externally
    pass

cdef extern from "cbs_structs.c":
    # C is include here so that it doesn't need to be compiled externally
    pass

cdef extern from "utilities.c":
    # C is include here so that it doesn't need to be compiled externally
    pass

cdef extern from "cbs.h":
    ctypedef struct segment_t:
        int32_t start
        int32_t end
        double t

    ctypedef struct index_pair_t:
        int32_t index
        double value

    ctypedef struct cbs_stat_t:
        int32_t max_start
        int32_t max_end
        double max_t


    #-------------------------------------------------------------------------------------
    #                                    cbs_structs.c
    #=====================================================================================

    #segment_t segment_init(int32_t start, int32_t end, double_t t) nogil

    #index_pair_t index_pair_init(int32_t index, double_t value) nogil

    #cbs_stat_t cbs_stat_init(int32_t max_start, int32_t max_end, double_t max_t, bool threshold) nogil

    cbs_stat_t calculate_cbs_stat(double x[], int length) nogil

    double calculate_tstat(double x[], int i, int length) nogil

    cbs_stat_t calculate_cbs(double x[], int length, int shuffles, double p) nogil


    #-------------------------------------------------------------------------------------
    #                                    cbs.c
    #=====================================================================================

    void rsegment(double x[], int start, int end, aiarray_t *L, int shuffles, double p) nogil

    aiarray_t *calculate_segment(double x[], int length, int shuffles, double p) nogil

    aiarray_t *calculate_validate(double x[], int length,  aiarray_t *L, int shuffles, double p) nogil


    #-------------------------------------------------------------------------------------
    #                                    cbs_labeled.c
    #=====================================================================================

    void rsegment_labeled(double x[], int start, int end, labeled_aiarray_t *segments, char *label, int shuffles, double p) nogil

    labeled_aiarray_t *calculate_segment_labeled(double x[], int length, labeled_aiarray_t *segments, char *label, int shuffles, double p) nogil

    labeled_aiarray_t *calculate_validate_labeled(double x[], int length,  labeled_aiarray_t *segments, labeled_aiarray_t *vsegments, char *label, int shuffles, double p) nogil


# Functions
cdef cbs_stat_t _determine_cbs_stat(double[::1] x)

cdef double _determine_t_stat(double[::1] x, int i)

cdef cbs_stat_t _determine_cbs(double[::1] x, int shuffles, double p)

cdef aiarray_t *_segment(double[::1] x, int shuffles, double p)
cdef void _segment_labeled(double[::1] values, labeled_aiarray_t *c_lsegments, char *label, int shuffles, double p)

cdef aiarray_t *_validate(double[::1] values, IntervalArray segments, int shuffles, double p)
cdef void _validate_labeled(double[::1] values, LabeledIntervalArray segments, labeled_aiarray_t *vsegments, char *label, int shuffles, double p)