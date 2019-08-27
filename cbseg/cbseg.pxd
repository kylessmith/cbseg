import numpy as np
cimport numpy as np
cimport cython
from libc.stdint cimport uint32_t, int32_t, int64_t
from ailist cimport AIList, ailist_t


cdef extern from "cbs.c":
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


    segment_t segment_init(int32_t start, int32_t end, double t) nogil

    index_pair_t index_pair_init(int32_t index, double value) nogil

    cbs_stat_t cbs_stat_init(int32_t max_start, int32_t max_end, double max_t) nogil

    cbs_stat_t calculate_cbs_stat(double x[], int length) nogil

    double calculate_tstat(double x[], int i, int length) nogil

    cbs_stat_t calculate_cbs(double x[], int length, int shuffles, double p) nogil

    void rsegment(double x[], int start, int end, ailist_t *L, int shuffles, double p) nogil

    ailist_t *calculate_segment(double x[], int length, int shuffles, double p) nogil

    ailist_t *calculate_validate(double x[], int length,  ailist_t *L, int shuffles, double p) nogil


cdef cbs_stat_t _determine_cbs_stat(double[::1] x)

cdef double _determine_t_stat(double[::1] x, int i)

cdef cbs_stat_t _determine_cbs(double[::1] x, int shuffles, double p)

cdef ailist_t *_segment(double[::1] x, int shuffles, double p)

cdef ailist_t *_validate(double[::1] x, ailist_t *L, int shuffles, double p)