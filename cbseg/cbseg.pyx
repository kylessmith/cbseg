import os
import numpy as np
cimport numpy as np
cimport cython
from libc.stdint cimport uint32_t, int32_t, int64_t
from ailist cimport AIList, ailist_t


def get_include():
    """
    Get file directory if C headers
    
    Arguments:
    ---------
        None
    Returns:
    ---------
        str (Directory to header files)
    """

    return os.path.split(os.path.realpath(__file__))[0]


cdef cbs_stat_t _determine_cbs_stat(double[::1] x):
    cdef int length = x.size
    cdef cbs_stat_t cbs_stat = calculate_cbs_stat(&x[0], length)

    return cbs_stat

def determine_cbs_stat(double[::1] x):
    cdef cbs_stat_t cbs_stat = _determine_cbs_stat(x)

    return cbs_stat


cdef double _determine_t_stat(double[::1] x, int i):
    cdef int length = x.size
    cdef double tstat = calculate_tstat(&x[0], i, length)

    return tstat

def determine_t_stat(double[::1] x, int i):
    cdef double tstat = _determine_t_stat(x, i)

    return tstat


cdef cbs_stat_t _determine_cbs(double[::1] x, int shuffles, double p):
    cdef int length = x.size
    cdef cbs_stat_t cbs_stat = calculate_cbs(&x[0], length, shuffles, p)

    return cbs_stat

def determine_cbs(double[::1] x, int shuffles=1000, double p=0.05):
    cdef cbs_stat_t cbs_stat = _determine_cbs(x, shuffles, p)

    return cbs_stat


cdef ailist_t *_segment(double[::1] x, int shuffles, double p):
    cdef int length = x.size
    cdef ailist_t *L = calculate_segment(&x[0], length, shuffles, p)

    return L

def segment(double[::1] x, int shuffles=1000, double p = 0.05):
    cdef ailist_t *L = _segment(x, shuffles, p)
    cdef AIList segments = AIList()
    segments.set_list(L)

    return segments


cdef ailist_t *_validate(double[::1] x, ailist_t *L, int shuffles, double p):
    cdef int length = x.size
    cdef ailist_t *SV = calculate_validate(&x[0], length, L, shuffles, p)

    return SV
    
def validate(double[::1] x, AIList segments, int shuffles=1000, double p = 0.05):
    cdef ailist_t *SV = _validate(x, segments.interval_list, shuffles, p)
    cdef AIList valid_segments = AIList()
    valid_segments.set_list(SV)

    return valid_segments