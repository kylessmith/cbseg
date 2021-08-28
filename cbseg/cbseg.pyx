import os
import pandas as pd
import numpy as np
cimport numpy as np
cimport cython
from libc.stdint cimport uint32_t, int32_t, int64_t
from ailist.IntervalArray_core cimport IntervalArray, aiarray_t
from ailist.LabeledIntervalArray_core cimport LabeledIntervalArray, labeled_aiarray_t, labeled_aiarray_init


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


cdef cbs_stat_t _determine_cbs_stat(double[::1] values):
    cdef int length = values.size
    cdef cbs_stat_t cbs_stat = calculate_cbs_stat(&values[0], length)

    return cbs_stat

def determine_cbs_stat(double[::1] values):
    cdef cbs_stat_t cbs_stat = _determine_cbs_stat(values)

    return cbs_stat


cdef double _determine_t_stat(double[::1] values, int i):
    cdef int length = values.size
    cdef double tstat = calculate_tstat(&values[0], i, length)

    return tstat

def determine_t_stat(double[::1] values, int i):
    cdef double tstat = _determine_t_stat(values, i)

    return tstat


cdef cbs_stat_t _determine_cbs(double[::1] values, int shuffles, double p):
    cdef int length = values.size
    cdef cbs_stat_t cbs_stat = calculate_cbs(&values[0], length, shuffles, p)

    return cbs_stat

def determine_cbs(double[::1] values, int shuffles=1000, double p=0.05):
    cdef cbs_stat_t cbs_stat = _determine_cbs(values, shuffles, p)

    return cbs_stat


cdef aiarray_t *_segment(double[::1] values, int shuffles, double p):
    cdef int length = values.size
    cdef aiarray_t *L = calculate_segment(&values[0], length, shuffles, p)

    return L

cdef void _segment_labeled(double[::1] values, labeled_aiarray_t *c_lsegments, char *label, int shuffles, double p):
    cdef int length = values.size
    calculate_segment_labeled(&values[0], length, c_lsegments, label, shuffles, p)

def segment(np.ndarray values, np.ndarray labels=None, int shuffles=1000, double p = 0.05):
    """
    Implementation of a Circular Biniary Segmentation
    algorithm to segment values into Intervals

    Parameters
    ----------
        values : numpy.ndarray
			Floats to segment
        labels : numpy.ndarray
            Labels for breakup segmentation
        shuffles : int
            Number of shuffles conduct
        p : double
            Pvalue cutoff

    Returns
    -------
        segments : IntervalArray or LabeledIntervalArray
			Segment intervals
    """

    # Call C segmentation function
    cdef char *label
    cdef aiarray_t *c_segments
    cdef labeled_aiarray_t *c_lsegments
    if labels is None:
        c_segments = _segment(values, shuffles, p)
    else:
        # Finf unique labels
        labels = labels.astype(bytes)
        unique_labels = pd.factorize(labels)
        c_lsegments = labeled_aiarray_init()
        for label in range(len(unique_labels)):
            label_values = values[labels==label]
            _segment_labeled(label_values, c_lsegments, label, shuffles, p)

    # Wrap C intervals
    cdef IntervalArray segments
    cdef LabeledIntervalArray lsegments
    if labels is None:
        segments = IntervalArray()
        segments.set_list(c_segments)

        return segments
    else:
        lsegments = LabeledIntervalArray()
        lsegments.set_list(c_lsegments)

        return lsegments


cdef aiarray_t *_validate(double[::1] values, IntervalArray segments, int shuffles, double p):
    cdef int length = values.size
    cdef aiarray_t *SV = calculate_validate(&values[0], length, segments.c_aiarray, shuffles, p)

    return SV

cdef void _validate_labeled(double[::1] values, LabeledIntervalArray segments, labeled_aiarray_t *vsegments, char *label, int shuffles, double p):
    cdef int length = values.size
    calculate_validate_labeled(&values[0], length, segments.ail, vsegments, label, shuffles, p)
    
def validate(double[::1] values, segments, int shuffles=1000, double p = 0.05):
    """
    Implementation of a Circular Biniary Segmentation validation
    algorithm to validate segment Intervals

    Parameters
    ----------
        values : numpy.ndarray
			Floats to segment
        labels : numpy.ndarray
            Labels for breakup segmentation
        segments : IntervalArray or LabeledIntervalArray
            Predicted segments to validate
        shuffles : int
            Number of shuffles conduct
        p : double
            Pvalue cutoff

    Returns
    -------
        valid_segments : IntervalArray or LabeledIntervalArray
			Valid segment intervals
    """

    # Call Validate C function
    cdef int i
    cdef char *label
    cdef aiarray_t *c_vsegments
    cdef labeled_aiarray_t *c_vlsegments
    cdef LabeledIntervalArray label_segments
    cdef double[::1] label_values
    if isinstance(segments, IntervalArray):
        c_vsegments = _validate(values, segments, shuffles, p)
    else:
        # Factorize labels
        i = 0
        for label in segments.unique_labels():
            label_segments = segments.get(label)
            label_values = values[i:len(label_segments)]
            c_vlsegments = labeled_aiarray_init()
            _validate_labeled(label_values, label_segments, c_vlsegments, label, shuffles, p)

    # Wrap C intervals
    cdef IntervalArray vsegments
    cdef LabeledIntervalArray vlsegments
    if isinstance(segments, IntervalArray):
        vsegments = IntervalArray()
        vsegments.set_list(c_vsegments)
        
        return vsegments
    else:
        vlsegments = LabeledIntervalArray()
        vlsegments.set_list(c_vlsegments)
        
        return vlsegments
