import numpy as np
import logging

log = logging.getLogger()
logging.basicConfig(level=logging.WARN)


def cbs_stat(x):
    """
    Given x, Compute the subinterval x[i0:i1] with the maximal segmentation statistic t. 
    Returns t, i0, i1
    """
    
    x0 = x - np.mean(x)
    n = len(x0)
    y = np.cumsum(x0)
    e0, e1 = np.argmin(y), np.argmax(y)
    i0, i1 = min(e0, e1), max(e0, e1)
    s0, s1 = y[i0], y[i1]

    return (s1 - s0)**2 * n / (i1 - i0 + 1) / (n + 1 - i1 + i0), i0, i1+1


def tstat(x, i):
    """
    Return the segmentation statistic t testing if i is a (one-sided)  breakpoint in x
    """

    n = len(x)
    s0 = np.mean(x[:i])
    s1 = np.mean(x[i:])

    return (n-i)*i/n*(s0-s1)**2


def cbs(x, shuffles=1000, p=0.05):
    """
    Given x, find the interval x[i0:i1] with maximal segmentation statistic t. Test that statistic against
    given (shuffles) number of random permutations with significance p.  Return True/False, t, i0, i1; True if
    interval is significant, false otherwise.
    """

    max_t, max_start, max_end = cbs_stat(x)

    if max_end-max_start == len(x):
        return False, max_t, max_start, max_end

    if max_start < 5:
        max_start = 0

    if len(x)-max_end < 5:
        max_end = len(x)

    thresh_count = 0
    alpha = shuffles*p
    xt = x.copy()
    for i in range(shuffles):
        np.random.shuffle(xt)
        threshold, s0, e0 = cbs_stat(xt)

        if threshold >= max_t:
            thresh_count += 1

        if thresh_count > alpha:
            return False, max_t, max_start, max_end

    return True, max_t, max_start, max_end


def rsegment(x, start, end, L=[], shuffles=1000, p=.05):
    """
    Recursively segment the interval x[start:end] returning a list L of pairs (i,j) where each (i,j) is a significant segment.
    """

    threshold, t, s, e = cbs(x[start:end], shuffles=shuffles, p=p)
    log.info('Proposed partition of {} to {} from {} to {} with t value {} is {}'.format(start, end, start+s, start+e, t, threshold))
    
    if (not threshold) | (e-s < 5) | (e-s == end-start):
        L.append((start, end))
    else:
        if s > 0:
            rsegment(x, start, start+s, L)
        
        if e-s > 0:
            rsegment(x, start+s, start+e, L)
        
        if start+e < end:
            rsegment(x, start+e, end, L)
    
    return L


def segment(x, shuffles=1000, p=.05):
    """
    Segment the array x, using significance test based on shuffles rearrangements and significance level p
    """

    start = 0
    end = len(x)
    L = []
    rsegment(x, start, end, L, shuffles=shuffles, p=p)

    return L


def validate(x, L, shuffles=1000, p=0.01):
    S = [x[0] for x in L]+[len(x)]
    SV = [0]
    left = 0
    for test, s in enumerate(S[1:-1]):
        t = tstat(x[S[left]:S[test+2]], S[test+1]-S[left])
        #print("   test:",test,", start:",S[left],", end:",S[test+2],", t:",t)
        log.info('Testing validity of {} in interval from {} to {} yields statistic {}'.format(S[test+1], S[left], S[test+2], t))
        threshold = 0
        thresh_count = 0
        site = S[test+1]-S[left]
        xt = x[S[left]:S[test+2]].copy()
        flag = True
        for k in range(shuffles):
            np.random.shuffle(xt)
            threshold = tstat(xt, site)
            if threshold > t:
                thresh_count += 1
            if thresh_count >= p*shuffles:
                flag = False
                log.info('Breakpoint {} rejected'.format(S[test+1]))
                break
        if flag:
            log.info('Breakpoint {} accepted'.format(S[test+1]))
            SV.append(S[test+1])
            left += 1
    SV.append(S[-1])

    return SV