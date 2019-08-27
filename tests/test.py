from cbseg import determine_cbs_stat, determine_t_stat, determine_cbs, segment, validate
import numpy as np
import pycbs
import time


print("Creating data")
np.random.seed(10)
x = np.random.random(300000)
x[10000:20000] = x[10000:20000] + 0.1
x[25000:27000] = x[25000:27000] - 1
print("mean:", np.mean(x), "\n\n")


print("CBS stat calulation")
start = time.perf_counter()
t, s, e = pycbs.cbs_stat(x)
print("py time:", time.perf_counter() - start)
print(s,e,t)

start = time.perf_counter()
stat = determine_cbs_stat(x)
print("c time:", time.perf_counter() - start)
print(stat['max_start'], stat['max_end'], stat['max_t'])
print("\n")


print("T stat calulation")
start = time.perf_counter()
t = pycbs.tstat(x, 100)
print("py time:", time.perf_counter() - start)
print(t)

start = time.perf_counter()
t = determine_t_stat(x, 100)
print("c time:", time.perf_counter() - start)
print(t)
print("\n")


print("CBS")
start = time.perf_counter()
thresh, t, s, e = pycbs.cbs(x)
print("py time:", time.perf_counter() - start)
print(s,e,t, thresh)

start = time.perf_counter()
stat = determine_cbs(x)
print("c time:", time.perf_counter() - start)
print(stat['max_start'], stat['max_end'], stat['max_t'])
print("\n")


print("Segment")
start = time.perf_counter()
L1 = pycbs.segment(x)
print("py time:", time.perf_counter() - start)
print(L1)

start = time.perf_counter()
L2 = segment(x)
print("c time:", time.perf_counter() - start)
L2.display()
print('\n')


print("Validate")
start = time.perf_counter()
SV1 = pycbs.validate(x, L1)
print("py time:", time.perf_counter() - start)
print(SV1)

start = time.perf_counter()
SV2 = validate(x, L2)
print("c time:", time.perf_counter() - start)
SV2.display()
print('\n')