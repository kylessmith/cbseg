# Circular Binary Segmentation

[![Build Status](https://travis-ci.org/kylessmith/cbseg.svg?branch=master)](https://travis-ci.org/kylessmith/cbseg) [![PyPI version](https://badge.fury.io/py/cbseg.svg)](https://badge.fury.io/py/cbseg)
[![Coffee](https://img.shields.io/badge/-buy_me_a%C2%A0coffee-gray?logo=buy-me-a-coffee&color=ff69b4)](https://www.buymeacoffee.com/kylessmith)

Circular Binary Segmentation (CBS) is an algorithm tests for change-points 
using a maximal t-statistic with a permutation reference distribution to obtain 
the corresponding p-value.


All citations should reference to [original paper][paper].

## Install

If you dont already have numpy and scipy installed, it is best to download
`Anaconda`, a python distribution that has them included.  
```
    https://continuum.io/downloads
```

Dependencies can be installed by:

```
    pip install -r requirements.txt
```

PyPI install, presuming you have all its requirements installed:
```
    NOT AVAILABLE...YET
	pip install cbseg
```

## Usage

```python
from cbseg import determine_cbs_stat, determine_t_stat, determine_cbs, segment, validate
import numpy as np

# Create data
np.random.seed(10)
x = np.random.random(300000)
x[10000:20000] = x[10000:20000] + 0.1
x[25000:27000] = x[25000:27000] - 1

# Calculate segments
L = segment(x, shuffles=1000, p=0.05)
for l in L:
   print(l)
# Interval(0-9998, 0, 0.8195453778139452)
# Interval(9998-20005, 0, 0.6922370217962109)
# Interval(20005-23672, 0, 0.37855730268017107)
# Interval(23672-24999, 0, 0.47507824637709734)
# Interval(24999-27000, 0, 1.106086021061645)
# Interval(27000-27765, 0, 0.20298421273826422)
# Interval(27765-300000, 0, 0.4615259909041863)

# Validate segments
SV = validate(x, L, shuffles=1000, p=0.05)
for l in SV:
   print(l)
# Interval(0-9998, 0, 59.217947014694886)
# Interval(9998-20005, 0, 24.51595009558689)
# Interval(20005-24999, 0, 1426.421130519047)
# Interval(24999-27000, 0, 214.73719913960934)
# Interval(27000-27765, 0, 1446.2331071442745)
# Interval(27765-300000, 0, 1446.2331071442745)

```


## Original paper

> E. S. Venkatraman, Adam B. Olshen; A faster circular binary segmentation algorithm for the analysis of array CGH data, Bioinformatics, https://doi.org/10.1093/bioinformatics/btl646


[paper]: https://academic.oup.com/bioinformatics/article/23/6/657/413480