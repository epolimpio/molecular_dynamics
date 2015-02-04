#!/usr/bin/env python
import matplotlib as matl
import matplotlib.pyplot as plt
import numpy as np
import pylab as P

with open("momentum.dat") as f:
    data = map(float, f)
print np.mean(data)

n, bins, patches = P.hist(data, 50, normed=1, histtype='stepfilled')
P.setp(patches, 'facecolor', 'g', 'alpha', 0.75)

y = P.normpdf( bins, 0, 1)
l = P.plot(bins, y, 'k--', linewidth=1.5)

P.show()