#!/user/bin/python
"""
    Demo of the histogram (hist) function with a few features.

    In addition to the basic histogram, this demo shows a few optional features:

        * Setting the number of data bins
        * The ``normed`` flag, which normalizes bin heights so that the integral of
          the histogram is 1. The resulting histogram is a probability density.
        * Setting the face color of the bars
        * Setting the opacity (alpha value).

    """
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma 
from scipy import optimize

# example data
mu = 100 # mean of distribution
sigma = 5 # standard deviation of distribution

def gaussian(x, mu, sig):
    return np.power(2 * np.pi , -0.5)*np.exp(-np.power(x - mu , 2.) / (2 * np.power(sig, 2.)))/sig

def poisson(x, mu):
    return np.power(mu , x)/gamma(x+1.)*np.exp(-mu)

num_bins = 50
# the histogram of the data
x = mu + sigma * np.random.randn(10000) 
n, bins, patches = plt.hist(x, num_bins, normed=1, facecolor='blue', alpha=0.5)
bin_middles=0.5*(bins[1:]+bins[:-1])
moments=(mu)
paras, cov=optimize.curve_fit(poisson,bin_middles,n,p0=moments)
# add a 'best fit' line
y = poisson(bins, paras[0])
plt.plot(bins, y, 'r--')
plt.xlabel('Smarts')
plt.ylabel('Probability')
plt.title(r'Histogram of IQ: $\mu={0}$, $\sigma={1}$'.format(paras[0],paras[0]**0.5))

# Tweak spacing to prevent clipping of ylabel
plt.subplots_adjust(left=0.15)
plt.show()
