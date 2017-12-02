#!/usr/bin/env python

'''
try poisson fitting and gaussian
'''
import numpy as np
import curvefit
from sys import argv
import pyfits
imA_name = argv[-1]
imA = pyfits.getdata(imA_name, dtype=float)
paras, cov =curvefit.poisson_fitting("poisson_fitting:{0}".format(imA_name),imA,VERBOSE=4)
#print(curvefit.poisson(2,2))
