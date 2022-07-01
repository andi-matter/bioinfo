# -*- coding: utf-8 -*-
"""
Created on Thu Jun 30 12:49:50 2022

@author: erns_ae
"""


import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.optimize import curve_fit

degrees = np.loadtxt("degrees_unbinned.txt")


density = True

# print(degrees)

bins = int(max(degrees)/7)

plt.hist(degrees, bins=bins, density=density, label="data")
# plt.xlim((0,2000))
# plt.show()

def powerlaw(x, gamma, a):
    global density 
    # an unnecessary distinction bc you need the scale factor anyway, just affects ylims in plots so what the heck ever
    if density:
        return  a * np.power(x, -gamma)
    else: return  a * np.power(x, -gamma)

def expon(x, gamma, a):
    return a * np.exp(x * -gamma)

degree_binned, bin_edges = np.histogram(degrees, bins=bins, density=density)

rev_deg_counts = degree_binned[::-1]
rev_deg = bin_edges[-2::-1]

popt, pcov = curve_fit(powerlaw, rev_deg, rev_deg_counts)

popt2, pcov2 = curve_fit(expon, rev_deg, rev_deg_counts)

gamma = popt[0]
a = popt[1]

gamma2 = popt2[0]
a2 = popt2[1]
print(popt)

xrange = np.arange(0, max(bin_edges), 0.5)
yrange = powerlaw(xrange, gamma, a)
yrange2 = expon(xrange, gamma2, a2)

plt.plot(xrange, yrange, label="powerlaw, gamma = {:.3f}, c = {:.3f}".format(gamma, a))
#plt.plot(xrange, yrange2, label="exponential")
plt.legend()

if density:
    plt.xlim((0,300))
    plt.ylim((0, 0.02))
plt.title("Degree PDF and histogram")
plt.savefig("histogram.pdf", dpi=72, bbox_inches="tight")