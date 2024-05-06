# -*- coding: utf-8 -*-

#import pandas as pd
import numpy as np
#import scipy as sp
import matplotlib.pyplot as plt
#from functools import partial
#from timeit import default_timer as timer
from scipy.integrate import quad

#%%
# =============================================================================
# Symmetric models
# =============================================================================

def asymmetric_gaussian(x, mu, fwhm, a):
    sigma = fwhm / ( 2.0 * ( 2.0 * np.log(2) )**0.5)
#    sigma = fwhm
    x_star = (x - mu)/sigma
    asym_fac = 1.0 + ( a * x_star )/( 1.0 + (1 + a**2) * x_star**2 )**0.5
    factor_one = 1.0 / (sigma * (2.0 * np.pi)**0.5 )
    # Factor 2
    factor_two = np.e**( - 0.5 * ( x_star / asym_fac )**2)
    return factor_one * factor_two

## Test
#mu = 0.0
#fwhm = 1.0
#a = 5.0
#W = 5.0
#parameters = [ mu, fwhm, a]
#angle_range = np.linspace(mu - W, mu + W, 1000)
#model_vec = [ asymmetric_gaussian(x, *parameters) for x in angle_range]
#plt.plot(angle_range, model_vec)
##
###model_vec = [ asymmetric_gaussian(x, mu, fwhm, 0.0) for x in angle_range]
###plt.plot(angle_range, model_vec)
###plt.figure()
###
###model_vec = [ asymmetric_gaussian(x, mu, fwhm, 0.0) - asymmetric_gaussian(x, *parameters) for x in angle_range]
###plt.plot(angle_range, model_vec)
###plt.figure()
###
##I = quad(asymmetric_gaussian, -2*W, 2*W, args= (mu, fwhm, a))
##print("\n Area: ", I, "\n")
#
#area_list = []
#a_range = np.linspace(0, 100, 1000) 
#for a in a_range:
#     I = quad(asymmetric_gaussian, mu - 2*W, mu +  2*W, args= (mu, fwhm, a))
#     area_list.append(I[0])
#
#plt.plot(a_range, area_list)
#plt.figure()

#%% 
def asymmetric_cauchy(x, mu, fwhm, a):
    gamma = fwhm/2.0
#    gamma = fwhm
    x_star = (x - mu)/ gamma
    asym_fac = 1.0 + ( a * x_star )/( 1.0 + (1 + a**2) * x_star**2 )**0.5
    numerator = 1.0 / ( gamma * np.pi)
    denominator = 1.0 + (x_star / asym_fac)**2
    return numerator/denominator

# Test
#mu = 0.0
#fwhm = 1.0
#a = 1.0
#W = 10.0
#parameters = [ mu, fwhm, a]
#angle_range = np.linspace(mu - W, mu + W, 1000)
#model_vec = [ asymmetric_cauchy(x, *parameters) for x in angle_range]
#plt.plot(angle_range, model_vec)
#
##asymmetric_cauchy(-1, *parameters)
##asymmetric_cauchy(+1, *parameters)
##asymmetric_cauchy(+1, *parameters) + asymmetric_cauchy(-1, *parameters)
##asymmetric_cauchy( 1, mu, fwhm, 0.0)
##asymmetric_cauchy(-1, mu, fwhm, 0.0)
##asymmetric_cauchy(+1, mu, fwhm, 0.0) + asymmetric_cauchy(-1, mu, fwhm, 0.0)
#
#model_vec = [ asymmetric_cauchy(x, mu, fwhm, a) for x in angle_range]
#plt.plot(angle_range, model_vec)
#
#I_R = quad(asymmetric_cauchy, 0, 30000, args= (mu, fwhm, a))
#print("\n Area Right: ", I_R, "\n")
#
#I_L = quad(asymmetric_cauchy, -30000, 0, args= (mu, fwhm, a))
#print("\n Area Right: ", I_L, "\n")
#
#print("\n Total Area: ", I_R[0] + I_L[0], "\n")


#%%
def asymmetric_voigt(x, mu, fwhm, eta, scale, a):
     gaussian_summand = asymmetric_gaussian(x, mu, fwhm, a)
#     sigma = fwhm * (2.0 * ( 2.0 * np.log(2) )**0.5 )
#     gaussian_summand = asymmetric_gaussian(x, mu, sigma, a)
     cauchy_summand = asymmetric_cauchy(x, mu, fwhm, a)
#     gamma = fwhm * 2.0
#     cauchy_summand = asymmetric_cauchy(x, mu, gamma, a)
     return scale * (eta * gaussian_summand + (1 - eta) * cauchy_summand)
## Test
mu = 0.0
fwhm = 1.0
eta = 0.5
scale = 1.0
a = 2.0
W = 3*fwhm

#             mu, fwhm, eta, scale, a
#parameters = [ mu, fwhm, eta, scale, a]
#parameters_sym = [ mu, fwhm, eta, scale, 0.0]
#angle_range = np.linspace(mu - 2*W, mu + 2*W, 1000)
#model_vec = [ asymmetric_voigt(x, *parameters) for x in angle_range]
#model_vec_sym = [ asymmetric_voigt(x, *parameters_sym) for x in angle_range]
#plt.plot(angle_range, model_vec)
#plt.plot(angle_range, model_vec_sym)
#plt.figure()
##
#I = quad(asymmetric_voigt, -2*W, 2*W, args= (mu, fwhm, eta, scale, a))
#print("\n Area: ", I, "\n")
#
#area_list = []
#a_range = np.linspace(0, 100, 1000) 
#for a in a_range:
#     I = quad(asymmetric_voigt, mu - 2*W, mu +  2*W, args= (mu, fwhm, eta, scale, a))
#     area_list.append(I[0])
#
#plt.plot(a_range, area_list)
#plt.figure()
#%%
def unpack_arguments(argumets_list):
     len(argumets_list)
     arguments_unpacked = []
     for i in range(len(argumets_list)): 
          arguments_unpacked = arguments_unpacked + [*argumets_list[i]]
     return arguments_unpacked
#%%
def asymmetric_voigt_area_transfer_mix(x, *args):
    # args = mus, fwhms, etas, scales
    number_of_peaks = int(len(args)/5)
    #print("number of voigts: ", number_of_voigts)
    add_voigts = 0.0
    # Add Voigts
    for i in range(0, number_of_peaks):
        add_voigts = add_voigts + asymmetric_voigt(x, args[i*5],
                                        args[1 + i*5], 
                                        args[2 + i*5], 
                                        args[3 + i*5],
                                        args[4 + i*5])
    return add_voigts    


## TEST
## mus, fwhms, etas, scales, a_s
#n = 7
#mus = [ float(i) for i in range(n)]
#num_peaks = int(len(mus))
#fwhms = num_peaks*[0.01]
#etas = num_peaks*[0.3]
#scales = num_peaks*[10.0]
#a_s = num_peaks * [10.0]
#
#arguments = list(zip(mus, fwhms, etas, scales, a_s))
#
##result = list( map(voigt, [ [1.0, 0.0, 1.0, 0.5, 1.0], [1.0, 5.0, 1.0, 0.5, 1.0] ]) )
#
#arguments_unpacked = unpack_arguments(arguments)
#     
#angle_range = np.linspace(mus[0] - 1, mus[0] + n, 10000)
#
##begin = timer()
#gamma_vec = [asymmetric_voigt_area_transfer_mix( x, *arguments_unpacked ) for x in angle_range]
##end = timer() # stop timer
##print("\n time for Voigt_Mix function: ", end - begin, "\n")
#
#plt.plot(angle_range, gamma_vec)