# -*- coding: utf-8 -*-

import numpy as np
#import matplotlib.pyplot as plt

#%%
"""
FUNCTION
"""
def normalized_asymmetric_gaussian(x, a):
    x_star = x
    asym_fac = 1.0 + ( a * x_star )/( 1.0 + (1 + a**2) * x_star**2 )**0.5
    factor_one = 1.0 / (2.0 * np.pi)**0.5
    # Factor 2
    factor_two = np.e**( - 0.5 * ( x_star / asym_fac )**2)
    return factor_one * factor_two

# Test
#a_0 = 2.0
#W = 20.0
#plot_window = 30
#
#angle_range = np.linspace(- plot_window, plot_window, 1000)
#model_vec = [ normalized_asymmetric_gaussian(x, a_0) for x in angle_range]
#plt.plot(angle_range, model_vec)
#plt.figure()

#%%
"""
FUNCTION
"""
def normalized_asymmetric_cauchy(x, a):
    x_star = x
    asym_fac = 1.0 + ( a * x_star )/( 1.0 + (1 + a**2) * x_star**2 )**0.5
    numerator = 1.0 / np.pi
    denominator = 1.0 + (x_star / asym_fac)**2
    return numerator/denominator

# Test
#a_0 = 1.0
#plot_window = 20.0
#angle_range = np.linspace( -plot_window, plot_window, 1000)
#model_vec = [ normalized_asymmetric_cauchy(x, a_0) for x in angle_range]
#plt.plot(angle_range, model_vec)

#%%
"""
FUNCTION
"""
def normalized_asymmetric_voigt(x, eta, scale, a):
     # Gauss summand
     gaussian_summand = normalized_asymmetric_gaussian(x, a)
     
     # Cauchy summand
     cauchy_summand = normalized_asymmetric_cauchy(x, a)
     
     return scale * (eta * gaussian_summand + (1 - eta) * cauchy_summand)

#####Test
#a_0 = 1.0
#eta_0 = 0.7
#scale_0 = 10.0
#plot_window = 20.0
#angle_range = np.linspace( -plot_window, plot_window, 1000)
#model_vec = [ normalized_asymmetric_voigt(x, eta_0, scale_0, a_0) for x in angle_range]
#plt.plot(angle_range, model_vec)