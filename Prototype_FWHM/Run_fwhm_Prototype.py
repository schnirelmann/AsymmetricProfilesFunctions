# -*- coding: utf-8 -*-

import Auxiliary_Functions_FWHM as affwhm

#%%

# Set parameters
mu_0 = -20.0
fwhm_0 = 10.0
a_0 = 10.0
eta_0 = 0.2
scale_0 = 10000.0
plot_window_0 = 30.0
epsilon = 0.0000000001 # required precision for fwhm

# Run FWHM Prototype
affwhm.fwhm_left_right_algorithm(mu_0, fwhm_0, eta_0, scale_0, a_0, epsilon, plot_window_0)