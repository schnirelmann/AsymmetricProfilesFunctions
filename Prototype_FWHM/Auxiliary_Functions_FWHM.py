# -*- coding: utf-8 -*-

import numpy as np
import Functions_Asymmetric_Peaks as fap
from functools import partial
import matplotlib.pyplot as plt

#%%
""" 
FUNCTION 
"""
def bisection_search(f, x_low, x_high, y_search, precision):
     """
     Bisection search over monotone function. 
     
     Input: parameters such that f(x_low) <= y_search <= f(x_high). f must be 
     monotone on the interval [x_low, x_high].
     
     Output: the first output is x_0 such that |f(x_0) - y_search| < precision
     and the second output is |f(x_0) - y_search| itself.
     """

     assert(f(x_low) < y_search)
     assert(f(x_high) > y_search)
     
     i = 0
     x_l = x_low
     x_h = x_high
     top_iterations = 1000
     
     mid_point = (x_l + x_h)/2.0
     f_mid = f(mid_point) - y_search
     
     while i < top_iterations and np.abs(f_mid) > precision:
          mid_point = (x_l + x_h)/2.0
          f_mid = f(mid_point) - y_search
          if np.sign(f_mid) == -1:
               x_l = mid_point
          else:
               x_h = mid_point
          i += 1
     
     return mid_point, np.abs(f_mid)
#%%
## Test
#mu_0 = 0.0
#fwhm_0 = 1.0
#a_0 = 5
#eta_0 = 0.8
#scale_0 = 1000.0
#plot_window_0 = 10.0
#epsilon = 0.00000001
#
#fix_voigt = partial(fap.asymmetric_voigt, mu = 0.0, fwhm = fwhm_0, eta = eta_0, scale = scale_0, a = a_0)
#x_low, x_high, y_search = 4.0, 0.0, fix_voigt(0.0)/2.0
#
#bisection_search(fix_voigt, x_low, x_high, y_search, epsilon)

#%%
""" 
FUNCTION 
"""
def find_x_below_halb_max_in_voigt(fwhm_0, eta_0, scale_0, a_0):
     """ Returns x, z such that voigt(x) and voigt (z) < voigt(0.0)/2 (peaks_height/2). The voigt function takes parameters mu= 0.0 , fwhm, eta, scale, a """
     fix_voigt = partial(fap.asymmetric_voigt, mu = 0.0, fwhm = fwhm_0, eta = eta_0, scale = scale_0, a = a_0)
     halb_max = fix_voigt(0.0)/2.0
     
     # Tentative values
     step = fwhm_0/2.0
     left_low  = -step
     right_low =  step
     v_ll = fix_voigt(left_low)
     v_rl = fix_voigt(right_low)
     
     jump = 2.0
     while v_ll >  halb_max:
          left_low = - jump * step
          v_ll = fix_voigt(left_low)
          jump += 1.0
     
     jump = 2.0
     while v_rl > halb_max:
          right_low = jump * step
          v_rl = fix_voigt(right_low)
          jump += 1.0
     
     return left_low, right_low
#%%
# Test
#mu_0 = 0.0
#fwhm_0 = 1.0
#a_0 = 5
#eta_0 = 0.8
#scale_0 = 1000.0
#plot_window_0 = 10.0
#epsilon = 0.00000001
#
#fix_voigt = partial(fap.asymmetric_voigt, mu = 0.0, fwhm = fwhm_0, eta = eta_0, scale = scale_0, a = a_0)
#x_low_left, x_low_right = find_x_below_halb_max_in_voigt(fwhm_0, eta_0, scale_0, a_0)
#x_high_left, x_high_right = 0.0, 0.0 
#halve_height = fix_voigt(0.0)/2.0
#
#bisection_search(fix_voigt, x_low_left, x_high_left, halve_height, epsilon)

#%%
""" 
FUNCTION 
"""
def print_and_plot_results(x_fwhm_left, precision_left, x_fwhm_right, precision_right, plot_window, mu_0, fwhm_0, eta_0, scale_0, a_0, halve_height, total_fwhm):
     print("\n ------------------------ Summary of Results --------------------")
     print("\n fwhm_left, precision : ", np.abs(x_fwhm_left), " ---- ", precision_left)
     print("\n fwhm_right, precision: ", np.abs(x_fwhm_right)," ---- ", precision_right)
     print("\n fwhm_total           : ", total_fwhm, "\n")
     
     fixed_voigt = partial(fap.asymmetric_voigt, mu = mu_0, fwhm= fwhm_0, eta = eta_0, scale = scale_0, a =  a_0)
     
     angle_range = np.linspace(mu_0 - plot_window, mu_0 + plot_window, int(plot_window) * 100)
     model_vec = [ fixed_voigt(x) for x in angle_range]
     plt.plot(angle_range, model_vec)
     #plt.plot(angle_range, model_vec_sym)
     plt.axhline(y= halve_height, color='black', linestyle='-')
     plt.axvline(x= mu_0, color='black', linestyle='-')
     plt.plot([x_fwhm_left + mu_0, x_fwhm_right + mu_0],[fixed_voigt(x_fwhm_left + mu_0), fixed_voigt(x_fwhm_right + mu_0)], "o", color = "r")
     plt.figure()
     return

#%%
""" 
MAIN FUNCTION 
"""
def fwhm_left_right_algorithm(mu_0, fwhm_0, eta_0, scale_0, a_0, precision, plot_window):
     
     # print parameters
     information_string = "\n ------------------ Parameters -------------------" "\n mu = " + str(mu_0) + "\n fwhm = " + str(fwhm_0) + "\n a = " + str(a_0) + "\n eta = " + str(eta_0) + "\n scale = " + str(scale_0) + "\n plot_window = " + str(plot_window) + "\n precision = " + str(precision)

     print(information_string)
     
     # points to the left and right of the halve peak's height
     x_low_left, x_low_right = find_x_below_halb_max_in_voigt(fwhm_0, eta_0, scale_0, a_0)
     x_high_left, x_high_right = 0.0, 0.0
     
     # Voigt with given parameters and peak position 0.0
     fixed_voigt = partial(fap.asymmetric_voigt, mu = 0.0, fwhm= fwhm_0, eta = eta_0, scale = scale_0, a =  a_0)
     # Halve height of peak
     halve_height = fixed_voigt(0.0)/2.0
     
     # Left FWHM
     x_fwhm_left, precision_left = bisection_search(fixed_voigt, x_low_left, x_high_left, halve_height, precision)
     
     # Right FWHM
     x_fwhm_right, precision_right = bisection_search(fixed_voigt, x_low_right, x_high_right, halve_height, precision)
     
     # Total FWHM
     total_fwhm = np.abs(x_fwhm_left) + np.abs(x_fwhm_right)
     
     # Plot and print results
     print_and_plot_results(x_fwhm_left, precision_left, x_fwhm_right, precision_right, plot_window, mu_0, fwhm_0, eta_0, scale_0, a_0, halve_height, total_fwhm)

     return x_fwhm_left, x_fwhm_right, total_fwhm
#%%
## Test
#
#mu_0 = -20.0
#fwhm_0 = 10.0
#a_0 = 10.0
#eta_0 = 0.2
#scale_0 = 10000.0
#plot_window_0 = 30.0
#epsilon = 0.0000000001
#
#fwhm_left_right_algorithm(mu_0, fwhm_0, eta_0, scale_0, a_0, epsilon, plot_window_0)