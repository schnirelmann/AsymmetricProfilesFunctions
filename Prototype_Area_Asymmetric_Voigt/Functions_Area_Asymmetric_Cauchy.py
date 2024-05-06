# -*- coding: utf-8 -*-

#import matplotlib.pyplot as plt
import numpy as np
import cmath

#%%
""" 
FUNCTION 
"""
def integral_function_x_3_5(a, root):
     numerator = 1/4 * a * (root**5 - root**3) *  np.log(1 - root)
     denominator = (3 * a**2 * root**4 - 3 * a**2 * root**2 + a**2 + root**6)
     return numerator/denominator

#%%
""" 
FUNCTION 
"""
def cauchy_integral_x_3_5(a):
     """
     Returns second part of the integral of asymmetric Cauchy.
     Attention: arithmetic with complex numbers must be enable!
     """
     # Factors to construct roots
     factor_0 = (2 * a**4 + 2 * a**2 + 2**(1/3) * (a**4 * (a**2 + 1)**2)**(1/3))**0.5

     factor_1 = 8 * a**4 + 8 * a**2 - 2 * 2**(1/3) * (a**4 * (a**2 + 1)**2)**(1/3)

     factor_2 = -64 * a**6 - 96 * a**4 - 32 * a**2
     
     # Root construction (they have imaginary parts)
     root_1 = -(-a**2 - factor_0/2**0.5 - 1/2 * (factor_1 - factor_2/(4 * 2**0.5 * factor_0))**0.5 )**0.5
     
     root_3 = -(-a**2 - factor_0/2**0.5 + 1/2 * (factor_1 - factor_2/(4 * 2**0.5 * factor_0))**0.5 )**0.5
     
     root_5 = -(-a**2 + factor_0/2**0.5 - 1/2 * (factor_1 + factor_2/(4 * 2**0.5 * factor_0))**0.5 )**0.5

     root_6 = (-a**2 + factor_0/2**0.5 - 1/2 * (factor_1 + factor_2/(4 * 2**0.5 * factor_0))**0.5 )**0.5
     
     # Result
     integral_result_x_3_5 = 2 * (integral_function_x_3_5(a, root_1).real + integral_function_x_3_5(a, root_3) +  integral_function_x_3_5(a, root_5).real + integral_function_x_3_5(a, root_6).real)
 
     return 1/np.pi * integral_result_x_3_5

#%%
""" 
FUNCTION 
"""
def integral_big_polynomial(a, root):
     # Numerator
     numerator = (3 * (a**2 + 1) * root**4 + (2 * a**4 + 3 * a**2 + 1) * root**6 + 3 * root**2 + 1) * np.log(-root)
     
     # Denominator
     denominator = (a**2 + 1)**2 * root**7 + 3 * (a**2 + 1)**2 * root**5 + 3 * (a**2 + 1) * root**3 + root
     
     # Result
     integral_result_big_polynomial = 1/8 * numerator/denominator 
     
     return integral_result_big_polynomial

#%%
""" 
FUNCTION 
"""
def cauchy_integral_big_polynomial(a):
     # Terms to calculate roots
     term_0 = (-4/(a**2 + 1) + (2.0 * 2**(1/3) * (np.abs(a))**(4/3))/(1 + a**2)**(4/3) + 4)**0.5
     
     term_1 = 96/(a**2 + 1.0) - 32/(a**2 + 1)**2 - 64.0
     
     term_2 = -8/(a**2 + 1) - (2.0 * 2**(1/3) * (np.abs(a))**(4/3))/ ( 1 + a**2 )**(4/3)
    
     # Calculation of roots
     root_1 = -complex(-0.5 * term_0 - 0.5 * complex(term_2 - (term_1)/(4 * term_0) + 8)**0.5 - 1.0)**0.5
     
     root_4 = complex(-0.5 * term_0 + 0.5 * complex(term_2 - (term_1)/(4 * term_0) + 8)**0.5 - 1)**0.5
     
     root_5 = -complex(0.5 * term_0 - 0.5 * complex(term_2 + (term_1)/(4 * term_0) + 8)**0.5 - 1)**0.5

     root_6 = complex(0.5 * term_0 - 0.5 * complex(term_2 + (term_1)/(4 * term_0) + 8)**0.5 - 1)**0.5
     
     # Calculate integral with roots     
     cauchy_integral_big_polynomial = 2 * ( integral_big_polynomial(a, root_1).real + integral_big_polynomial(a, root_4).real + integral_big_polynomial(a, root_5).real + integral_big_polynomial(a, root_6).real)

     return 1/np.pi * cauchy_integral_big_polynomial

#%%
""" 
FUNCTION 
"""
def asymmetric_cauchy_integral_right_left_total(a):
     """
     Attention "a" must be different of 0. 
     ToDo: add exception when a = 0.
     """
     # Auxiliary terms
     integral_big_polynomial = cauchy_integral_big_polynomial(a)
     integral_small_polynomial = cauchy_integral_x_3_5(a)
     
     # Asymmetric Cauchy Right 
     right_area_asymmetric_cauchy = np.abs(integral_big_polynomial + integral_small_polynomial)
     
     # Asymmetric Cauchy Left Area
     left_area_asymmetric_cauchy = np.abs(integral_big_polynomial - integral_small_polynomial)
     
     # Asymmetric Cauchy Total Area
     total_area_asymmetic_cauchy = 2 * np.abs(integral_big_polynomial)
     
     return left_area_asymmetric_cauchy, right_area_asymmetric_cauchy, total_area_asymmetic_cauchy

#%%
# Test
#test = asymmetric_cauchy_integral_right_left_total(1.0)
#
#print("\n test: ", a , test)
