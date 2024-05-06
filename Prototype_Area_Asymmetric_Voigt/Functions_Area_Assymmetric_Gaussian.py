# -*- coding: utf-8 -*-

from scipy.integrate import quad
import Functions_Normalized_Asymmetric_Peaks as fnap


#%%
"""
FUNCTION
"""
def area_asymmetric_gaussian_left_right_total(a, W):
     """
     Returns a numeric aproximation using "quad" in scipy.integrate of the left, right and total area for the normalized asymmetric gaussian (the area remains the same also if it's not normalized). The suggested value for W is 20 when |a| < 100) 
     """
     
     area_left_tuple = quad(fnap.normalized_asymmetric_gaussian, -W, 0, args= (a))
     area_right_tuple = quad(fnap.normalized_asymmetric_gaussian, 0, W, args= (a))
     area_total = area_left_tuple[0] + area_right_tuple[0]
     precision_total = area_left_tuple[1] + area_right_tuple[1]
     area_total_tuple = (area_total, precision_total)
     
     return area_left_tuple, area_right_tuple, area_total_tuple

# Test
#area_asymmetric_gaussian_left_right_total(a = 1.0, W = 20)

