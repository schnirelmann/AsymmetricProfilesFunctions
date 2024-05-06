# -*- coding: utf-8 -*-

import Functions_Area_Assymmetric_Gaussian as faag
import Functions_Area_Asymmetric_Cauchy as faac

#%%
"""
FUNCTION
"""
def area_voigt(a, eta, scale, integration_window_gaussian):
     ######## Gaussian areas
     left_area_gaussian_t, right_area_gaussian_t, total_area_gaussian_t = faag.area_asymmetric_gaussian_left_right_total(a, integration_window_gaussian)
     left_area_gaussian = left_area_gaussian_t[0]
     right_area_gaussian = right_area_gaussian_t[0]
     total_area_gaussian = total_area_gaussian_t[0]
     total_precision_gaussian = total_area_gaussian_t[1]
     # Gaussian areas tuple
     g_areas_t = (left_area_gaussian, right_area_gaussian, total_area_gaussian, total_precision_gaussian)

     ######## Cauchy areas
     left_area_cauchy, right_area_cauchy, total_area_cauchy = faac.asymmetric_cauchy_integral_right_left_total(a)
     # Cauchy areas tuple
     c_areas_t = (left_area_cauchy, right_area_cauchy, total_area_cauchy)
     
     ######## Voigt areas
     left_area_voigt = scale * (eta * left_area_gaussian + (1 - eta) * left_area_cauchy)
     right_area_voigt = scale * (eta * right_area_gaussian + (1 - eta) * right_area_cauchy)
     total_area_voigt = left_area_voigt + right_area_voigt
     total_precision_voigt = (scale * eta) * total_precision_gaussian
     # Voigt areas tuple
     v_areas_t = (left_area_voigt, right_area_voigt, total_area_voigt, total_precision_voigt)
     
     return g_areas_t, c_areas_t, v_areas_t

#%%

def testing_function_area_asymmetric_voigt(a, eta, scale, integration_window):

     g_areas_t, c_areas_t, v_areas_t = area_voigt(a, eta, scale, integration_window)
     
     print("\n ------------------- Parameters -----------------------------")
     print(" a          = ", a) 
     print(" eta        = ", eta)
     print(" scale      = ", scale)
     print(" W_integral = ", integration_window)
          
     print("\n ------------------- Areas Normalized Gaussian --------------")
     print("\n Area_Gauss_left  : ", g_areas_t[0])
     print("\n Area_Gauss_right : ", g_areas_t[1])
     print("\n Area_Gauss_total : ", g_areas_t[2])
     print("\n Total_Precision  : ", g_areas_t[3])
     
     
     print("\n ------------------- Areas Normalized Cauchy ----------------")
     print("\n Area_Cauchy_left : ", c_areas_t[0])
     print("\n Area_Cauchy_right: ", c_areas_t[1])
     print("\n Area_Cauchy_total: ", c_areas_t[2])
     print("\n Total_Precision  : ", 0.0 )
     
     print("\n ------------------- Areas Voigt ----------------------------")
     print("\n Area_Voigt_left  : ", v_areas_t[0])
     print("\n Area_Voigt_right : ", v_areas_t[1])
     print("\n Area_Voigt_total : ", v_areas_t[2])
     print("\n Total_Precision  : ", v_areas_t[3])

     return

##### Test 
#a_0 = 10.0
#eta_0 = 0.5
#scale_0 = 10000.0
#integration_window_gaussian = 30
#
#testing_function_area_asymmetric_voigt(a_0, eta_0, scale_0, integration_window_gaussian)
