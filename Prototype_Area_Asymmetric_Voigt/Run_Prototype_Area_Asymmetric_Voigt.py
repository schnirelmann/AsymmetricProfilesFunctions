# -*- coding: utf-8 -*-

import Functions_Area_Asymmetric_Voigt as faav

#%%

a_0 = 10.0
eta_0 = 0.5
scale_0 = 10000.0
integration_window_gaussian = 30

#%%
# Calculate and Print Areas
faav.testing_function_area_asymmetric_voigt(a_0, eta_0, scale_0, integration_window_gaussian)