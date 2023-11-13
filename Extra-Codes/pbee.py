# This module contains functions to use in Performance-Based Earthquake Engineering (PBEE) analysis.
# Author: Alexis Contreras R.

# Importing modules
import numpy as np
from scipy.stats import norm
from scipy.stats import lognorm

def lognorm_hazard_curve(IM,lambda_IM):

