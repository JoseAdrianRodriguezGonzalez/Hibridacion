import numpy as np
import  matplotlib.pyplot as plt
import plotly.graph_objects as go
from scipy.special import sph_harm, genlaguerre, factorial
class electron:
    def __init__(self):
        pass
    def radial(self,r,n=1,l=0):
        if n <= 0 or l < 0 or l >= n:
            raise ValueError("Values of n and l must be postive or being at the form of n > 0, l >= 0 y l < n")
        a0_radius = 1
        prefactor = np.sqrt( ((2 / n * a0_radius) ** 3 * (factorial(n - l - 1))) / (2 * n * (factorial(n + l))) )
        laguerre = genlaguerre(n - l - 1, 2 * l + 1)
        p = 2 * r / (n * a0_radius)
        return  prefactor * np.exp(-p / 2) * (p ** l) * laguerre(p)