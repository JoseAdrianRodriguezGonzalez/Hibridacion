import numpy as np
import  matplotlib.pyplot as plt
import plotly.graph_objects as go
from scipy.special import sph_harm, genlaguerre, factorial
class electron:
    def __init__(self):
        pass
    def radial(self,r,n,l):
        if n <= 0 or l < 0 or l >= n:
            raise ValueError("Values of n and l must be postive or being at the form of n > 0, l >= 0 y l < n",n,l)
        a0_radius = 1
        prefactor = np.sqrt( ((2 / n * a0_radius) ** 3 * (factorial(n - l - 1))) / (2 * n * (factorial(n + l))) )
        laguerre = genlaguerre(n - l - 1, 2 * l + 1)
        p = 2 * r / (n * a0_radius)
        return  prefactor * np.exp(-p / 2) * (p ** l) * laguerre(p)
    def RealSpherical(self,l,m):
        #Creade grid of phi and theta angles for ploting surface mesh
        phi, theta = np.linspace(0, np.pi, 100), np.linspace(0, 2*np.pi, 100)
        phi, theta = np.meshgrid(phi, theta)

        #Calcualte spherical harmonic with given m and l
        Y = sph_harm(abs(m), l, theta, phi)
        if m<0:
            Y=np.sqrt(2)*(-1)**m*Y.imag
        elif m>0:
            Y=np.sqrt(2)*(-1)**m*Y.real
        R=abs(Y)
        # Let's normalize color scale
        fcolors    = Y.real
        fmax, fmin = fcolors.max(), fcolors.min()
        fcolors    = (fcolors - fmin)/(fmax - fmin)
        if m>l:
            raise ValueError("l and m must be postive, it means, that l can't be smaller than m")
        return R,phi,theta,fcolors,l,m
class plot:
    def __init__(self, *args):
        super(plot, self).__init__(*args)
    def plot_radial_(self,d,n,l):
        r = np.linspace(0, d, 1000)
        R_function = electron.radial(electron,r, n, l)
        #plt.plot(r, Rnl)
        #plt.plot(r, Rnl**2)
        plt.plot(r, r**2 * R_function**2)
        plt.title(f'$R_{{n, l}}(r)$ distancia = {d}')
        plt.xlabel(r'$r [a_0]$')
        plt.ylabel(r'$R_nl(r) r^2$')
        plt.show()
    def plot_spherical_real(self,package=None,R=None,theta=None,phi=None,fcolors=None,l=None,m=None):
        if package:
            fig = go.Figure(data=[go.Surface(x=package[0]*np.sin(package[1]) * np.cos(package[2]),
                                     y=package[0]*np.sin(package[1]) * np.sin(package[2]),
                                     z=package[0]*np.cos(package[1]),
                                     surfacecolor=package[3],
                                     colorscale='balance')])

    # Show the plot
            fig.update_layout(title=fr'$Y_{package[4], package[5]}$', autosize=False,
                            width=700, height=700,
                            margin=dict(l=65, r=50, b=65, t=90),paper_bgcolor="black",
                            scene=dict(
                            xaxis=dict(visible=False),
                            yaxis=dict(visible=False),
                            zaxis=dict(visible=False)),
                            font=dict(
                            color="white"))
            fig.show()
        else:
            fig = go.Figure(data=[go.Surface(x=R*np.sin(phi) * np.cos(theta),
                                     y=R*np.sin(phi) * np.sin(theta),
                                     z=R*np.cos(phi),
                                     surfacecolor=fcolors,
                                     colorscale='balance')])

    # Show the plot
            fig.update_layout(title=fr'$Y_{l, m}$', autosize=False,
                            width=700, height=700,
                            margin=dict(l=65, r=50, b=65, t=90),paper_bgcolor="black",
                            scene=dict(
                            xaxis=dict(visible=False),
                            yaxis=dict(visible=False),
                            zaxis=dict(visible=False)),
                            font=dict(
                            color="white"))
            fig.show()
                