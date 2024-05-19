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
    def ImaginarySpherical(self,l,m):
       #Creade grid of phi and theta angles for ploting surface mesh
        phi, theta = np.linspace(0, np.pi, 100), np.linspace(0, 2*np.pi, 100)
        phi, theta = np.meshgrid(phi, theta)

        #Calcualte spherical harmonic with given m and l
        Ylm = sph_harm(m, l, theta, phi)
        R=abs(Ylm)

        # Let's normalize color scale
        fcolors    = Ylm.imag
        fmax, fmin = fcolors.max(), fcolors.min()
        fcolors    = (fcolors - fmin)/(fmax - fmin)
        if m>l:
                raise ValueError("l and m must be postive, it means, that l can't be smaller than m")
        return R,phi,theta,fcolors,l,m 
    def normalized_wf(self,d,n,l,m):
        r     = np.linspace(0, d, 10000)
        pr = electron.radial(self,r, n, l)**2 * r**2 * (r[1]-r[0])
        max_r = r[ np.where(np.cumsum(pr) >0.95)[0][0] ]

        # Set coordinates grid to assign a certain probability to each point (x, y) in the plane
        x = y = np.linspace(-max_r, max_r, 501)
        x, y = np.meshgrid(x, y)
        r = np.sqrt((x ** 2 + y ** 2))

        # Ψnlm(r,θ,φ) = Rnl(r).Ylm(θ,φ)
        psi = electron.radial(self,r, n, l) * sph_harm(m, l, 0, np.arctan(x / (y + 1e-7)))
        psi_sq = np.abs(psi)**2
        return psi_sq,max_r,n,l,m
    def Cartesian_definition(self):
      xyz = np.linspace(-10, 10, 51)
      x,y,z = np.meshgrid(xyz, xyz, xyz, sparse=False)
      return x,y,z
    def CartesianToSpherical(self,x,y,z):
        r = np.sqrt(x**2 + y**2 + z**2)
        phi = np.arctan2(y+1e-10, x)
        theta = np.where( np.isclose(r, 0.0), np.zeros_like(r), np.arccos(z/r) )
        return r,phi,theta
    def normalized_wf3D(self,n,l,m):
        x,y,z = electron.Cartesian_definition(self)
        r,phi,theta = electron.CartesianToSpherical(self,x,y,z)
        psi = electron.radial(self,r, n, l) * np.real(sph_harm(m, l, phi, theta))
        return psi


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
    def plot_spherical_imaginary(self,package=None,R=None,theta=None,phi=None,fcolors=None,l=None,m=None):

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
    def plot_wf_2d(self,package=None,psi_sq=None,max_r=None,n=None,l=None,m=None):
        fig, ax = plt.subplots()
        if package:
            ax.contour(package[0], 40, cmap='RdBu', extent=[-package[1], package[1],-package[1], package[1]])
            ax.set_title(r'$|\psi_{{({0}, {1}, {2})}}|^2$'.format(package[2], package[3], package[4]), fontsize=15)
            ax.set_aspect('equal')
            ax.set_xlabel(r'$x$')
            ax.set_ylabel(r'$y$')
            ax.set_xlim(-package[1], package[1])
            ax.set_ylim(-package[1], package[1])
            plt.show()
        else:
            ax.contour(psi_sq, 40, cmap='RdBu', extent=[-max_r, max_r,-max_r, max_r])
            ax.set_title(r'$|\psi_{{({0}, {1}, {2})}}|^2$'.format(n, l, m), fontsize=15)
            ax.set_aspect('equal')
            ax.set_xlabel(r'$x$')
            ax.set_ylabel(r'$y$')
            ax.set_xlim(-max_r, max_r)
            ax.set_ylim(-max_r, max_r)
            plt.show()
    def plot_wf_3d(self,psi):
        x,y,z = electron.Cartesian_definition(electron)
        fig= go.Figure(data=go.Isosurface(
        x=x.flatten(),
        y=y.flatten(),
        z=z.flatten(),
        value=abs(psi).flatten(),
        colorscale='RdBu',
        isomin=-.75*abs(psi).min(),
        isomax=.75*abs(psi).max(),
        surface_count=6,
        opacity=0.5,
        caps=dict(x_show=False,y_show=False,z_show=False)
        ))
        fig.update_layout(paper_bgcolor="black",scene=dict(
                        xaxis=dict(visible=False),
                        yaxis=dict(visible=False),
                        zaxis=dict(visible=False)),
                        font=dict(
                        color="white"))
        fig.show()
