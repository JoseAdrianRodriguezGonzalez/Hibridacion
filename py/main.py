from class_electron import electron
from class_electron import plot
e=electron()
print(e.radial(1,2,0))
print(e.RealSpherical(1,-1))
p=plot()
p.plot_radial_(2,4,0)
p.plot_spherical_real(e.RealSpherical(1,0))