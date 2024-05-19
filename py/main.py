from class_electron import electron
from class_electron import plot
e=electron()
p=plot()
p.plot_radial_(2,4,0)
p.plot_spherical_real(e.RealSpherical(1,0))
p.plot_spherical_imaginary(e.ImaginarySpherical(2,-1))
p.plot_wf_2d(e.normalized_wf(500,2,1,-1))
p.plot_wf_3d(e.normalized_wf3D(3,2,0))