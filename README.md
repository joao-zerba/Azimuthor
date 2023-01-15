# azimuthor
Azimuthal Integrator Software

<p>Azimuthor is a graphical user interface (GUI) for azimuthal average of Xray/neutron/electron 
scattering data using PyFAI python library. The idea of this GUI is to make the integration 
process more intuitive and, in this way, accessible and practical to the user.</p>
<p>The input are the geometry parameters, the mask and the scattering image. These
parameters can be obtained using FIT2D or pyFAI tools.</p>
<p>The mask file should be in .msk or npy formats. The scattering images file should be in
HDF5, H5, CXI, NPY and Tiff.</p>
<p>The software can work in single and batch modes, the last one is suitalble when dealing
with multiple files. It has an interactive display in which the integrated curves may be
vizualised and selected for exportation. The data are stored in hdf5 and dat files.</p>

Developers: 
<p>Paulo Ricardo Garcia <pauloricardoafg@yahoo.com.br></p>
<p>João Paulo Castro Zerba <joao.p.c.zerba@outlook.com></p>


**References:
<p>Fast Azimuthal Integration using Python. Available on https://pyfai.readthedocs.io/en/master/.
Accessed on May 12th,2022</p>
<p>Ashiotis, G., Deschildre, A., Nawaz, Z., Wright, J. P., Karkoulis, D., Picca, F. E., & Kieffer, J. (2015).
The fast azimuthal integration Python library: pyFAI. Journal of applied crystallography, 48(2), 510-519.</p>
<p>Hammersley, A. P. (2016). FIT2D: a multi-purpose data reduction, analysis and visualization
program. Journal of Applied Crystallography, 49(2), 646-652.</p>
<p>Kieffer, J., & Karkoulis, D. (2013, March). PyFAI, a versatile library for azimuthal regrouping. In Journal
of Physics: Conference Series (Vol. 425, No. 20, p. 202012). IOP Publishing.</p>
