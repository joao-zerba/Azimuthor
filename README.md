# azimuthor
Azimuthal Integrator Software

Azimuthor is a graphical user interface (GUI) for azimuthal average of Xray/neutron/electron 
scattering data using PyFAI python library. The idea of this GUI is to make the integration 
process more intuitive and, in this way, accessible and practical to the user. 
The input are the geometry parameters, the mask and the scattering image. These
parameters can be obtained using FIT2D or pyFAI tools
The mask file should be in .msk or npy formats. The scattering images file should be in
HDF5, H5, CXI, NPY and Tiff.
The software can work in single and batch modes, the last one is suitalble when dealing
with multiple files. It has an interactive display in which the integrated curves may be
vizualised and selected for exportation. The data are stored in hdf5 and dat files

Developers: 
Paulo Ricardo Garcia <pauloricardoafg@yahoo.com.br>
João Paulo Castro Zerba <joao.p.c.zerba@outlook.com>
