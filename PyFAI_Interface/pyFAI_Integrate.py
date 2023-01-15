 #!/usr/bin/env python3
# *-* coding: utf-8 *-*

#Authors: 
# Paulo Ricardo Garcia, paulo.garcia@lnls.br
# Joao Paulo Castro Zerba, joao.zerba@lnls.br
# Software to perform azimuthal integration of scattering images using pyFAI python library.
#Date: 22/11/2021

import h5py
import pyFAI, fabio, os, sys
import numpy as np
from pyFAI import azimuthalIntegrator
from os.path import expanduser
	
from PyQt5 import QtCore, QtWidgets, uic


class PyFAI_interface:
	def __init__(self, Int_pars, dat_folder, flip_msk):
		'''
		This module does the curves integration by using pyFAI library.
		'''
		# This check is to init the module from gui to be able to
		# read the poni file in order to write in the geometric fields 
		# in the graphical interface.
		self.flip_mask = flip_msk
		if type(Int_pars) == dict:
			# Loading all parameters from dictionary.
			self.loadingDictParams(Int_pars)
			# Set folder path for the dat files of the integration.
			self.dat_folder = dat_folder
			# Initializing the Intensity list of a curve.
			self.Iq = []
			# Initializing the q-values list of a curve.
			self.q = []
			# Initializing the standard deviation list of a curve.
			self.std = []
			# Initializing the incremented curve name.
			self.Output_names = []
			# Initializing the curve index.
			self.I_Index = 0

			self.Load_Geom()# Loading Geometric Parameters.
		
		else:
			pass

	def loadingDictParams(self, Int_pars):
		# This method loads all parameters sent when launching the integration.

		self.folder_path =Int_pars['folderpath']
		image_range = Int_pars['imgrange']
		self.initial_r = Int_pars['qrangeinitial']
		self.final_r = Int_pars[ 'qrangefinal']
		self.initial_azimuth = float(Int_pars['azimuthinitial'])
		self.final_azimuth = float(Int_pars['azimuthfinal'])
		self.poni_file = Int_pars['ponifile']
		self.mask_file = Int_pars['maskfile']
		self.radial_unity = Int_pars['qrange_unit']
		self.qbins =  int(Int_pars['qbins']) 
		self.hdf_group = Int_pars['hdfpath']
		self.Erro_Model = Int_pars['error_model']
		self.SolidAngle_correction = Int_pars['solidangle']
		self.Normalization_Factor = float(Int_pars['normfactor'])
		self.Polarization_factor = Int_pars['polarization']
		self.Spline = Int_pars['splinefile']
		if not self.poni_file:
			self.S_to_D_distance  = float(Int_pars['sampdetdist'])*1E+3
			self.X_c = float(Int_pars['xcenter'])
			self.Y_c = float(Int_pars['ycenter'])
			self.Xp_size =  float(Int_pars['xsize'])
			self.Yp_size =  float(Int_pars['ysize'])
			self.Tilt =  float(Int_pars['tilt'])
			self.Tilt_plan =   float(Int_pars['tiltplan'])
		
			self.Wvl = float(Int_pars['wave'])*1E-10
		
		
		if self.initial_r and self.final_r:
			self.QRange_bool  = True
			self.initial_r = float(self.initial_r)
			self.final_r = float(self.final_r)
			
		else:
			self.QRange_bool  = False
		if self.Polarization_factor:
			self.Polarization_factor = int(self.Polarization_factor)
			
		if self.Spline:
			self.Spline = int(self.Spline)
			
		# Converting images intervals to integer.
		# The intervals can be like x1-x2;x3-x4 ....
		self.images_interval =  [0] * len(image_range)
		for i in range(len(image_range)):
			if image_range[i] == 'All':
				self.images_interval[i] = [-1]
			else:
				self.images_interval[i] = self.Extract_intervals(image_range[i])
			filename = str(Int_pars["folderpath"]) + "/"+str(Int_pars["datafile"][i])
			self.Check_hdf_dimm(filename, self.hdf_group)
			Nimg = self.Hdf_Size( filename, self.hdf_group)
			
            
			for j in range(len(self.images_interval[i])):
				self.images_interval[i][j] = self.Interval_fine_tooth_comb(self.images_interval[i][j], Nimg)
				

	def Interval_fine_tooth_comb(self,Interval, nimg):
		#function to correct the image intervals if necessary

		Array_bool = isinstance(Interval,list)
		#check if the interval var is an array
		#the only case in which it is not an array is when
		#integrating all the images, in this case it is just = -1

		if Array_bool and len(Interval)>1:

			if Interval[0] > Interval[1]:
				temp = Interval[1]
				Interval[1] = Interval[0]
				Interval[0] = temp

			if Interval[1] > nimg -1:
				Interval[1] = nimg -1
			if Interval[0] < 0:
				Intevral[0] = 0

		elif Array_bool == True and  Interval[0] != -1 and len(Interval) == 1:
			if Interval[0] > nimg -1:
				Interval[0] = nimg -1
			elif Interval[0] < 0:
				Interval[0] = 0
			else:
				pass
		else:
			pass
		return Interval 
		
	def get_Images_number(self,img_number):
		#To import the number of images from MainWindow class
		#It just works in the single mode

		self.NImg = img_number 

	def find_string_character(self,s, ch):
		# This method helps to extract the image intervals.
		# It will find and return the semicolon indexes postions.
		return [i for i, letter in enumerate(s) if letter == ch]

	def Extract_intervals(self, img_range_str):

		#---------------------------------------------------------------------------------------------------#
		#Separating strings based on ';'
     
		semicolon_indexes = self.find_string_character(img_range_str, ";")
		
		if len(semicolon_indexes) == 0:
			img_range = [0]*1
			img_range[0] = img_range_str
		else:
			img_range = [0]*(len(semicolon_indexes) +1)
			img_range[0] = img_range_str[0: semicolon_indexes[0]]
			for i in range(1,len( semicolon_indexes)):
				img_range[i] = img_range_str[ semicolon_indexes[i-1]+1: semicolon_indexes[i]]
			img_range[len(img_range)-1] = img_range_str[ semicolon_indexes[-1]+1:len(img_range_str)]
    
   		#-----------------------------------------------------------------------------------------------------#

		# Breaking each string based on '-' and converting to float
		numerical_img_range = [0] * len(img_range)
		for i in range (len(img_range)):
			dash_index = self.find_string_character(img_range[i], "-")

			if dash_index :
				numerical_img_range[i] = [int(img_range[i][0:dash_index[0]]) ,int(img_range[i][dash_index[0] +1:len(img_range[i])])]
			else:
				numerical_img_range[i] = [int(img_range[i])]
		
     	#-----------------------------------------------------------------------------------------------------#    	
		return numerical_img_range

	def Plotto_OneGraph(self, q,I,sigma, XScale, YScale, q_unity, File_name):
		# This function generates the png files of
		# the scattering intensity graphs
		# Method not in use.

		import matplotlib.pyplot as plt	
		import matplotlib.ticker as mticker
		from matplotlib.ticker import FormatStrFormatter 
		from cycler import cycler		
	
		if q_unity == "q_A^-1":
			q_unity = "q 1/A"
		if q_unity == "q_nm^-1":
			q_unity = "q 1/nm"

		plt.figure()
		if XScale == 2:
			plt.xscale('log')
		if YScale == 2:
			plt.yscale('log')
	
		plt.errorbar(q,I,yerr = sigma,color = 'k',fmt='-o', markersize=3)
		plt.xlabel(q_unity)
		plt.ylabel('Intensity (a. u.)')
		plt.savefig(File_name+'.png', format = 'png')
		plt.close()

	def Group_checker(self, f, group):
		# Do this in the interface level?
		# This method checks the group path of the dataset in hdf5.
		# if the group is correctly given, it is ok otherwise it 
		# changes to the standard entry/data/data

		if group in f:
			# print("True group")
			pass
		elif "entry/data/data" in f:
			# print("True entry/data/data")
			group = "entry/data/data"
		elif "entry/data" in f:
			group = "entry/data"
		else:
			pass

		return group

	def Check_hdf_dimm(self, file_path, group):
		# This method is called in Qthread to check ndim and this ndim is then sent to callpyfai.
        # file_path is full file path .

		Name, self.Ext = os.path.splitext(file_path)
		
		try:
			if self.Ext == ".hdf5" or self.Ext == ".h5" or self.Ext == ".cxi":
				with h5py.File(file_path, "r", swmr = True) as f:

					group = self.Group_checker(f, group)
					ndim = len(f[group].shape)
					# print("ndim = ", ndim)

			elif self.Ext == ".npy":
				data_temp = np.load(file_path)
				if len(data_temp.shape) == 3:
					ndim = len(data_temp)
				elif len(data_temp.shape) == 2:
					ndim = 1
				else:
					pass
			elif self.Ext == ".tif":
					ndim = 1
				
			return ndim
		except Exception as e:
			print("Error trying to open file. Please check if it exists in the path folder.\nError: ", e)
			return -1
			pass
		
	def Hdf_Size(self, file_name, group):
		# This method is called from qthread to check sizes of hdf and npy files.

		try:
			if self.Ext == ".hdf5" or self.Ext == ".h5" or self.Ext == ".cxi":
				with h5py.File(file_name, "r", swmr = True) as f:
					group = self.Group_checker(f, group)#updates the group path
					return f[group].shape[0]

			elif self.Ext == ".npy":
				data_temp = np.load(file_name)
				if len(data_temp.shape) == 3:
					return len(data_temp)
				elif len(data_temp.shape) == 2:
					return 1
				else:
					pass
			elif self.Ext == ".tif":
				return  1
			else:
				pass

		except Exception as e:
			print("\t Could not read the file!",e)

	def Load_file(self, file_name, group, ndim, j):
		# It reads the hdf or npy file and store the images in numpy array.

		# print("Load file: ",file_name, group, ndim, j)
		try:
			if self.Ext == ".hdf5" or self.Ext == ".h5" or self.Ext == ".cxi":
				with h5py.File(file_name, "r", swmr = True) as f:
					
					# To check group path. There are two types in use in the lab "data" and "entry/data/data".
					# updates the group path.
					
					group = self.Group_checker(f, group)
					if ndim == 1 or ndim == 2:# when shape is 2.
						hdfData = np.array(f[group])
						# print("shape hdfdata = ", hdfData.shape)
						return hdfData
					elif ndim == 3:
						hdfData = np.array(f[group][j])
						# print("shape hdfdata = ", hdfData.shape)
					elif ndim == 4:
						hdfData = np.array(f[group][j][0])
						# print("shape hdfdata = ", hdfData.shape)
					else:
						pass

			elif self.Ext == ".npy":
				
				hdfData = np.load(file_name)
				if len(hdfData.shape) == 3:
					return hdfData[j]
				else:
					return hdfData
			elif self.Ext ==".tif":
				from PIL import Image
				tifData = np.array(Image.open(file_name))
				return tifData

			else:
				pass
			return hdfData#TODO check for errors

		except Exception as e:
			print("\t Could not read the hdf file! ",e)
			return None		
	
	def Load_mask(self, Mask_name):
	
		# Load the mask file.
		if self.flip_mask:	
			return np.flip((fabio.open(Mask_name)).data,0)
		else:
                         Name, Ext = os.path.splitext(Mask_name)
                         if Ext == '.msk':
                             return (fabio.open(Mask_name)).data
                         else:
                             return np.load(Mask_name)
                            
		
	def Load_geometry_IntegrationPars(self, Wvl, S_to_D_distance, X_c, Y_c, Tilt, Tilt_plan, Xp_size, Yp_size, Spline_filename):
		# Loading geometrical integration parameters with pyFAI library.

		try:
			IntP = azimuthalIntegrator.AzimuthalIntegrator(dist=0, poni1=0, poni2=0, rot1=0, rot2=0, rot3=0, pixel1=None, pixel2=None, splineFile=None, detector=None, wavelength=Wvl)	 
			IntP.setFit2D(S_to_D_distance, X_c, Y_c, tilt=Tilt, tiltPlanRotation=Tilt_plan, pixelX=Xp_size, pixelY=Yp_size, splineFile=Spline_filename)
		
		except Exception as e:
			print("\t Could not load geometry from IntegrationPars.py!", e)
		
		return IntP

	def Load_geometry_poni(self, poni_file):
		# This method loads the poni file generated using pyFAI calib.

		try:
			poni = pyFAI.load(poni_file)
 		 
		except Exception as e:
			print("\t Could not read the poni file!", e)
		
		return poni

	def Call_pyfai(self, img_idx, file_name, file_path, scatt, Mask_Arr, ndim):
		# This method is the one called for integrating the curves in the run method
		# in qthread of the gui script.

		#print("\nCall_pyfai file_name = ", file_name,"\n")	
		#print("\nCall_pyfai file_path = ", file_path,"\n")	
	
		It_Name =  file_name+'_'+ '{0:05d}'.format(img_idx) + '.dat'
		OutPut_name = self.dat_folder + '/' + It_Name 

		I_data = self.Load_file(file_path, self.hdf_group, ndim, img_idx)
		self.Iq.append(0)
		self.q.append(0)
		self.std.append(0)
		self.Output_names.append(It_Name)
		
		# Use the argument  filename = OutPut_name, if it is wanted to have the dat files saved.
		# If the q range is defined ...
		if self.QRange_bool:
			self.q[self.I_Index ],self.Iq[self.I_Index ],self.std[self.I_Index ] = scatt.integrate1d(I_data, self.qbins, mask = Mask_Arr, 
			unit = self.radial_unity,radial_range = (self.initial_r, self.final_r), azimuth_range = (self.initial_azimuth,self.final_azimuth), 
			error_model = self.Erro_Model, correctSolidAngle = self.SolidAngle_correction, 
			normalization_factor= self.Normalization_Factor )
		# Use the argument  filename = OutPut_name, if it is wanted to have the dat files saved.
		# ... otherwise it will be automatically set	
		else:
			self.q[self.I_Index ],self.Iq[self.I_Index ],self.std[self.I_Index ] = scatt.integrate1d(I_data, self.qbins, mask = Mask_Arr, 
			unit = self.radial_unity,azimuth_range = (self.initial_azimuth,self.final_azimuth), 
			error_model = self.Erro_Model, correctSolidAngle = self.SolidAngle_correction,
			normalization_factor= self.Normalization_Factor )

		#Finally, plotting the results
		#self.Plotto_OneGraph(self.q[self.I_Index ],self.Iq[self.I_Index ],self.std[self.I_Index ], 2, 2, self.radial_unity, It_Name )
		self.I_Index +=1
		#list with the first entry being the % value, and the second the file number index in datafile
		#print('i pyfai = ', i)
				
	def Load_Geom(self):
		# Loading mask file
		try:
			self.Mask_Arr = self.Load_mask(self.mask_file)
		except Exception as e:
			print("\t Could not read the mask file!", e)
	
		# Loading geometry from poni file ...
		if self.poni_file:	
			self.scatt = self.Load_geometry_poni(self.poni_file)
		# ... or using the geometric parameters defined in the graphical interface
		else:
			self.scatt = self.Load_geometry_IntegrationPars(self.Wvl, self.S_to_D_distance, self.X_c, self.Y_c, self.Tilt, self.Tilt_plan, self.Xp_size, self.Yp_size, self.Spline)

	#def datExporter(self):
		
