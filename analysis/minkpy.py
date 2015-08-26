#importing dependencies
from __future__ import division     #enable division in python 2.7
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import ctypes
import c2raytools as c2t
import pickle
import multiprocessing as mp
import os

# defining the minkowski c-wrapper
def mink(input_array, box_size, nr_bins, low, high, output_array):
    input_array = input_array.astype('float32')
    lib_file = '../minkowski_files/MINKOWSKI2_PY/Minkowski_python/minkowski_python.so'
    lib = ctypes.cdll.LoadLibrary(lib_file)
    func = lib.minkowski
    func(ctypes.c_void_p(input_array.ctypes.data), ctypes.c_void_p(box_size.ctypes.data),\
         ctypes.c_int(nr_bins), ctypes.c_float(low) , ctypes.c_float(high) ,\
         ctypes.c_void_p(output_array.ctypes.data))

    return output_array


############################################################
# defining functions to use in the minkowski_python function
############################################################

#defining smoothing function
def smoothing_function(input_array, sigma):
    print 'performing smoothing with sigma = ' + str(sigma)
    for layer in xrange(0,len(input_array)-1):
        input_array[:,:,layer] = c2t.smooth_gauss(input_array[:,:,layer], sigma) #gaussian 2D smoothing

    tophat = c2t.tophat_kernel_3d((1,1,sigma))
    input_array = c2t.smooth_with_kernel(input_array, kernel = tophat) #tophat smoothing in line of sight
  
    return input_array

#defining function that removes edge of box
def remove_edge(input_array, x=0 , y=0, z=0):
        
    input_array = input_array[x:len(input_array)-x,
                y:len(input_array)-y,
                z:len(input_array)-z]       #removes edge x,y,z pixels from edge of box
            
    return input_array

#defining z correction function
def z_correction(input_array, z):
    
    print 'removing z-dependence'
    input_array = input_array*(10/(1+float(z)))**(1/2)

    return input_array

#defining complementary field function
def complement(input_array):  
    
    print 'transforming to complementary field'
    input_array = input_array*(-1)

    return input_array

#defines function that creates gaussian noise cube
def gaussian_noise(size, mean, std, smoothing, x, y, z):
        
        print 'create noise array with std = ' + str(std) + ' and mean = ' + str(mean) 
        random_field = np.random.normal(mean, std*np.sqrt(smoothing**3),(size,size,size))

        random_field_smoothed = smoothing_function(random_field,smoothing)
        sigma_0sm =  np.std(random_field_smoothed[x:len(random_field_smoothed)-x,
                                                  y:len(random_field_smoothed)-y,
                                                  z:len(random_field_smoothed)-z])

        new_random_field = (std/sigma_0sm)*random_field

        return new_random_field

#defines function that adds the noise cube
def add_gaussian_noise(input_array, noise_array):

    return (input_array+noise_array)

#defines function that saves images
def save_image(input_array, z, name, filename):

    plt.imshow(input_array[256,:,:], interpolation='none')
    plt.colorbar()
    plt.savefig(filename + '_images/' + str(z) + '_x_' + str(name) +'.png')
    plt.clf()

    plt.imshow(input_array[:,256,:], interpolation='none')
    plt.colorbar()
    plt.savefig(filename + '_images/' + str(z) + '_y_' + str(name) +'.png')
    plt.clf()

    plt.imshow(input_array[:,:,256], interpolation='none')
    plt.colorbar()
    plt.savefig(filename + '_images/' + str(z) + '_z_' + str(name) +'.png')
    plt.clf()


#####################################################################
# some preperatory steps before running the minkowski_python function
#####################################################################

#load redshifts from file
with open ('../redshift_files/xfrac_redshifts.txt', "r") as myfile:
    redshifts=myfile.read().splitlines()

#set upper value and lower value of redshifts
z_upper_limit = 11.313
z_lower_limit = 5

# selects redshift indexes between the limits
indexes = sp.where( (np.array([float(z) for z in redshifts]) <= z_upper_limit) & \
                   (np.array([float(z) for z in redshifts]) >= z_lower_limit))

#use the above indexes to select the desired redhifts
redshifts = redshifts[indexes[0][0]:indexes[0][-1] +1]

# load aditional data (global ion fraction etc.)
xfracs = np.loadtxt('../minkowski_files/global_ion_fraction.dat'
                  ,usecols = (0,1,2,3,4))

#store information in variables
xfracs = xfracs[:,:][0:len(redshifts)]
ion_frac_vol = xfracs[:,1][::-1]      # average ionization fraction of all cells per volume
ion_frac_mass = xfracs[:,2][::-1]     # average ionization fraction of all cells per mass
neutral_frac_vol = xfracs[:,3][::-1]  # average nutral fraction of all cells per volume
neutral_frac_mass = xfracs[:,4][::-1] # average nutral fraction of all cells per mass

#create directory to store images
data_filename = 'xibox_smth4_edge10'
os.mkdir(data_filename + '_images')

#create dictionary object to store output data
data_dict = {}

#create noise array. Comment/uncomment depending on choice
noise_array = gaussian_noise(size=600, mean=0, std=10, smoothing=10, x=10, y=10, z=10)
#noise_array = None

##############################################
# defining the multiprocess minkowski function
##############################################

def minkowski_python(z,i, noise_array, data_filename):
    
    # create local dictionary object
    local_data_dict = {}

    # the data we want to analyse, a 3D array of floats. Comment/uncomment depending on format
    #input_array = c2t.read_cbin('/media/gorgel/My Passport 3/simon_master/dT_' + str(z) + '.cbin')                   #cbin format
                                                                                             
    input_array = c2t.XfracFile('/media/gorgel/My Passport 3/simon_master/xfrac3d_' + str(z) + '.bin')               # bin format
    input_array = input_array.xi
	
    #select options by comment/uncomment
    #save_image(input_array, z, 'raw', data_filename)

    print 'running minkowski program for' + ' z = ' +str(z)
    
    #input_array = add_gaussian_noise(input_array, noise_array)      #add noise
    #input_array = z_correction(input_array, z)                              #remove z.dependence
    #save_image(input_array, z, 'noise', data_filename)
    #input_array = smoothing_function(input_array, sigma=4)                   #smooth data
    #save_image(input_array, z, 'smth', data_filename)
    #input_array = remove_edge(input_array, 10, 10, 10)
    #save_image(input_array, z, 'edge', data_filename)
    #input_array = subtract_mean(input_array)
    #input_array = complement(input_array)                                   #calculate on complementary field
    #save_image(input_array, z, 'compl', data_filename)

    #minkowski settings
    size_array = len(input_array)
    box_dimensions = np.array([size_array, size_array , size_array])        # the size of the data we want to analyse, a 1D array of three integers (same as the -x -y -z options)
    nr_bins = int(65)                               # integer: number of bins to use (same as -b option)
    low = float(0)                                  # float: lowest value of threshold (same as -l option)
    high = float(1)                                # float: highest values of threshold (same as -h option)
    output_array = np.zeros((2,nr_bins+1,5), dtype=np.float32)      #output array, (koenderink+crofton, number of bins +1, 5 = th,V0,V1,V2,V3)

    #start minkowski program
    print "starting minkowski calculation"
    mink(input_array, box_dimensions, nr_bins, high, low, output_array)

    # collect and order output data in dictionary object
    local_data_dict.update({z : {'thresholds': output_array[0][:,4], 'CV0': output_array[0][:,0], 'CV1': output_array[0][:,1], \
                         'CV2': output_array[0][:,2] , 'CV3': output_array[0][:,3], \
                      'KV0': output_array[1][:,0], 'KV1': output_array[1][:,1], \
                         'KV2': output_array[1][:,2] , 'KV3': output_array[1][:,3], \
                      'ion_frac_vol' : ion_frac_vol[i], 'ion_frac_mass' :ion_frac_mass[i],\
                      'neutral_frac_vol' : neutral_frac_vol[i], 'neutral_frac_mass' : neutral_frac_mass[i]}})

    return local_data_dict



#define number of cores for multiprocess. processes is the number of cores you want to distribute the work over
pool = mp.Pool(processes=2,maxtasksperchild=1)

# starts the minkowski miltiprocess function and saves result
results = [pool.apply_async(minkowski_python, args=(z, i, noise_array, data_filename)) for z,i in zip(redshifts, range(0,len(redshifts)))]
outputs = [p.get() for p in results]

#put output in data dictionary
for output in outputs:
	data_dict.update(output)

#save data to file (pickeld binary that can be opened with the pickle module)
with open(data_filename + '.pkl', 'wb') as outfile:
    pickle.dump(data_dict, outfile, pickle.HIGHEST_PROTOCOL)
