
# coding: utf-8

# # Joblib for Daniel:
# 
# Trying to implement parallelism into Daniels problem.

# 
# ## Some Tests with random values
# I don't know if these quite match your data types

# In[135]:

import time
import tempfile
import shutil
import os
import numpy as np

from joblib import Parallel, delayed
from joblib import load, dump


# In[136]:

def griddata(gridpoints, tlayer, teff_logg_feh, method='linear', rescale=True):
    """ Do what ever it does"""
    # put a short wait.
    time.sleep(0.5)
    return np.sum(tlayer) * teff_logg_feh[0] + teff_logg_feh[1] + teff_logg_feh[2]   # thing to test inputs


# In[137]:

def inside_loop(newatm, models, layer, column, gridpoints, teff_logg_feh):
    tlayer = np.zeros( len(models))
    for inx, model in enumerate(models):
        tlayer[indx] = model[layer, column]
    #print(" for layer = {0}, column = {1}".format(layer, column))
    print("[Worker {0:d}] Layer {1:d} and Column {2:d} is about to griddata".format(os.getpid(), layer, column))
    newatm[layer,column] = griddata(gridpoints, tlayer, teff_logg_feh, method='linear', rescale=True)
    


# In[138]:

layers = range(3)
columns = range(2)
gridpoints = 5
teff = 1000
logg = 1
feh = -0.01
model1 = np.array([[1,2],[3,4],[5,6]])
model2 = np.array([[7,8],[9,10],[11,12]])
models = [model1, model2, model1*2, model2*2] # random models


# In[139]:

#%%timeit 
newatm = np.zeros([len(layers), len(columns)])
generator = (inside_loop(newatm, models, layer, column, gridpoints, (teff, logg, feh)) for  layer in layers for column in columns)

for i in generator:
    #print(newatm)
    pass
print(newatm)


# In[ ]:

#%%timeit
# Turning parallel
newatm = np.zeros([len(layers), len(columns)])
print("newatm before parallel", newatm)
Parallel(n_jobs=-1, verbose=1) (delayed(inside_loop)(newatm, models, layer, column, gridpoints, (teff, logg, feh)) for layer in layers for column in columns)

time.sleep(0.5)
print("newatm after parallel", newatm)

# This runs in parallel but it does not return any data yet.

#Need to memmap the results


# ## Parallel over both loops with memapping
# Look here to implement the memmap to your solution:

# In[128]:


def inside_loop(newatm, models, layer, column, gridpoints, teff_logg_feh):
    tlayer = np.zeros( len(models))
    for inx, model in enumerate(models):
        tlayer[indx] = model[layer, column]
    newatm[layer,column] = griddata(gridpoints, tlayer, teff_logg_feh, method='linear', rescale=True)
    
def griddata(gridpoints, tlayer, teff_logg_feh, method='linear', rescale=True):
    """ Do what ever it does"""
    time.sleep(0.5)
    return True   # thing to test inputs

folder = tempfile.mkdtemp()
newatm_name = os.path.join(folder, 'newatm')
try:
    # Pre-allocate a writeable shared memory map as a container for the
    # results of the parallel computation
    newatm = np.memmap(newatm_name, dtype=model.dtype, shape=model.shape, mode='w+')  # need to adjsut the shape

    print("newatm before parallel", newatm)
    Parallel(n_jobs=-1, verbose=1) (delayed(inside_loop)(newatm, models, layer, column, gridpoints, (teff, logg, feh)) for layer in layers for column in columns)

    time.sleep(0.5)
    print("newatm after parallel", newatm)

finally:
    # deleting temp files after testing the reuslt in example 
        try:
            shutil.rmtree(folder)
        except:
            print("Failed to delete: " + folder)
            


# In[ ]:




# # Direct copy of Joblib memmaping example

# In[ ]:

def sum_row(input, output, i):
    """Compute the sum of a row in input and store it in output"""
    sum_ = input[i, :].sum()
    print("[Worker {0:d}] Sum for row {1:d} is {2:f}".format(os.getpid(), i, sum_))
    output[i] = sum_

if __name__ == "__main__":
    rng = np.random.RandomState(42)
    folder = tempfile.mkdtemp()
    samples_name = os.path.join(folder, 'samples')
    sums_name = os.path.join(folder, 'sums')
    try:
        # Generate some data and an allocate an output buffer
        samples = rng.normal(size=(10, int(1e6)))

        # Pre-allocate a writeable shared memory map as a container for the
        # results of the parallel computation
        sums = np.memmap(sums_name, dtype=samples.dtype,
                         shape=samples.shape[0], mode='w+')
        print("samples shape", samples.shape)
        # Dump the input data to disk to free the memory
        dump(samples, samples_name)

        # Release the reference on the original in memory array and replace it
        # by a reference to the memmap array so that the garbage collector can
        # release the memory before forking. gc.collect() is internally called
        # in Parallel just before forking.
        samples = load(samples_name, mmap_mode='r')

        # Fork the worker processes to perform computation concurrently
        Parallel(n_jobs=4)(delayed(sum_row)(samples, sums, i)
                           for i in range(samples.shape[0]))

        # Compare the results from the output buffer with the ground truth
        print("Expected sums computed in the parent process:")
        expected_result = samples.sum(axis=1)
        print(expected_result)

        print("Actual sums computed by the worker processes:")
        print(sums)

        assert np.allclose(expected_result, sums)
    finally:
        try:
            shutil.rmtree(folder)
        except:
            print("Failed to delete: " + folder)


# In[ ]:



