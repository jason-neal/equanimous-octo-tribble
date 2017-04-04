
# coding: utf-8

# # Testing numpy Stride
# For snr calculation windowing

# In[21]:

import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
from __future__ import division
from numpy.lib import stride_tricks
import pandas as pd
import seaborn as sns
get_ipython().magic('matplotlib inline')


# In[22]:

fname = "Test_spectra.fits"
data = fits.getdata(fname)
hdr = fits.getheader(fname)

wl = data["Wavelength"]
I = data["Extracted_DRACS"]
# print(type(I))
print(I.dtype)
wl = np.array(wl, dtype="float64")  # Turn >f4 into float64
I = np.array(I, dtype="float64")  # Turn >f4 into float64
print(I.dtype)
print(I)


# In[ ]:

binsize = 100

# Try using stride on np.array

# striding
nums = np.arange(len(I), dtype="int")

print("itemsize", nums.itemsize, "dtype", nums.dtype)
hop_length = 1
# stride_tests with numbers
frame_length = binsize
num_frames = 1 + (len(nums) - frame_length) / hop_length
row_stride = nums.itemsize * hop_length  # *hopesize
print(frame_length)
print(num_frames)
print(row_stride)

col_stride = nums.itemsize

nums_strided = stride_tricks.as_strided(nums, shape=(num_frames, frame_length), strides=(row_stride, col_stride))

print("nums", nums)
print("nums_strided =", nums_strided)

# row wise transform
row_sum = np.sum(nums_strided, axis=1)
# print(row_sum)
snr = 1 / np.std(nums_strided, axis=1)
print(snr)


# In[ ]:

# with I
frame_length = binsize
num_frames = 1 + (len(I) - frame_length) / hop_length
row_stride = I.itemsize * hop_length  # *hopesize
print(frame_length)
print(num_frames)
print(row_stride)

col_stride = I.itemsize
I_strided = stride_tricks.as_strided(I, shape=(num_frames, frame_length), strides=(row_stride, col_stride))

# print("nums", I)
# print("nums_strided =", I_strided)

snr = 1 / np.std(I_strided, axis=1)
print(snr)


# In[ ]:

plt.plot(snr)
plt.show()


# In[23]:


def strided_snr(data, frame_length, hop_length=1):
    num_frames = 1 + (len(data) - frame_length)/hop_length
    row_stride = data.itemsize * hop_length  # *hopesize
    col_stride = data.itemsize
    data_strided = stride_tricks.as_strided(data, shape=(num_frames, frame_length), strides=(row_stride, col_stride))

    print("length of data_strided", len(data_strided))
    snr = 1/np.std( data_strided, axis=1)
    # print("frame_length", frame_length)
    # print("num_frames", num_frames)
    # print("len(snr)", len(snr))
    # print(snr)

    # zeropad to make uniform length of spectra
    missing_size = len(data) - len(snr)
    print("missing size", missing_size)
    before = missing_size // 2
    end = missing_size // 2
    if missing_size % 2 is not 0:
        print("missing size is not even")
    padded_snr = np.pad(snr, (before, end), "constant")
    # print("padded length", len(padded_snr))
    # print(padded_snr)
    return padded_snr


def strided_sum(data, frame_length, hop_length=1):
    num_frames = 1 + (len(data) - frame_length) / hop_length
    row_stride = data.itemsize * hop_length  # *hopesize
    col_stride = data.itemsize
    data_strided = stride_tricks.as_strided(data, shape=(num_frames, frame_length), strides=(row_stride, col_stride))

    print("length of data_strided", len(data_strided))
    print("binsize", frame_length)
    print("hop_length", hop_length)
    print(data_strided)
    total = np.sum(data_strided, axis=1)
    # print("frame_length", frame_length)
    # print("num_frames", num_frames)
    # print("len(snr)", len(snr))
    # print(snr)

    # zeropad to make uniform length of spectra
    missing_size = len(data) - len(total)
    pad_size = (len(data) - len(total)) // 2
    # print("missing size", missing_size)
    before = missing_size // 2
    end = missing_size // 2
    if missing_size % 2 is not 0:
        print("missing size is not even")
    padded_total = np.pad(total, (pad_size, pad_size), "constant")
    # print("padded length", len(padded_snr))
    # print(padded_snr)
    return padded_total


# This doesn't seem to work that well with pandas not sure why
# store_array = np.empty((1024, len(bins)), dtype=data.dtype)
# for i, bin in enumerate(bins):
#   store_array[:, i] = strided_snr(I, bin)



# In[30]:

# loop over the different bin sizes
bins = np.arange(3, 51, 2)
hopper = 1
store_list = []
for i, b in enumerate(bins):
    store_list.append(strided_snr(I, b, hop_length=hopper))
print("done")


# In[31]:

# print(store_array)
print(store_list)


# In[32]:

# turn into a pandas dataframe
# dataframe = pd.DataFrame(data=store_array, columns=range(1024), index=bins)
# dataframe = pd.DataFrame(store_array, index=bins, columns=list(range(1024)))
# print(dataframe)

# print(dataframe.dtypes)


# In[33]:

df_list = pd.DataFrame(store_list, index=bins, columns=np.round(wl, 2))
print(df_list)


# In[36]:

sns.set()
cmap = sns.diverging_palette(220, 10, as_cmap=True)
ax = sns.heatmap(store_list, cmap=cmap, xticklabels=200, vmax=300, vmin=10)
# ax = sns.heatmap(df_list)
# plt.xticks(np.arange(int(np.min(wl)), int(np.max(wl) + 1), 1.0))
ax.set(ylabel="Binsize", xlabel="Wavelenght")


# In[37]:

# seaborn heatmap plot
sns.set()
cmap = sns.diverging_palette(220, 10, as_cmap=True)

ax = sns.heatmap(df_list, xticklabels=200, vmax=300, vmin=10)
# ax = sns.heatmap(df_list)
# plt.xticks(np.arange(int(np.min(wl)), int(np.max(wl) + 1), 1.0))
ax.set(ylabel="Binsize",
       xlabel="Wavelenght")


# In[35]:

# ax = sns.heatmap(store_list)
wl[50]-wl[0]


# In[ ]:




# # test on known data

# In[17]:

data = np.arange(20)

binsizes = range(1, 6, 2)
store = []
# opt = np.get_printoptions()
# np.set_printoptions(threshold='nan')



for b in binsizes:
    store.append(strided_sum(data, b))

# np.set_printoptions(**opt)


# In[18]:

SNRrand = pd.DataFrame(store, index=binsizes)
print(SNRrand)


# In[19]:

sns.set()
# cmap = sns.diverging_palette(220, 10, as_cmap=True)

ax = sns.heatmap(SNRrand, xticklabels=20)
# ax = sns.heatmap(df_list)
# plt.xticks(np.arange(int(np.min(wl)), int(np.max(wl) + 1), 1.0))
ax.set(ylabel="Binsize",
       xlabel="Wavelenght")


# In[ ]:




# In[ ]:
