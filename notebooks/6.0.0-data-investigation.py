# %% import libraries
import rasterio
import os
import numpy as np
%matplotlib inline

# %% data directories
data_dir = r"C:\Users\admin\Desktop\GitHub\riskaudit\data\obsvariables\GOOGLE_EARTH\clipped"
files = [ f for f in os.listdir(data_dir) if f.endswith('.tif') ]

# %% open the raster file
fp = os.path.join(data_dir, files[0])
raster = rasterio.open(fp)

# %%
raster.crs

# %%
# Affine transform (how raster is scaled, rotated, skewed, and/or translated)
raster.transform

# %%
# Dimensions
print(raster.width)
print(raster.height)
raster.count
# %%
# Bounds of the file
raster.bounds

# %%
raster.driver

# %%
# No data values for all channels
raster.nodatavals

# %%
# All Metadata for the whole raster dataset
raster.meta

# %%
# Read the raster band as separate variable
band1 = raster.read(1)
print(type(band1))
print(band1.dtype)

# %%
from rasterio.plot import show
# %%
show((raster, 1))
# %%
import matplotlib.pyplot as plt
%matplotlib inline
# %%
fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, nrows=1, figsize=(10, 4), sharey=True)

# %%
show((raster, 1), cmap='Reds', ax=ax1)
show((raster, 2), cmap='Greens', ax=ax2)
show((raster, 3), cmap='Blues', ax=ax3)
# %%
# Add titles
ax1.set_title("Red")
ax2.set_title("Green")
ax3.set_title("Blue")
# %%
fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, nrows=1, figsize=(10, 4), sharey=True)
show((raster, 1), cmap='Reds', ax=ax1)
show((raster, 2), cmap='Greens', ax=ax2)
show((raster, 3), cmap='Blues', ax=ax3)
ax1.set_title("Red")
ax2.set_title("Green")
ax3.set_title("Blue")
# %%
# Read the grid values into numpy arrays
red = raster.read(1)
green = raster.read(2)
blue = raster.read(3)

# Function to normalize the grid values
def normalize(array):
    """Normalizes numpy arrays into scale 0.0 - 1.0"""
    array_min, array_max = array.min(), array.max()
    return ((array - array_min)/(array_max - array_min))

# Normalize the bands
redn = normalize(red)
greenn = normalize(green)
bluen = normalize(blue)

# Create RGB natural color composite
rgb = np.dstack((redn, greenn, bluen))

# Let's see how our color composite looks like
plt.imshow(rgb)

# %%
