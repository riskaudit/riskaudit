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
