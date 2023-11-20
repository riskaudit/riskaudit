# %% install packages
import rasterio
from rasterio.mask import mask
import geopandas as gpd
from osgeo import gdal
from osgeo import ogr
from osgeo import gdalconst

# %%
os.getcwd()

# %%
epsgOut = 4326
shp_fn =  r"C:\Users\admin\Desktop\GitHub\riskaudit\data\groundtruth\PHL_QC_EMI_AUGUST2022\barangay_boundary\NAME_V2_North Fair.gpkg"
dataShp = gpd.read_file(shp_fn)
dataShp = dataShp.to_crs(epsg=epsgOut)
shp_temp = ogr.Open(dataShp.to_json())
mb_l = shp_temp.GetLayer()

