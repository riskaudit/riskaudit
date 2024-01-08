# %% Run the following cell to initialize the API. The output will contain instructions on how to grant this notebook access to Earth Engine using your account.
# https://gorelick.medium.com/fast-er-downloads-a2abd512aa26
import ee
import multiprocessing
ee.Initialize(opt_url='https://earthengine-highvolume.googleapis.com')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm, gamma, f, chi2
import pandas as pd
import IPython.display as disp
import json
import csv 
import os
import datetime
import requests
import shutil
from retry import retry
from datetime import datetime
from datetime import timedelta
import time
from multiprocessing import cpu_count
from multiprocessing.pool import ThreadPool
%matplotlib inline
def ymdList(imgcol):
    def iter_func(image, newlist):
        date = ee.Number.parse(image.date().format("YYYYMMdd"));
        newlist = ee.List(newlist);
        return ee.List(newlist.add(date).sort())
    ymd = imgcol.iterate(iter_func, ee.List([]))
    return list(ee.List(ymd).reduce(ee.Reducer.frequencyHistogram()).getInfo().keys())
@retry(tries=10, delay=1, backoff=2)
def download_url(args):


    t0 = time.time()

    url = downloader(args[0],args[2])
    fn = args[1] 

    try:
        r = requests.get(url)
        with open(fn, 'wb') as f:
            f.write(r.content)
        return(url, time.time() - t0)
    except Exception as e:
        print('Exception in download_url():', e)
@retry(tries=10, delay=1, backoff=2)
def downloader(ee_object,region): 
    try:
        #download image
        if isinstance(ee_object, ee.image.Image):
            # print('Its Image')
            url = ee_object.getDownloadUrl({
                    'scale': 10,
                    'crs': 'EPSG:4326',
                    'region': region,
                    'format': 'GEO_TIFF'
                })
            return url
        
        #download imagecollection
        elif isinstance(ee_object, ee.imagecollection.ImageCollection):
            print('Its ImageCollection')
            ee_object_new = ee_object.mosaic()
            url = ee_object_new.getDownloadUrl({
                    'scale': 10,
                    'crs': 'EPSG:4326',
                    'region': region,
                    'format': 'GEO_TIFF'
                })
            return url
    except:
        print("Could not download")
@retry(tries=10, delay=1, backoff=2)
def download_parallel(args):
    cpus = cpu_count()
    results = ThreadPool(cpus - 1).imap_unordered(download_url, args)
    for result in results:
        print('url:', result[0], 'time (s):', result[1])
t0 = time.time()
from datetime import datetime
from time import mktime

meteor_path =  '/Users/joshuadimasaka/Desktop/PhD/GitHub/riskaudit/data/groundtruth/METEOR_PROJECT_2002'       
output_path = '/Users/joshuadimasaka/Desktop/PhD/GitHub/riskaudit/data/obsvariables/METEOR_PROJECT_2002/SENTINEL-2-MSI_LVL2A'

country_list = os.listdir(meteor_path); country_list.sort()
if '.DS_Store' in country_list: country_list.remove('.DS_Store')

# %%
for ic in range(len(country_list)): #range(2, 41): # len(country_list)):
    # icountry = country_list[ic]
    # icountry = country_list[custom_list[ic]]
    icountry = country_list[ic]
    geoJSON_path = meteor_path + '/' + icountry + '/tiles/extents'
    filenamelist = os.listdir(geoJSON_path); filenamelist.sort()
    if '.DS_Store' in filenamelist: filenamelist.remove('.DS_Store')
    for ifilename in range(len(filenamelist)): # range(custom start index,len(filenamelist)):
        # filename = filenamelist[0]
        filename = filenamelist[ifilename]
        result_path = output_path+'/'+icountry+'/'+filename[:-9]
        if not os.path.exists(result_path):
            os.makedirs(result_path)

        print(filename)
        with open(geoJSON_path+'/'+filename) as fa:
            geoJSON = json.load(fa)      
        coords = geoJSON['features'][0]['geometry']['coordinates']
        aoi = ee.Geometry.Polygon(coords)
        region = aoi.toGeoJSONString()

        # get direction and orbit number
        startDATE = ee.Date('2015-01-01')
        endDATE = ee.Date('2023-12-31')
        im_coll = (ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
                    .filterBounds(aoi)
                    .filterDate(startDATE,endDATE)
                    .sort('system:time_start'))      

        ymdlistvariable = ymdList(im_coll)
        ymd_year = [el[:4] for el in ymdlistvariable]
        uniq_year = list(map(int, list(set(ymd_year))))
        uniq_year.sort()

        ims = []
        fns = []
        rgns = []
        for i in range(len(uniq_year)):
            startDATE = ee.Date(str(uniq_year[i]) + '-01-01')
            endDATE = ee.Date(str(uniq_year[i]) + '-12-31')

            ims.append(im_coll.filterDate(startDATE,endDATE).select('B1').mean().clip(aoi))
            ims.append(im_coll.filterDate(startDATE,endDATE).select('B2').mean().clip(aoi))
            ims.append(im_coll.filterDate(startDATE,endDATE).select('B3').mean().clip(aoi))
            ims.append(im_coll.filterDate(startDATE,endDATE).select('B4').mean().clip(aoi))
            ims.append(im_coll.filterDate(startDATE,endDATE).select('B5').mean().clip(aoi))
            ims.append(im_coll.filterDate(startDATE,endDATE).select('B6').mean().clip(aoi))
            ims.append(im_coll.filterDate(startDATE,endDATE).select('B7').mean().clip(aoi))
            ims.append(im_coll.filterDate(startDATE,endDATE).select('B8').mean().clip(aoi))
            ims.append(im_coll.filterDate(startDATE,endDATE).select('B8A').mean().clip(aoi))
            ims.append(im_coll.filterDate(startDATE,endDATE).select('B11').mean().clip(aoi))
            ims.append(im_coll.filterDate(startDATE,endDATE).select('B12').mean().clip(aoi))

            fns.append(str(result_path+'/'+str(uniq_year[i])+"_B1_aerosol.tif"))
            fns.append(str(result_path+'/'+str(uniq_year[i])+"_B2_blue.tif"))
            fns.append(str(result_path+'/'+str(uniq_year[i])+"_B3_green.tif"))
            fns.append(str(result_path+'/'+str(uniq_year[i])+"_B4_red.tif"))
            fns.append(str(result_path+'/'+str(uniq_year[i])+"_B5_red1.tif"))
            fns.append(str(result_path+'/'+str(uniq_year[i])+"_B6_red2.tif"))
            fns.append(str(result_path+'/'+str(uniq_year[i])+"_B7_red3.tif"))
            fns.append(str(result_path+'/'+str(uniq_year[i])+"_B8_nir.tif"))
            fns.append(str(result_path+'/'+str(uniq_year[i])+"_B8A_red4.tif"))
            fns.append(str(result_path+'/'+str(uniq_year[i])+"_B11_swir1.tif"))
            fns.append(str(result_path+'/'+str(uniq_year[i])+"_B12_swir2.tif"))

            rgns.append(region)
            rgns.append(region)
            rgns.append(region)
            rgns.append(region)
            rgns.append(region)
            rgns.append(region)
            rgns.append(region)
            rgns.append(region)
            rgns.append(region)
            rgns.append(region)
            rgns.append(region)

        download_parallel(zip(ims, fns, rgns))
# %%
