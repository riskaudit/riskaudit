# %% Run the following cell to initialize the API. The output will contain instructions on how to grant this notebook access to Earth Engine using your account.
# https://gorelick.medium.com/fast-er-downloads-a2abd512aa26
import ee
import multiprocessing
ee.Authenticate()
ee.Initialize()
ee.Initialize(opt_url='https://earthengine-highvolume.googleapis.com')
# %%
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
# %% create function to download the image or imagecollection as you desire
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
                    'scale': 30,
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
                    'scale': 30,
                    'crs': 'EPSG:4326',
                    'region': region,
                    'format': 'GEO_TIFF'
                })
            return url
    except:
        print("Could not download")
def download_parallel(args):
    cpus = cpu_count()
    results = ThreadPool(cpus - 1).imap_unordered(download_url, args)
    for result in results:
        print('url:', result[0], 'time (s):', result[1])
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t0 = time.time()
from datetime import datetime
from time import mktime

meteor_path =  '/Users/joshuadimasaka/Desktop/PhD/GitHub/riskaudit/data/groundtruth/METEOR_PROJECT_2002'       
output_path = '/Users/joshuadimasaka/Desktop/PhD/GitHub/riskaudit/data/obsvariables/METEOR_PROJECT_2002/HLSL30_LANDSAT89'

country_list = os.listdir(meteor_path); country_list.sort()
if '.DS_Store' in country_list: country_list.remove('.DS_Store')

# %%
for icountry in country_list:
    icountry = country_list[0]

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
        startDATE = ee.Date('2013-01-11')
        endDATE = ee.Date('2023-12-31')
        im_coll = (ee.ImageCollection('NASA/HLS/HLSL30/v002')
                    .filter(ee.Filter.lt('CLOUD_COVERAGE', 0.001))
                    .filter(ee.Filter.lt('SPATIAL_COVERAGE', 100))
                    .filterBounds(aoi)
                    .filterDate(startDATE,endDATE)
                    .sort('system:time_start'))    
        im_list = im_coll.toList(im_coll.size())
        acq_times = im_coll.aggregate_array('system:time_start').getInfo()
        ymdlistvariable = ymdList(im_coll)
        ymd_year = [el[:4] for el in ymdlistvariable]
        indexes = [(x, ymd_year.index(x)) for x in set(ymd_year)]

        ims = []
        fns = []
        rgns = []
        for i in range(len(indexes)):
            print(i)

            ims.append(ee.Image(im_list.get(indexes[i][1])).select('B1').clip(aoi))
            ims.append(ee.Image(im_list.get(indexes[i][1])).select('B2').clip(aoi))
            ims.append(ee.Image(im_list.get(indexes[i][1])).select('B3').clip(aoi))
            ims.append(ee.Image(im_list.get(indexes[i][1])).select('B4').clip(aoi))
            ims.append(ee.Image(im_list.get(indexes[i][1])).select('B5').clip(aoi))
            ims.append(ee.Image(im_list.get(indexes[i][1])).select('B6').clip(aoi))
            ims.append(ee.Image(im_list.get(indexes[i][1])).select('B7').clip(aoi))
            ims.append(ee.Image(im_list.get(indexes[i][1])).select('B9').clip(aoi))
            ims.append(ee.Image(im_list.get(indexes[i][1])).select('B10').clip(aoi))
            ims.append(ee.Image(im_list.get(indexes[i][1])).select('B11').clip(aoi))
            ims.append(ee.Image(im_list.get(indexes[i][1])).select('Fmask').clip(aoi))


            fns.append(str(result_path+'/'+ymdlistvariable[indexes[i][1]]+"_B1_coastalaerosol.tif"))
            fns.append(str(result_path+'/'+ymdlistvariable[indexes[i][1]]+"_B2_blue.tif"))
            fns.append(str(result_path+'/'+ymdlistvariable[indexes[i][1]]+"_B3_green.tif"))
            fns.append(str(result_path+'/'+ymdlistvariable[indexes[i][1]]+"_B4_red.tif"))
            fns.append(str(result_path+'/'+ymdlistvariable[indexes[i][1]]+"_B5_nir.tif"))
            fns.append(str(result_path+'/'+ymdlistvariable[indexes[i][1]]+"_B6_swir1.tif"))
            fns.append(str(result_path+'/'+ymdlistvariable[indexes[i][1]]+"_B7_swir2.tif"))
            fns.append(str(result_path+'/'+ymdlistvariable[indexes[i][1]]+"_B9_cirrus.tif"))
            fns.append(str(result_path+'/'+ymdlistvariable[indexes[i][1]]+"_B10_tirs1.tif"))
            fns.append(str(result_path+'/'+ymdlistvariable[indexes[i][1]]+"_B11_tirs2.tif"))
            fns.append(str(result_path+'/'+ymdlistvariable[indexes[i][1]]+"_Fmask.tif"))


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