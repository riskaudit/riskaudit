# %% Run the following cell to initialize the API. The output will contain instructions on how to grant this notebook access to Earth Engine using your account.
import ee
ee.Authenticate()
ee.Initialize()
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
from datetime import datetime
from datetime import timedelta
import time
# from multiprocessing import cpu_count
# from multiprocessing.pool import ThreadPool
%matplotlib inline
# %% create function to download the image or imagecollection as you desire
def downloader(ee_object,region): 
    try:
        #download image
        if isinstance(ee_object, ee.image.Image):
            print('Its Image')
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
# %%
def download_url(path, newfile):
    t0 = time.time()
    url = path
    fn = newfile
    try:
        r = requests.get(url)
        with open(fn, 'wb') as f:
            f.write(r.content)
        return(url, time.time() - t0)
    except Exception as e:
        print('Exception in download_url():', e)
def ymdList(imgcol):
    def iter_func(image, newlist):
        date = ee.Number.parse(image.date().format("YYYYMMdd"));
        newlist = ee.List(newlist);
        return ee.List(newlist.add(date).sort())
    ymd = imgcol.iterate(iter_func, ee.List([]))
    return list(ee.List(ymd).reduce(ee.Reducer.frequencyHistogram()).getInfo().keys())
# def download_parallel(args):
#     cpus = cpu_count()
#     results = ThreadPool(cpus - 1).imap_unordered(download_url, args)
#     for result in results:
#         print('url:', result[0], 'time (s):', result[1])
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t0 = time.time()
from datetime import datetime
from time import mktime

meteor_path =  '/Users/joshuadimasaka/Desktop/PhD/GitHub/riskaudit/data/groundtruth/METEOR_PROJECT_2002'       
output_path = '/Users/joshuadimasaka/Desktop/PhD/GitHub/riskaudit/data/obsvariables/METEOR_PROJECT_2002/SENTINEL1-DUAL_POL_GRD_HIGH_RES'

country_list = os.listdir(meteor_path); country_list.sort()
if '.DS_Store' in country_list: country_list.remove('.DS_Store')
for icountry in country_list:
    icountry = country_list[0]
    geoJSON_path = meteor_path + '/' + icountry + '/tiles/extents'
    filenamelist = os.listdir(geoJSON_path); filenamelist.sort()
    if '.DS_Store' in filenamelist: filenamelist.remove('.DS_Store')
    for filename in filenamelist:
        # filename = filenamelist[0]
        result_path = output_path+'/'+icountry+'/'+filename[:-9]
        if not os.path.exists(result_path):
            os.makedirs(result_path)

        print(filename)
        with open(geoJSON_path+'/'+filename) as fa:
            geoJSON = json.load(fa)      
        coords = geoJSON['features'][0]['geometry']['coordinates']
        aoi = ee.Geometry.Polygon(coords)
        region = aoi.toGeoJSONString()
        startDATE = ee.Date('2014-06-14')
        endDATE = ee.Date('2023-12-31')
        im_coll = (ee.ImageCollection('COPERNICUS/S1_GRD_FLOAT')
                    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
                    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
                    .filter(ee.Filter.eq('instrumentMode', 'IW'))
                    .filter(ee.Filter.eq('resolution', 'H'))
                    .filterBounds(aoi)
                    .filterDate(startDATE,endDATE)
                    .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
                    .sort('system:time_start'))    
        orbitN = im_coll.aggregate_array('relativeOrbitNumber_start').getInfo() 
        im_coll = im_coll.filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
        im_list = im_coll.toList(im_coll.size())
        acq_times = im_coll.aggregate_array('system:time_start').getInfo()
        ymdlistvariable = ymdList(im_coll)
        ymd_year = [el[:4] for el in ymdlistvariable]
        indexes = [(x, ymd_year.index(x)) for x in set(ymd_year)]

        # urls = []
        # fns = []
        for i in range(len(indexes)):
            print(i)
            im1 = ee.Image(im_list.get(indexes[i][1])).select('VH').clip(aoi)
            im2 = ee.Image(im_list.get(indexes[i][1])).select('VV').clip(aoi)
            # urls.append(downloader(im1,region))
            # urls.append(downloader(im2,region))
            # fns.append(str(result_path+'/'+ymdlistvariable[indexes[i][1]]+'_'+str(orbitN[0])+"_VH.tif"))
            # fns.append(str(result_path+'/'+ymdlistvariable[indexes[i][1]]+'_'+str(orbitN[0])+"_VV.tif"))

            t0 = time.time()
            result = download_url(
                downloader(im1,region),
                str(result_path+'/'+ymdlistvariable[indexes[i][1]]+'_'+str(orbitN[0])+"_VH.tif")
                )
            print('url:', result[0], 'time:', result[1])

            t0 = time.time()
            result = download_url(
                downloader(im2,region),
                str(result_path+'/'+ymdlistvariable[indexes[i][1]]+'_'+str(orbitN[0])+"_VV.tif")
                )
            print('url:', result[0], 'time:', result[1])

        # download_parallel(zip(urls, fns))     
# %%
