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
output_path = '/Users/joshuadimasaka/Desktop/PhD/GitHub/riskaudit/data/obsvariables/METEOR_PROJECT_2002/SENTINEL1-DUAL_POL_GRD_HIGH_RES'

country_list = os.listdir(meteor_path); country_list.sort()
if '.DS_Store' in country_list: country_list.remove('.DS_Store')

# %%
# custom_list = [13,14,15,18]
custom_list = [37,38,39,40,41,44,45]
for ic in range(len(country_list)): #range(2, 41): # len(country_list)):
    # icountry = country_list[ic]
    # icountry = country_list[custom_list[ic]]
    ims = []
    fns = []
    rgns = []
    icountry = country_list[custom_list[ic]]
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
        im_coll1 = (ee.ImageCollection('COPERNICUS/S1_GRD_FLOAT')
                    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
                    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
                    .filter(ee.Filter.eq('instrumentMode', 'IW'))
                    .filter(ee.Filter.eq('resolution', 'H'))
                    .filterBounds(aoi)
                    .filterDate(startDATE,endDATE)
                    .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING'))
                    .sort('system:time_start'))      
        im_coll2 = (ee.ImageCollection('COPERNICUS/S1_GRD_FLOAT')
                    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
                    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
                    .filter(ee.Filter.eq('instrumentMode', 'IW'))
                    .filter(ee.Filter.eq('resolution', 'H'))
                    .filterBounds(aoi)
                    .filterDate(startDATE,endDATE)
                    .filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'))
                    .sort('system:time_start'))  
        uniq_year1 = list(map(int, list(set([el[:4] for el in ymdList(im_coll1)]))))
        uniq_year2 = list(map(int, list(set([el[:4] for el in ymdList(im_coll2)]))))
        year_from_1 = []
        year_from_2 = []
        if not im_coll1.aggregate_array('relativeOrbitNumber_start').getInfo():
            im_coll = im_coll2
        elif len(uniq_year1) >= len(uniq_year2):
            im_coll = im_coll1
            if len(list(set(uniq_year2).difference(uniq_year1))) > 0:
                year_from_2 = list(set(uniq_year2).difference(uniq_year1))
        elif len(uniq_year1) < len(uniq_year2):
            im_coll = im_coll2
            if len(list(set(uniq_year1).difference(uniq_year2))) > 0:
                year_from_1 = list(set(uniq_year1).difference(uniq_year2))

        orbitN = im_coll.aggregate_array('relativeOrbitNumber_start').getInfo() 
        im_coll = im_coll.filter(ee.Filter.eq('relativeOrbitNumber_start', orbitN[0]))

        ymdlistvariable = ymdList(im_coll)
        ymd_year = [el[:4] for el in ymdlistvariable]
        uniq_year = list(map(int, list(set(ymd_year))))
        uniq_year.sort()

        # ims = []
        # fns = []
        # rgns = []
        for i in range(len(uniq_year)):
            startDATE = ee.Date(str(uniq_year[i]) + '-01-01')
            endDATE = ee.Date(str(uniq_year[i]) + '-12-31')
            im1 = im_coll.filterDate(startDATE,endDATE).select('VH').mean().clip(aoi)
            im2 = im_coll.filterDate(startDATE,endDATE).select('VV').mean().clip(aoi)
            ims.append(im1)
            ims.append(im2)
            fns.append(str(result_path+'/'+str(uniq_year[i])+'_'+str(orbitN[0])+"_VH.tif"))
            fns.append(str(result_path+'/'+str(uniq_year[i])+'_'+str(orbitN[0])+"_VV.tif"))
            rgns.append(region)
            rgns.append(region)

        if len(year_from_2) > 0:
            for i in range(len(year_from_2)):
                startDATE = ee.Date(str(year_from_2[i]) + '-01-01')
                endDATE = ee.Date(str(year_from_2[i]) + '-12-31')
                im1 = im_coll2.filterDate(startDATE,endDATE).select('VH').mean().clip(aoi)
                im2 = im_coll2.filterDate(startDATE,endDATE).select('VV').mean().clip(aoi)
                ims.append(im1)
                ims.append(im2)
                orbitN = im_coll2.aggregate_array('relativeOrbitNumber_start').getInfo() 
                fns.append(str(result_path+'/'+str(year_from_2[i])+'_'+str(orbitN[0])+"_VH.tif"))
                fns.append(str(result_path+'/'+str(year_from_2[i])+'_'+str(orbitN[0])+"_VV.tif"))
                rgns.append(region)
                rgns.append(region)

        if len(year_from_1) > 0:
            for i in range(len(year_from_1)):
                startDATE = ee.Date(str(year_from_1[i]) + '-01-01')
                endDATE = ee.Date(str(year_from_1[i]) + '-12-31')
                im1 = im_coll1.filterDate(startDATE,endDATE).select('VH').mean().clip(aoi)
                im2 = im_coll1.filterDate(startDATE,endDATE).select('VV').mean().clip(aoi)
                ims.append(im1)
                ims.append(im2)
                orbitN = im_coll1.aggregate_array('relativeOrbitNumber_start').getInfo() 
                fns.append(str(result_path+'/'+str(year_from_1[i])+'_'+str(orbitN[0])+"_VH.tif"))
                fns.append(str(result_path+'/'+str(year_from_1[i])+'_'+str(orbitN[0])+"_VV.tif"))
                rgns.append(region)
                rgns.append(region)

    download_parallel(zip(ims, fns, rgns))
# %%
