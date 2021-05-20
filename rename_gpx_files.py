#!/usr/bin/env python
# -*- coding: latin 1 -*-


"""

tar pakkaus kompressiolla

tar -cvzf foo.tar.gz *

purku

gunzip foo.tar.gz
tar -xvf foo.tar



"""


import pdb
import time

import os
import glob
import shutil
import pdb
import sys
import xml.etree.ElementTree as ET

import png  # pip install pypng
import io
import gzip
import datetime

# pip install fitparse
import fitparse

gpx_files_folder = '/home/antti/tmp/dumpperi'
target_folder = '/home/antti/map_tiles/renamed_sports_tracker_activities'
# resoluutiotavoite


def parse_fit_file(filename):


    if filename.endswith('.gz'):
        f = gzip.open(filename, 'rb')
        buf = io.BytesIO(f.read())
        buf.seek(0)
        fitfile = fitparse.FitFile(buf)
    else:
        fitfile = fitparse.FifFile(filename)

    record = fitfile.get_messages("record").__next__()
    first_timestamp = record.get('timestamp').value

    
    return first_timestamp
    

def parse_gpx_file(filename):


    ET.register_namespace('', "http://www.topografix.com/GPX/1/1")
    ET.register_namespace('', "http://www.topografix.com/GPX/1/0")
    ET.register_namespace('gpxx', "http://www.garmin.com/xmlschemas/GpxExtensions/v3")
    ET.register_namespace('gpxtrkx',"http://www.garmin.com/xmlschemas/TrackStatsExtension/v1" )
    ET.register_namespace('wptx1', "http://www.garmin.com/xmlschemas/WaypointExtension/v1")
    ET.register_namespace('gpxtpx', "http://www.garmin.com/xmlschemas/TrackPointExtension/v1")

    # Strava dumpista osta gpx tiedostoista gzip:llä pakattuja,joka puretaan ensin file objektiin
  
    if filename.endswith('.gz'):
        f = gzip.open(filename, 'rb')
        buf = io.BytesIO(f.read())
        buf.seek(0)
        element = ET.parse(buf)
    else:
        # parse file content to list of trees
        element=ET.parse(filename)
    
    
    gpx_version = element._root.attrib.get('version')[-1]
    

    # ideoita
    # - kaivele kaikista fileistä träkit ('{http://www.topografix.com/GPX/1/1}trk')
    # - kaivele träkeistä segmentit ('{http://www.topografix.com/GPX/1/1}trkseg')
    # - kaivele segmenteistä track pointit ('{http://www.topografix.com/GPX/1/1}trkpt')
    #
    # esim: trees[0].getroot().findall('{http://www.topografix.com/GPX/1/1}trk')[0]
    # .findall('{http://www.topografix.com/GPX/1/1}trkseg')[0]
    # .findall('{http://www.topografix.com/GPX/1/1}trkpt')

    # root - trk - seg - trkpt

    # Tarkoituksena on tunkea kaikki track pointit yhteen ja samaan segmenttiin yhteen ainoaan träkkiin


    
    
    point_elements=[]
    tracks=element.getroot().findall('{http://www.topografix.com/GPX/1/%s}trk'%gpx_version)
    for track in tracks:
        segments=track.findall('{http://www.topografix.com/GPX/1/%s}trkseg'%gpx_version)
        for segment in segments:
            points=segment.findall('{http://www.topografix.com/GPX/1/%s}trkpt'%gpx_version)
            point_elements.extend(points)
            break
        break
                        
    # pisteen aikaleiman saa pihalle oheiselle, voi käyttää aikafiltteröintii
    date_str= points[0].findall('{http://www.topografix.com/GPX/1/%s}time'%gpx_version)[0].text

    # aika voi olla desimaalisekunneilla, rapsitaan pois ('2020-11-29T09:20:21.400Z')

    if date_str.find('.') > -1:
        date_str = date_str.split('.')[0]+'Z'

    
    first_timestamp = datetime.datetime.strptime(date_str, '%Y-%m-%dT%H:%M:%SZ')
    
    return first_timestamp



all_gpx = glob.glob(gpx_files_folder + '/*.gpx')
all_gpx_gz = glob.glob(gpx_files_folder + '/*.gpx.gz')
all_fit = glob.glob(gpx_files_folder + '/*.fit.gz')


all_files = all_gpx + all_gpx_gz + all_fit


for gpx_file in all_files:
    #print("Parsing file: {}".format(gpx_file))

    if gpx_file.endswith('gpx') or gpx_file.endswith('gpx.gz'):
        try:
            first_timestamp = parse_gpx_file(gpx_file)
        except:
            print("-------FAILED GPX FILE: {}".format(gpx_file))
            continue

    elif gpx_file.endswith('fit.gz'):
        try:
            first_timestamp = parse_fit_file(gpx_file)
        except:
            print("-------FAILED FIT FILE: {}".format(gpx_file))
            continue

    else:
        print('Ei tunnettu tiedostopaate, file: {}'.format(gpx_file))
        continue
        
        
    
    timestring = first_timestamp.strftime('%Y%m%d_%H%M')
    
    #pdb.set_trace()   
    fileending = gpx_file[gpx_file.find('.'):]

    target_name = os.path.join(target_folder,timestring + fileending)
    
    # copy file to new name
    print("{} ==> {}".format(gpx_file, target_name))
    shutil.copy(gpx_file, target_name)


