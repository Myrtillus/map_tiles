#!/usr/bin/env python
# -*- coding: latin 1 -*-


import scipy.sparse
from scipy.sparse import csr_matrix
import scipy.signal

import math
import numpy
import pdb
import time

import os
import glob
import pdb
import sys
import xml.etree.ElementTree as ET

import png  # pip install pypng
import io
import gzip
import datetime

# pip install fitparse
import fitparse

tile_folder_path = '/home/ubuntu/map_tiles/tiles'



gpx_files_folder ='/home/ubuntu/map_tiles/activities'

# resoluutiotavoite
resolution_ew_deg = 0.0001475/3
resolution_ns_deg = 0.0000615/3

# laskettavan alueen rajat
west_edge = 23.3
east_edge = 24.3
north_edge = 61.7
south_edge = 61.3

# laskettavat zoomi levelit tiilille
zoom_levels=[11,12,13,14,15,16,17,18]


# Aikaväli, miltä gpx jäljet hyväksytään mukaan
timerange_start = datetime.datetime(2000,1,1)
timerange_end = datetime.datetime(2022,1,1)


# Aikaikkunat värjäyksille
# 0-6 kk keltainen
# 6-12 kk oranssi
# 12 - 18 kk tumman punainen
# 18+ kk harmaa


track_age_limits = [
    [6*30, 7],
    [12*30, 6],
    [18*30, 5],
]

palette=[
        (0,0,0,0),        # läpinäkyvä tausta (0)
        (0,0,0,0),        # varaus
        (0,0,0,0),        # varaus
        (0,0,0,0),        # varaus
        (177, 178, 179, 150),   # hylätyt trackit (4)
        (161, 35, 71, 0xff),   # kolmos luokka (5)
        (240, 132, 0, 0xff),   # kakkosluokka (6)
        (236, 240, 0, 0xff),   # uusimmat (7)
        ]   

#################################################################



# convert latlon to tilenumbers
def deg2num(lat_deg, lon_deg, zoom):
  lat_rad = math.radians(lat_deg)
  n = 2.0 ** zoom
  xtile = int((lon_deg + 180.0) / 360.0 * n)
  ytile = int((1.0 - math.asinh(math.tan(lat_rad)) / math.pi) / 2.0 * n)
  return (xtile, ytile)

# convert tile numbers to lat lon  
def num2deg(xtile, ytile, zoom):
  n = 2.0 ** zoom
  lon_deg = xtile / n * 360.0 - 180.0
  lat_rad = math.atan(math.sinh(math.pi * (1 - 2 * ytile / n)))
  lat_deg = math.degrees(lat_rad)
  return (lat_deg, lon_deg)

def scale(im, nR, nC):
  nR0 = len(im)     # source number of rows 
  nC0 = len(im[0])  # source number of columns 
  return [[ im[int(nR0 * r / nR)][int(nC0 * c / nC)]  for c in range(nC)] for r in range(nR)]


# ----------------------------------

def parse_fit_file(filename):


    # https://www.fitfileviewer.com/
    # fit fileen koordinaatistomuutos
    # https://gis.stackexchange.com/questions/371656/garmin-fit-coodinate-system



    if filename.endswith('.gz'):
        f = gzip.open(filename, 'rb')
        buf = io.BytesIO(f.read())
        buf.seek(0)
        fitfile = fitparse.FitFile(buf)
    else:
        fitfile = fitparse.FifFile(filename)
    
    gps_points=[]

    first_record = True    
    for record in fitfile.get_messages("record"):

        if first_record:
            first_record = False
            first_timestamp = record.get('timestamp').value
            if first_timestamp < timerange_start or first_timestamp > timerange_end:
                print('Track is outside desired timerange')
                return (None, None)

        try:
            lat = record.get('position_lat').value/11930465.0
            lon = record.get('position_long').value/11930465.0
        except:
            pass
            continue

        gps_points.append((float(lat), float(lon)))

    if gps_points == []:
        print("Could not obtain any points")
        return (None, None)
    return (gps_points, datetime.datetime.now() - first_timestamp)        
    

# ----------------------------------


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

    # tsekataan, että track halutulta aikaväliltä

    # pisteen aikaleiman saa pihalle oheiselle, voi käyttää aikafiltteröintii
    date_str= points[0].findall('{http://www.topografix.com/GPX/1/%s}time'%gpx_version)[0].text

    # aika voi olla desimaalisekunneilla, rapsitaan pois ('2020-11-29T09:20:21.400Z')

    if date_str.find('.') > -1:
        date_str = date_str.split('.')[0]+'Z'

    
    first_timestamp = datetime.datetime.strptime(date_str, '%Y-%m-%dT%H:%M:%SZ')
    if first_timestamp < timerange_start or first_timestamp > timerange_end:
        print('Track is outside desired timerange')
        return (None, None)


    # x ja y eli pituus ja leveyspiirit
    gps_points = [(float(y.get('lat')),float(y.get('lon'))) for y in [x.attrib for x in point_elements]]
    
    return (gps_points, datetime.datetime.now() - first_timestamp)


# ----------------------------------


def latlon2slot(latlon_pair, ns_slot_count, ew_slot_count):
    # konvertoi lat lon koordinaatit varastomatriisin slotti-indekseihin
	# slottinumero kulkee etelästä pohjoiseen ja lännestä itään alla olevassa
	# varsinaisessa matriisissa ns indeksi kulkee ylhäältä alas
	ns_slot = int((latlon_pair[0]-south_edge) / resolution_ns_deg)
	ew_slot = int((latlon_pair[1]-west_edge) / resolution_ew_deg)

	return(ns_slot,ew_slot)




# ----------------------------------
# ----------------------------------
# ----------------------------------
# ----------------------------------




span_ew = east_edge - west_edge
span_ns = north_edge - south_edge

ew_slot_count = int(float(span_ew) / float(resolution_ew_deg))
ns_slot_count = int(float(span_ns) / float(resolution_ns_deg))

# Harva scipy matriisi, johon kasataan pisteinformaatiota

slot_matrix = csr_matrix((ns_slot_count, ew_slot_count), dtype=int)

# Keräillään gpx ja fit tiedostojen nimet polkuineen

all_gpx = glob.glob(gpx_files_folder + '/*.gpx')
all_gpx_gz = glob.glob(gpx_files_folder + '/*.gpx.gz')
all_fit = glob.glob(gpx_files_folder + '/*.fit.gz')

all_files = all_gpx + all_gpx_gz + all_fit


# Looppaillaan tiedostot lävitse

some_data_found = False

for gpx_file in all_files:
    print("Parsing file: {}".format(gpx_file))

    if gpx_file.endswith('gpx') or gpx_file.endswith('gpx.gz'):
        try:
            (latlon_pairs, track_age) = parse_gpx_file(gpx_file)
        except:
            print("-------FAILED GPX FILE: {}".format(gpx_file))
            continue

    elif gpx_file.endswith('fit.gz'):
        try:
            (latlon_pairs, track_age) = parse_fit_file(gpx_file)
        except:
            print("-------FAILED FIT FILE: {}".format(gpx_file))
            continue

    else:
        print('Ei tunnettu tiedostopaate, file: {}'.format(gpx_file))
        continue
        
    # parserit palauttaa Nonen, jos ei dataa tai halutun aikavälin ulkopuolella  
    if latlon_pairs is None:
        continue

    some_data_found = True


    # Ripotellaan gpx tiedoston pisteet hilamatriisiin

    slot_indeces = []		
    for item in latlon_pairs:
	    if float(item[0]) < south_edge or float(item[0]) > north_edge:
		    continue
	    if float(item[1]) < west_edge or float(item[1]) > east_edge:
		    continue
	    slot_indeces.append(latlon2slot([float(x) for x in item], ns_slot_count, ew_slot_count))  #latlon2slot return (ns_slot,ew_slot)
	
	# otetaan pihalle indeksit, jotka menevät matriisin ulkopuolelle
    remove=[]
    for item in slot_indeces:
        if item[0] <0:
            remove.append(item)
            continue
        if item[0]>ns_slot_count -1:
            remove.append(item)
            continue
        if item[1] < 0:
            remove.append(item)
            continue
        if item[1] > ew_slot_count -1:
            remove.append(item)
            continue

    slot_indeces = [x for x in slot_indeces if x not in remove]	
	    
    if slot_indeces == []:
	    #print('gps points outside the given area')
	    continue

    # store to matrix 
    # matriisi on pystysuunnassa ylösalaisin. Pienet slotinumerot edustavat eteläistä aluetta ja suuret numerot pohjoista.
    # matriisin vasenlaita länsilaitaa ja oikea itälaitaa.

    # trackille palautui ikä, joten haarukoidaan, käytettävä väripaletin indeksi
    paletti_index_for_track = [x[1] for x in track_age_limits if track_age < datetime.timedelta(x[0])]
    if paletti_index_for_track == []:
        paletti_index_for_track = 4
    else:
        paletti_index_for_track = max(paletti_index_for_track)


    # matriisin shapessa eka on rivien määrä (NS suunta) ja toinen on leveys
    # slicing koordinaateista eka osoittaa riviä ja toinen saraketta

    (I,J) = zip(*slot_indeces)    # I = ns slot indeksi, J = ew indeksi

    # poista duplikaatti parit, jotta eivät summaudu matriisiin. Näin yhteen hilamatriisin elementtiin
    # tulee vain yksi arvo, joka edustaa pisteen ikää, vaikka gpx fileessä pisteitä olisi enemmänkin

    out = list(set(zip(I,J)))

    # tuupataan temppi matriisiin

    (I,J) = zip(*out)
    V = numpy.ones((len(J)))*paletti_index_for_track 
    tmp_slot_matrix = csr_matrix((V, (I, J)), shape = (ns_slot_count, ew_slot_count), dtype=numpy.int8)

    # Yhdistetään tmp matriisin ja varsinaisen varastomatriin tiedot siten, että elementeittäin valitaan
    # aina suurempi arvo. Tehtävä boolean kikkailun kautta, koska suoraa metodia ei ollut tiedossa
    # harvalle matriisille.

    boolean = slot_matrix < tmp_slot_matrix
    slot_matrix = slot_matrix - slot_matrix.multiply(boolean) + tmp_slot_matrix.multiply(boolean)




if not some_data_found:
    print("Ei mitään gps dataa renderöitäväksi")
    sys.exit()

####################################
# Aloitetaan tiilien renderöinti
####################################


for zoom_level in zoom_levels:
    print("Lasketaan taso: {}".format(zoom_level))
    if not os.path.isdir(os.path.join(tile_folder_path,str(zoom_level))):
        os.mkdir(os.path.join(tile_folder_path,str(zoom_level)))

    # laske tiilien numerot reunoilla valitulla zoomitasolla

    # vasen ylälaita
    (xtile_start, ytile_start) = deg2num(north_edge, west_edge, zoom_level)
    
    # oikea alalaita
    (xtile_end, ytile_end) = deg2num(south_edge, east_edge, zoom_level)

    # looppaile tiilet x-suunnassa
    for tile_ykoord in range(ytile_start, ytile_end+1):
        # looppaile tiilet y-suunnassa
        
        for tile_xkoord in range(xtile_start,xtile_end+1):
            # laske tiilen alueen koordinaattireunat
            #def num2deg(xtile, ytile, zoom)
            #  return (lat_deg, lon_deg)
            
            (lat_north_edge, lon_west_edge) = num2deg(tile_xkoord,tile_ykoord, zoom_level)
            (lat_south_edge, lon_east_edge) = num2deg(tile_xkoord+1,tile_ykoord+1, zoom_level)

            # laske indeksit slot matriisiin, josta löytyvät kertymät
            #def latlon2slot(latlon_pair):
            #	return(ns_slot,ew_slot)
            (ns_start_slot, ew_start_slot) = latlon2slot((lat_north_edge, lon_west_edge), ns_slot_count, ew_slot_count)
            (ns_end_slot, ew_end_slot) = latlon2slot((lat_south_edge, lon_east_edge), ns_slot_count, ew_slot_count)

           
            # Aluerajat eivät osu tiilien saumakohtiin
            # Alueraja osuus esim pohjoisreunassa etelämmäksi kuin tiilen reuna. Kun taas lasketaan
            # toisinpäin tiilinumerosta sloteihin, niin tule negatiivisia slot numeroita.
            # Harvasta datamatriisista on siivutettava tavaraa olemassa olevalta alueelta
            # yli ylitse menevien indeksien verran on lisättävä nollarivejä tai nolla sarakkeita.
            
            # x-padding haettavalle matriisille
            if ew_start_slot < 0:
                left_padding = abs(ew_start_slot)
                ew_start_slot = 0
            else:
                left_padding = 0
                
            if ew_end_slot > slot_matrix.shape[1]:
                right_padding = ew_end_slot - slot_matrix.shape[1]
                ew_end_slot = slot_matrix.shape[1]
            else:
                right_padding = 0

            # y-padding haettavalle matriisille
            if ns_end_slot < 0:
                bottom_padding = abs(ns_end_slot)
                ns_end_slot = 0
            else:
                bottom_padding = 0
                
            if ns_start_slot > slot_matrix.shape[0]:
                top_padding = ns_start_slot - slot_matrix.shape[0]
                ns_start_slot = slot_matrix.shape[0]
            else:
                top_padding = 0

            # haetaan matriisi slice, keikautettu haku pystysuunnassa
            pickup_matrix = slot_matrix[ns_end_slot:ns_start_slot,:][:,ew_start_slot:ew_end_slot].todense()

            # padataan matriisi
            pickup_matrix = numpy.pad(pickup_matrix,((bottom_padding, top_padding),(left_padding, right_padding)),mode='constant', constant_values = 0)
            
            # flipataan y-suunnassa
            pickup_matrix = numpy.flipud(pickup_matrix)

            # tarkistetaan, että onko matriisissa mitään plotattavaa
            if numpy.sum(pickup_matrix) > 0:

                        
                # Jos yksittäiset pisteet sellaisenaan plotataan tiileen, niin niitä ei tahdo erottaa.
                # Näin ollen pisteen aluetta levennetään. Käytännössä se tehdään
                # matriisin kaksi suuntaiselle konvoluutiolla, jossa neliön muotoisen yksikkömatriisin
                # koon avulla voidaan säätää pisteen leveyttä tiilessä.
            
                multipliers ={
                9:60,
                10:30,
                11:15,
                12:10,
                13:7,
                14:4,
                15:3,
                16:2,
                17:1,
                18:1,
                }
                
                multiplier = multipliers.get(zoom_level)

                if multiplier:
                #if False:       


                    # Ennen konvoluutiota Matriisin elementeissä on kokonaislukuja, jotka kuvaavat
                    # pisteen ikää ja osoittavat tiettyyn väripaletin elementtiin.
                    # Tästä matriisista on ensin erotettava eri kokonaisluvut omiin matriiseihinsa
                    # ja merkattava niissä lukuja sisältävä elementit ykkösellä.
                    # Tämän jälkeen jokaiselle matriisille ajetaan erikseen konvoluutio ja korvataan
                    # konvoluution jälkeen nollasta poikkeavat elementit alkuperäisellä kokonaisluvulla
                    # Lopuksi vielä yhdistetään eri ikäluokkien matriisit yhteen siten, että valitaan aina
                    # elementeittäin suurin kokonaisluku, jolloin uusin jää voimaan.


                    # Haetaan indeksi eri ikäluokan elementteihin
                    index4 = numpy.where(pickup_matrix == 4)
                    index5 = numpy.where(pickup_matrix == 5)
                    index6 = numpy.where(pickup_matrix == 6)
                    index7 = numpy.where(pickup_matrix == 7)


                    # muodostetaan apumatriisit, joissa ikäluokittain ykkönen kohdissa, joissa on oli vastaava ikäluku
                    pickup_matrix_level4 = numpy.zeros(pickup_matrix.shape)
                    pickup_matrix_level4[index4] = 1

                    pickup_matrix_level5 = numpy.zeros(pickup_matrix.shape)
                    pickup_matrix_level5[index5] = 1

                    pickup_matrix_level6 = numpy.zeros(pickup_matrix.shape)
                    pickup_matrix_level6[index6] = 1

                    pickup_matrix_level7 = numpy.zeros(pickup_matrix.shape)
                    pickup_matrix_level7[index7] = 1


                    # pyöräytetään konvoluutio

                    conv_matrix = numpy.ones([multiplier, multiplier])

                    pickup_matrix_level4 = scipy.signal.convolve2d(pickup_matrix_level4, conv_matrix,mode="same")
                    pickup_matrix_level5 = scipy.signal.convolve2d(pickup_matrix_level5, conv_matrix,mode="same")
                    pickup_matrix_level6 = scipy.signal.convolve2d(pickup_matrix_level6, conv_matrix,mode="same")
                    pickup_matrix_level7 = scipy.signal.convolve2d(pickup_matrix_level7, conv_matrix,mode="same")


                    # meni flotareiksi, takaisin int8
                    # rajataan 0-1 välille ja kerrotaan levelillä
                    pickup_matrix_level4 = numpy.array(pickup_matrix_level4,dtype=int)
                    pickup_matrix_level4 = numpy.clip(pickup_matrix_level4,0,1) * 4

                    pickup_matrix_level5 = numpy.array(pickup_matrix_level5,dtype=int)
                    pickup_matrix_level5 = numpy.clip(pickup_matrix_level5,0,1) * 5

                    pickup_matrix_level6 = numpy.array(pickup_matrix_level6,dtype=int)
                    pickup_matrix_level6 = numpy.clip(pickup_matrix_level6,0,1) * 6

                    pickup_matrix_level7 = numpy.array(pickup_matrix_level7,dtype=int)
                    pickup_matrix_level7 = numpy.clip(pickup_matrix_level7,0,1) * 7


                    # otetaan eri level matriiseista elementeittäin maksimi renderöintiä varten

                    pickup_matrix = numpy.maximum(pickup_matrix_level4, pickup_matrix_level5)
                    pickup_matrix = numpy.maximum(pickup_matrix, pickup_matrix_level6)
                    pickup_matrix = numpy.maximum(pickup_matrix, pickup_matrix_level7)


                # varmistaan, että ollaa integerinä
                pickup_matrix = numpy.array(pickup_matrix,dtype=int)

                # resize the image matrix to 256 x 256 pixels
                pickup_matrix = scale(pickup_matrix, 256, 256)
                
                
                # talletaan oikeaan hakemistoon /zoom/x/y.png
                filepath = os.path.join(tile_folder_path,str(zoom_level), str(tile_xkoord))
                if not os.path.isdir(filepath):
                    os.mkdir(filepath)
                filename = "{}.png".format(tile_ykoord)
                
                print("saving: {}".format(os.path.join(filepath,filename)))
                
                with open(os.path.join(filepath,filename), 'wb') as f:
                    w = png.Writer(256,256, palette=palette, bitdepth=8)   
                    w.write(f, pickup_matrix)
                
                continue
                
                

