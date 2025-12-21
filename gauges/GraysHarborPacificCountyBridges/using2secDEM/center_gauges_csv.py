"""
Defines a function convert that does this following:
    Read in an excel spreadsheet with gauge locations in deg,min,sec
    convert to decimal degrees
    add 1000 to gauge numbers
    center the gauges in 1/3" finite volume cells
    write out a .csv file with columns gaugeno,longitude,latitude
    create a .kml file with the original locations and the recentered
        locations shown as red and yellow pins
"""

from pylab import *
import pandas as pd
import re,os,sys

CHT = os.environ['CHT']
sys.path.insert(0,os.path.join(CHT,'common_code'))
from center_points import adjust

# mesh spacing and a cell edge on a 1/3" grid, for centering gauges:
# (based on the domain and grid resolutions specified in setrun.py)
dx = dy = 1./(3*3600.)
xshift = 0.5*dx
yshift = 0.5*dy

#xshift = 0.
#yshift = 0.

x_edge = -137.0 - xshift
y_edge = 47. - yshift
center_gauges = False
orig_gauges_in_kml = True

def convert(fname):
    print('\n============\nConverting %s to kml' % fname)
    bridges = pd.read_csv(fname)

    bridges1 = bridges[bridges['Flagged for Inundation'] == 1]
    bridges0 = bridges[bridges['Flagged for Inundation'] == 0]

    fname_kml = os.path.splitext(fname)[0] + '_flags.kml'
    with open(fname_kml,'w') as kml_file:
        kml_file.write("""<?xml version="1.0" encoding="UTF-8"?>
            <kml xmlns="http://www.opengis.net/kml/2.2"
            xmlns:gx="http://www.google.com/kml/ext/2.2">
            <Document><name>%s</name>

            <Style id="Red">
            <IconStyle><color>FF0000FF</color></IconStyle>
            </Style>

            <Style id="Yellow">
            <IconStyle><color>FF00FFFF</color></IconStyle>
            </Style>
            """  % fname)

        for bridge_id,y,x,flag in zip(bridges['ID'],bridges['Latitude'],
                              bridges['Longitude'],
                              bridges['Flagged for Inundation']):
            print('\n------\n', bridge_id,y,x)
            gaugeno = bridge_id

            #if orig_gauges_in_kml:
            if flag:
                kml_file.write("""
                    <Placemark><name>%s</name>
                    <description>%.5f, %.5f</description>
                    <styleUrl>#%s</styleUrl>
                    <Point>
                    <coordinates>
                    %.9f, %.9f, 0.0
                    </coordinates>
                    </Point>
                    </Placemark>
                    """ % (gaugeno,x,y,'Red',x,y))


            else:
                #if center_gauges:
                # adjust so x,y are at cell centers of 1/3" grid:
                x = adjust(x, x_edge, dx, verbose=True)
                y = adjust(y, y_edge, dy, verbose=True)
                kml_file.write("""
                    <Placemark><name>%s</name>
                    <description>%.5f, %.5f</description>
                    <styleUrl>#%s</styleUrl>
                    <Point>
                    <coordinates>
                    %.9f, %.9f, 0.0
                    </coordinates>
                    </Point>
                    </Placemark>
                    """ % (gaugeno,x,y,'Yellow',x,y))

            print('gaugeno = %4i  y = %.7f  x = %.7f' % (gaugeno,y,x))


        kml_file.write("\n</Document>\n</kml>")

    print('\nCreated %s' % fname_kml)

if __name__ == '__main__':

    fname = 'GraysHarborPacificCounty.csv'
    convert(fname)
