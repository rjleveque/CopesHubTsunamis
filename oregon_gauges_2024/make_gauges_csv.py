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
import re
from center_points import adjust

# mesh spacing and a cell edge on a 1/3" grid, for centering gauges:
# (based on the domain and grid resolutions specified in setrun.py)
dx = dy = 1./(3*3600.)
#xshift = 0.5*dx
#yshift = 0.5*dy
xshift = 0.
yshift = 0.
x_edge = -137.0 - xshift
y_edge = 47. - yshift
center_gauges = True

def convert(fname):
    print('\n============\nConverting %s to csv and kml' % fname)
    df = pd.read_excel(fname + '.xlsx')
    fname_csv = fname + '.csv'
    fname_kml = fname + '.kml'
    with open(fname_csv,'w') as csv_file, open(fname_kml,'w') as kml_file:
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

        for vg,yll,xll in zip(df['VG'],df['Lat'],df['Long']):
            print('\n------\n', vg,yll,xll)
            gaugeno = vg + 1000
            tokens = re.split("""[°'"T]""",yll)
            y = float(tokens[0]) + float(tokens[1])/60. + float(tokens[2])/3600.
            tokens = re.split("""[°'"T]""",xll)
            x = float(tokens[0]) + float(tokens[1])/60. + float(tokens[2])/3600.
            x = -x

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
            
            if center_gauges:
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
                    """ % ('',x,y,'Yellow',x,y))
            
            print('gaugeno = %4i  y = %.7f  x = %.7f' % (gaugeno,y,x))

            csv_file.write('%5i,  %.7f,  %.7f\n' % (gaugeno,x,y))

        kml_file.write("\n</Document>\n</kml>")
        
    print('\nCreated %s and %s' % (fname_csv,fname_kml))

if __name__ == '__main__':

    #convert('VGListAstoria')
    
    if 1:
        for fname in ['VGListCoosBay','VGListAstoria',
                      'VGListNewport','VGListSeaside']:
            convert(fname)
