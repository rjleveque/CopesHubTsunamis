from pylab import *


gauge_list = []

lines = open('gnss_list.csv').readlines()

for line in lines[1:]:
    tokens = line.split(',')
    #print('tokens: ',tokens)
    try:
        name = tokens[0].strip()
        x = float(tokens[1])
        y = float(tokens[2])
        gauge = [name, x, y]
        gauge_list.append(gauge)
    except:
        pass

print('gauge_list = ',gauge_list)
    
#red_gauges = ['MIZU','USUD']
red_gauges = [g[0] for g in gauge_list]
    
kml_fname = 'gnss_csz62.kml'

with open(kml_fname,'w') as kml_file:
       
    kml_file.write("""<?xml version="1.0" encoding="UTF-8"?>
    <kml xmlns="http://www.opengis.net/kml/2.2"
    xmlns:gx="http://www.google.com/kml/ext/2.2">
    <Document><name>gnss_CSZ62</name>
    
    <Style id="Red">
    <IconStyle><color>FF0000FF</color></IconStyle>
    </Style>
    
    <Style id="Yellow">
    <IconStyle><color>FF00FFFF</color></IconStyle>
    </Style>
    """)
        
    for gauge in gauge_list:
        name = gauge[0]
        x = gauge[1]
        y = gauge[2]
        if name in red_gauges:
            style_id = 'Red'
        else:
            style_id = 'Yellow'
            
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
        """ % (name,x,y,style_id,x,y))
    
    kml_file.write("\n</Document>\n</kml>")
        
print('Created ', kml_fname)
