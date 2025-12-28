from pylab import *
import os

not_done = [99]
for i in range(550):    #i goes from 0 to 449
    not_done=not_done+[99]

fal=[False]
for i in range(550):    #i goes from 0 to 449
    fal=fal+[False]

CSZL1=fal; L1up_x1pt095S=fal; 
WGSDOGAMI_2pt5=fal; WGSDOGAMI_2pt5_x1pt03=fal; 

CSZL1=array(CSZL1); L1up_x1pt095S=array(L1up_x1pt095S); 
WGSDOGAMI_2pt5=array(WGSDOGAMI_2pt5); WGSDOGAMI_2pt5_x1pt03=array(WGSDOGAMI_2pt5_x1pt03);

#CSZL1 compliance at gauges 1-169:
##  This would be at indices 0-168;
CSZL1[0:169]=True

#L1up_x1pt095S compliance at gauges 1-185:
##  This would be at indices 0-184;
L1up_x1pt095S[0:185]=True

#WGS-DOGAMI compliance at gauges 1-142; 163-184; 202-306; 379-514;
##  This would be at indices 0-141; 162-183; 201-305; 378-513;
WGSDOGAMI_2pt5[0:142]=True; WGSDOGAMI_2pt5[162:184]=True; WGSDOGAMI_2pt5[201:306]=True;
WGSDOGAMI_2pt5[378:514]=True;

#WGS-DOGAMI_x1.03 compliance at gauges 1-306; 379-514;
##  This would be at indices 0-305; 378-513;
WGSDOGAMI_2pt5_x1pt03[0:306]=True;
WGSDOGAMI_2pt5_x1pt03[378:514]=True;

#### Pick up the Latitudes in array lat below
root_dir     = os.environ['CHT']
asce_powell_dir = root_dir + '/info/asce_powell_gauges'
pathname = asce_powell_dir + '/7-28_loyce_north-to-south-32km.txt'
ID,lon,lat = np.loadtxt(pathname, delimiter=',',\
                                  skiprows=1,usecols=(0,1,2),unpack=True)

print ('Gauge-no,   Latitude,  CSZ_L1,  L1up_1.095S,  WGSD2.5,  WGSD2.5_1.03')
for i in range(551):
    ig = i+1
    print('%8i, %10.6f, %7s, %12s, %8s, %13s' %(ig,lat[i],CSZL1[i],\
           L1up_x1pt095S[i],WGSDOGAMI_2pt5[i],WGSDOGAMI_2pt5_x1pt03[i]) )
