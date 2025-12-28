from pylab import *
import os

not_done = [99]
for i in range(550):    #i goes from 0 to 449
    not_done=not_done+[99]

fal=[False]
for i in range(550):    #i goes from 0 to 449
    fal=fal+[False]

BLM10=fal; BLD10=fal; BRD10=fal; BLMD=fal; FLD10=fal; FRD10=fal;

#Not upping BLMD yet.
BLM10_2pt0=fal; BLD10_1pt5=fal; BRD10_1pt5=fal; FLD10_1pt5=fal; FRD10_1pt5=fal;  

BLM10=array(BLM10); BLD10=array(BLD10); BRD10=array(BRD10); BLMD=array(BLMD); FLD10=array(FLD10); FRD10=array(FRD10);
BRD10_1pt5=array(BRD10_1pt5); BLM10_2pt0=array(BLM10_2pt0); BLD10_1pt5=array(BLD10_1pt5); FLD10_1pt5=array(FLD10_1pt5);
FRD10_1pt5=array(FRD10_1pt5);

########## Original Audrey time dependent ruptures
#BLMD is buried-locking-mur13-deep; BLM10 is buried-locking-str10-middle;
#BLD10 is buried-locking_str10-deep; BRD10 is buried-random_str10-deep;
#FLD10 is ft-locking-str10-deep; FRD10 is ft-random-str10-deep; 

#BLM10 had no compliance gauges, so leave that alone as False for all
#BLD10 had no compliance gauges, so leave that alone as False for all
#BLMD  had no compliance gauges, so leave that alone as False for all

#BRD10 had compliance at gauges 88-89; 366-389; 491-498;
##  This would be indices       87-88; 365-388; 490-497;
BRD10[87:89]=True; BRD10[365:389]=True; BRD10[490:498]=True;

#FLD10 had compliance at gauges 481-498;
##  This would be indices 480-497;
FLD10[480:498]=True;

#FRD10 had compliance at gauges 493-498;
##  This would be indices 492-498;
FRD10[492:498]=True;
###########

####  Info for several communities ####
##  Westport:  Latitude: 46.908, roughly Gauge 145 at latitude 46.906550, Gauge 144 is at latitude 46.920830 
#              BLD10_1pt5 is compliant there
#              BRD10_1pt5 is compliant there
##

########## Amplified Audrey gauges 
#BLM10_2pt0 compliance at gauges 1-122; 349-391; 480-504;
##  This would be indices 0-121; 348-390; 479-503;
BLM10_2pt0[0:122]=True; BLM10_2pt0[348:391]=True; BLM10_2pt0[479:504]=True; 

#BLD10_1pt5 compliance at gauges 70-103; 138-272; 349-389; 489-498;
##  This would be indices        69-102; 137-271; 348-388; 488-497;
BLD10_1pt5[69:103]=True; BLD10_1pt5[137:272]=True; BLD10_1pt5[348:389]=True; 
BLD10_1pt5[488:498]=True;

#BRD10_1pt5 had compliance at gauges 69-163; 351-440; 486-551;
##  This would be indices            68-162; 350-439; 485-550;
BRD10_1pt5[68:163]=True; BRD10_1pt5[350:440]=True; BRD10_1pt5[485:551]=True;

#FLD10_1pt5 compliance at gauges 207-274; 336-389; 466-551;
##  This would be indices  206-273; 335-388; 465-550;   
FLD10_1pt5[206:274]=True; FLD10_1pt5[335:389]=True; FLD10_1pt5[465:551]=True; 

#FRD10_1pt5 compliance at gauges 336-383; 490-551;
##  This would be indices 335-382; 489-550;
FRD10_1pt5[335:383]=True; FRD10_1pt5[489:551]=True;

#### Pick up the Latitudes in array lat below
root_dir     = os.environ['CHT']
asce_powell_dir = root_dir + '/info/asce_powell_gauges'
pathname = asce_powell_dir + '/7-28_loyce_north-to-south-32km.txt'
ID,lon,lat = np.loadtxt(pathname, delimiter=',',\
                                  skiprows=1,usecols=(0,1,2),unpack=True)

print ('G-no,   Latitude, BL10M, BL10M, BL10D, BL10D, BR10D, BR10D,  BL13D, FL10D, FL10D, FR10D, FR10D')
print ('    ,           ,      ,  2.0 ,      ,  1.5 ,      ,  1.5 ,       ,      ,  1.5 ,      ,  1.5 ')
for i in range(551):
    ig = i+1
    print('%4i, %10.6f, %5s, %5s, %5s, %5s, %5s, %5s, %6s, %5s, %5s, %5s, %5s' %(ig,lat[i],BLM10[i],BLM10_2pt0[i],\
           BLD10[i],BLD10_1pt5[i],BRD10[i],BRD10_1pt5[i],BLMD[i],FLD10[i],FLD10_1pt5[i],FRD10[i],FRD10_1pt5[i]))
