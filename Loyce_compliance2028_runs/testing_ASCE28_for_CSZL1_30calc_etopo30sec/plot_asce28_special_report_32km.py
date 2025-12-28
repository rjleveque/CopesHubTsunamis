"""
Testing the relevant asce gauges to see if the GeoClaw
output is compliant. Will use Westport since this will be
just pre to 2028 being adopted and we want to know some things!
The gauge locations changed from before.

 Westport 2025, Gauges 98 to 143
 Westport 2028, Gauges 127 to 164 

"""

# for making plots on remote machine:
import matplotlib
matplotlib.use('Agg')  # Use an image backend

from pylab import *
import os
import copy
import clawpack.pyclaw.gauges as gauges
from clawpack.geoclaw import dtopotools
import interptools
import params  #sets root_dir, loc, event, fgmax_extent

# SETTING directories where input data is found:
rundir = os.getcwd()
root_dir        = params.root_dir
assert os.path.isdir(root_dir), '*** did not find params.root_dir = %s' \
        % root_dir

out_dir         = params.out_dir
assert os.path.isdir(out_dir), '*** did not find output directory = %s' \
        % out_dir

plot_dir        = params.plot_dir
if not os.path.isdir(plot_dir):
    os.mkdir(plot_dir)

other_figures_dir    = plot_dir + '/_other_figures'
if not os.path.isdir(other_figures_dir):
    os.mkdir(other_figures_dir)

dtopo_dir = params.dtopo_dir
source_name = params.event
dtopofile = params.dtopofile
dtopotype = 3  # format of dtopo file

### 2028
asce_powell_dir = root_dir + '/info/asce_powell_gauges'

## Directories used
print('rundir = %s' % rundir)
print('root_dir is:  ',root_dir)
print('out_dir is: ',out_dir)
print('plot_dir is: ',plot_dir)
print('other_figures_dir is: ',other_figures_dir)
print('dtopo_dir is: ',dtopo_dir)
print('dtopo file is: ',dtopofile)
print('asce_powell_dir is: ',asce_powell_dir)
print(' ')


### 2028, Loyces file below has the lines from north (starting at 1) to south (ending at 551) 
###       from north to south increasing numbering.  ID[0]=1, ID[-1]=551.  The ID has the gaugeno
###       and it is always one more than the index number into the arrays ID, lon, lat, etc.

pathname = asce_powell_dir + '/7-28_loyce_north-to-south-32km.txt'
ID,lon,lat,amp,ID_32N,ID_32S,seg32_del_lat,seg32_km,seg32_mi,seg32_amp_mean = np.loadtxt(pathname, delimiter=',',\
                                  skiprows=1,usecols=(0,1,2,3,4,5,6,7,8,9),unpack=True)
print(' Colummns in the 7-28_loyce_north-to-south-32km.txt files are: ')
print(' ID,lon,lat,amp,ID_32N,ID_32S,seg32_del_lat,seg32_km,seg32_mi,seg32_amp_mean')
print(' ID[0]=0 (northernmost) and ID[-1]=551 (southernmost)')
print(' ')

#All gauges 2028   -- gauge 1 is farthest north and 551 is farthest south
gstart=1; gend=551;

##############  Particular Structures have structure_lat, gstart32, and gend32 appended manually ############
##############  This is meant for VES where we know its latitude and it is not necessary the latitude #######
##############  of one of the 2028 ASCE gauges.  We also do a compliance graph report for each structure ####
##
##  Always +-32.2 km of structure for relevant report
delt = 32.2/111.
structure_name=[]; structure_lat=[]; 
gstart32=[]; gend32=[];
xlim1=[]; xlim2=[]; xlim3=[]; xlim4=[];
xticks1=[]; xticks2=[]; xticks3=[]; xticks4=[]; xticks5=[]; xticks6=[];

gauges_to_plot=[]  #Only plotting the Westport ones, see below
if 1:   ##New 2028 Gauge numbering
    ### For Westport 2025 VES
    structure_name.append('WP')
    latitude =46.908
    structure_lat.append(latitude)
    gstart32.append(127); gend32.append(164);   
    #gauges_to_plot=gauges_to_plot + list(range(127,165))
    ## The gaugeno is 1 more than the index into the arrays from the asce file
    #xlim1.append(asce_dict[str(125)]['lat'])
    #xlim2.append(asce_dict[str(149)]['lat'])
    #xticks1.append(asce_dict[str(127)]['lat'])
    #xticks2.append(asce_dict[str(137)]['lat'])
    #xticks3.append(asce_dict[str(147)]['lat'])
    #xlim3.append(asce_dict[str(146)]['lat'])
    #xlim4.append(asce_dict[str(166)]['lat'])
    #xticks4.append(asce_dict[str(148)]['lat'])
    #xticks5.append(asce_dict[str(156)]['lat'])
    #xticks6.append(asce_dict[str(164)]['lat'])

    xlim1.append(lat[124])
    xlim2.append(lat[148])
    xticks1.append(lat[126])
    xticks2.append(lat[136])
    xticks3.append(lat[146])
    xlim3.append(lat[145])
    xlim4.append(lat[165])
    xticks4.append(lat[147])
    xticks5.append(lat[155])
    xticks6.append(lat[163])

#### Could add more structures here if desired
num_structures = len(structure_lat)

#### Checking all the structure information to see if we got it correct.
#### Experimenting to try to pick up the gauges that are within delt
#### latitude difference to our gauges of interest. Checking here to
#### see if we got it right.
for j in range(num_structures):
    N_extent = structure_lat[j] + delt
    S_extent = structure_lat[j] - delt
    print (' ')
    print ('STRUCTURE NAME: ',structure_name[j])
    print ('Latitude of this structure was: ',structure_lat[j])
    print ('N_extent and S_extent +-32km in latitude were: %.7f  %.7f ' %(N_extent,S_extent))
    print ('The latitude of the southern most chosen asce gauge, %.1f, for this structure was %.7f '\
            %(gend32[j],lat[gend32[j]-1] ) )
            #%(gend32[j],asce_dict[str(gend32[j])]['lat'] ) )
    print ('The latitude of the northern most chosen asce gauge, %.1f, for this structure was %.7f '\
            %(gstart32[j],lat[gstart32[j]-1] ) )
            #%(gstart32[j],asce_dict[str(gstart32[j])]['lat'] ) )
    print (' ')

geoclaw_dict={}
figure(400, figsize=(8,8))
figure(500, figsize=(8,8))

# Read in deformation file
dtopo = dtopotools.DTopography()
dtopo.read(dtopofile, dtopotype)
print ('Dtopo information: ')
print ('dtopo min of dtopo.X was  %.6f' %dtopo.X.min() ) 
print ('dtopo max of dtopo.X was  %.6f' %dtopo.X.max() ) 
print ('dtopo min of dtopo.Y was  %.6f' %dtopo.Y.min() ) 
print ('dtopo max of dtopo.Y was  %.6f' %dtopo.Y.max() ) 
print ('shape of dtopo.dZ was ',    shape(dtopo.dZ) )
print ('shape of dtopo.dZ[-1,:,:] was ', shape(dtopo.dZ[-1,:,:]) )
print ('dtopo max of dtopo.dZ[-1,:,:] was %.6f' %dtopo.dZ[-1,:,:].max() )
print ('dtopo min of dtopo.dZ[-1,:,:] was %.6f' %dtopo.dZ[-1,:,:].min() )
print (' ')

num_gauges = gend-gstart+1
geo_amp = ones(num_gauges)
for gaugeno in range(gstart,gend+1):
    if 0:
        print ('----')
        print ('Geoclaw Results for ASCE Gauge no: ',gaugeno)
    gauge = gauges.GaugeSolution(gaugeno, out_dir)    #reading from out_dir, not _output
    xlong,ylat = gauge.location

    ## Determine the subsidence or uplift at this xlong, ylat
    try:
        dzi_gauge = interptools.interp(xlong,ylat,dtopo.X,dtopo.Y,dtopo.dZ[-1,:,:])
        if 0:
            print ('(long,lat) = (',xlong,',',ylat,')' )
            print ('Subsidence/uplift at the gauge was:  %s m'  % dzi_gauge )
    except:
        dzi_gauge = 0.
        print ('(long,lat) = (',xlong,',',ylat,')' )
        print ('The gauge does not overlap dtopo, so no subsidence/uplift' )

    t = gauge.t / 60.   # convert to minutes
    q = gauge.q
    h = q[0,:]
    eta = q[3,:]

    ##   B changes with time as we refine levels, but our asce gauges are each
    ##   in a given region which is turned on at time 0 and is a constant refined
    ##   region, either 6sec resolution, or in some case 30 sec resolution, so
    ##   B should be constant, and safe to pick up B[0] as the post-quake
    ##   bathymetry on the finest level.  The initial bathymetry called B0 on
    ##   the finest level at the gauge is increased or decreased by dzi_gauge, the
    ##   subsidence or uplift to produce B[0].  That is, B0+dzi_gauge=B[0].
    ##
    ##   If the assumption about B being constant is false, it will show up on
    ##   the plots!
    ##
    ##   Well, for ascee gauges, they are not all within the finest resolution
    ##   region, so be careful.

    lent=len(t)
    B  = eta - h               #Note: B should be B[0] for all time, see above  
                               #      if the gauge is in the finest level.
                               #      and turned on at all times.

    ########## IMPORTANT NOTE #########
    #
    # With GeoClaw 5.5.0, the first line in the gauge file was after 5 sec, so
    # the statements below and the note above about B being B[0] for all time
    # would be correct. 
    #Bpost_quake_finest = B[0]
    #B0 = Bpost_quake_finest - dzi_gauge
    #
    # But with GeoClaw 5.8.0, the first line in the gauge file was at 0 sec, as
    # desired from the time the gauges were to be turned on, so we must change
    # those two statements above to be:
    #
    B0 = -h[0]     #the first h in the file was really h0, pre-quake, all ASCE
                   #gauges are in the water
    Bpost_quake_finest = B0 + dzi_gauge
    #
    # and this Bpost_quake_finest should match B[-1].  Actually, if it does that
    # means the interpolation of the dtopo matches pretty well what was used to
    # compute B.

    B0_long = B0*ones(lent)

    if 0:
        print (' Initial bathymetry at the gauge before quake was: ',B0)
        print (' Final bathymetry at the gauge right after quake B0+dzi_gauge was: ',Bpost_quake_finest)
        print (' Final bathymetry at the gauge at end of job run was: ',B[-1])
        print (' The two final bathys above should be the same')
        print (' ')
    ##
    hmax=h.max()
    etamax=eta.max()
    etamax_pquake = hmax + B[-1]
    
    geo_amp[gaugeno-1]=etamax

    figure(400)
    clf()
    if 0:
    #if gaugeno in gauges_to_plot: #Make this 1 if we want all plots.  Could be 551 of them.
        subplot(211)
        plot(t, B0_long,'y')
        plot(t, B,'g')
        plot([t[0],t[0]],[B0,B[-1]],'g--')
        xlabel('')
        ylabel('Initial (Y) and Varying Bathymetry (G) (m)')
        title('Gauge %i long=%s,lat=%s,dzi=%5.2f' %(gaugeno,xlong,ylat,dzi_gauge))

        ## eta is always the height of water above the fixed datum called MHW.
        subplot(212)
        plot(t, eta, 'b')
        #plot(t, B,'g')
        xlabel('')
        ylabel('Eta (m): h+B')
        title('')

        tight_layout()

        fname_gauge = 'ASCE_Gauge%s.png' % str(gaugeno).zfill(5)
        fname = source_name + '_' + fname_gauge
        fwhere = other_figures_dir + '/' + fname
        savefig(fwhere, bbox_inches='tight')
        if 0:
            print ('Created %s' % fname )

for jkm in range(num_structures):
    print (' ')
    print ('RESULTS for +- 32.2km of structure site: ', structure_name[jkm])
    print ('                     ASCE                                                                            ')
    print ('gauge no     latitude     longitude    GeoClaw eta     ASCE eta     .8*(ASCE eta)    .8*ASCE/Geo  ')
    starting = gstart32[jkm]; ending = gend32[jkm];
    num_gauges=(ending+1.0 - starting)

    ## Index is always one less than the gaugeno into the long, lat, and amp arrays
    ascee_long_save=lon[starting-1:ending]
    ascee_lat_save=lat[starting-1:ending]
    ascee_save=amp[starting-1:ending]
    ascee_pt8_save = ascee_save*.8
    mean_ascee = mean(ascee_save)
    geoc_save=geo_amp[starting-1:ending]; 
    mean_geoc = mean(geoc_save)
    mean_ratio = mean_ascee/mean_geoc;
    
    sum_geoc=0.0; 
    for gaugeno in range(starting,ending+1):  
        ## Now print the gauge no, geoclaw eta, asce eta
        geoc=geo_amp[gaugeno-1]
        ascee=amp[gaugeno-1]
        ascee_lat=lat[gaugeno-1]
        ascee_long=lon[gaugeno-1]
        ascee_pt8= 0.8*ascee
        ratio3 = .8*ascee/geoc
        print ('%8i %12.6f %13.6f %11.3f %14.3f %13.3f %15.3f' %(gaugeno,ascee_lat,ascee_long,\
                                                       geoc,ascee,ascee_pt8,ratio3))
    print (' ')
    print ('   The number of gauges from ',starting,' to ',ending,' inclusive was: ',num_gauges)
    print ('   ASCE Mean was: %.3f ' %mean_ascee)
    print ('   GeoClaw Mean was: %.3f' %mean_geoc)
    print ('   Ratio of ASCE Mean to GeoClaw Mean was: %.4f ' %mean_ratio)

    mean_ascee_save=mean_ascee*ones(ending+1-starting)
    mean_geoc_save=mean_geoc*ones(ending+1-starting)
    error_squared = (ascee_save - geoc_save)**2
    rootmeansquare_error = sqrt(mean(error_squared))
    print ('   The rootmeansquare_error for deviation of ASCE to GEO values was: %.3f' %rootmeansquare_error)

    figure(200,figsize=(8,12))
    clf()
    if 1:
        geoc_save1=geoc_save[0:21]; geoc_save2=geoc_save[21:];
        ascee_save1=ascee_save[0:21]; ascee_save2=ascee_save[21:];
        ascee_lat_save1=ascee_lat_save[0:21]; ascee_lat_save2=ascee_lat_save[21:];
        ascee_pt8_save1=ascee_pt8_save[0:21]
        ascee_pt8_save2=ascee_pt8_save[21:];
        mean_geoc_save1=mean_geoc_save[0:21]
        mean_geoc_save2=mean_geoc_save[21:];
        mean_ascee_save1=mean_ascee_save[0:21]
        mean_ascee_save2=mean_ascee_save[21:];
        ax = subplot(211)
        plot(ascee_lat_save1, geoc_save1, 'go', label='GeoClaw eta')
        plot(ascee_lat_save1, ascee_save1, 'ro', label='ASCE eta')
        plot(ascee_lat_save1, ascee_pt8_save1, 'ro', fillstyle='none',\
             label='0.8*(ASCE eta)')
        plot(ascee_lat_save1, mean_geoc_save1, 'g-',\
             label='GeoClaw mean=%7.3f' % mean_geoc )
        plot(ascee_lat_save1, mean_ascee_save1, 'r-',\
             label='ASCE mean=%7.3f' % mean_ascee )
        legend(loc='lower center', framealpha=1)
        ylabel('Eta (m) - distance above MHW')
        xlim(xlim1[jkm],xlim2[jkm])
        ax.set_xticks([xticks1[jkm],xticks2[jkm],xticks3[jkm]])  # labeled ticks
        ax.set_xticks(ascee_lat_save1, minor=True)
        grid(True, axis='x', which='both')
        ylim(0,12)
        title('Source: %s' % source_name )

        ax = subplot(212)
        plot(ascee_lat_save2, geoc_save2, 'go', label='GeoClaw eta')
        plot(ascee_lat_save2, ascee_save2, 'ro', label='ASCE eta')
        plot(ascee_lat_save2, ascee_pt8_save2, 'ro', fillstyle='none',\
             label='0.8*(ASCE eta)')
        plot(ascee_lat_save2, mean_geoc_save2, 'g-',\
             label='GeoClaw mean=%7.3f' % mean_geoc )
        plot(ascee_lat_save2, mean_ascee_save2, 'r-',\
             label='ASCE mean=%7.3f' % mean_ascee )
        xlabel('Gauge Latitude')
        ylabel('Eta (m) - distance above MHW')
        xlim(xlim3[jkm],xlim4[jkm])
        ylim(0,12)
        ax.set_xticks([xticks4[jkm],xticks5[jkm],xticks6[jkm]])  # labeled ticks
        ax.set_xticks(ascee_lat_save2, minor=True)
        grid(True, axis='x', which='both')
        title('')
        tight_layout()

        fwhere = source_name + '_' + structure_name[jkm] + '_Compatible_32.png'
        savefig(fwhere)
        if 0:
            print('Created %s' % fwhere)
        print (' ---------- END of output for STRUCTURE: ',structure_name[jkm])
        print (' ')

seg32_geo_amp_mean=zeros(gend-gstart+1)
satisfied_pt8=zeros(gend-gstart+1)
if 1:  #Print report for all the gauges from gstart to gend. Thinking of each gauge  as a structure, try to find the
       #two gauges (one 32.2km north of the gauge and one 32.2km south of the gauge, and print asce mean and geoclaw mean
    print (' ')
    print ('RESULTS for all gauges from %s to %s : ' %(gstart,gend) )
    print (' ')
    print ('gauge_no,   latitude,   longitude,   GeoClaw_eta,    ASCE_eta,     0.8*(ASCE_eta),  .8*ASCE/Geo,  g_N,  g_S,   del_km,  ASCE_mean,  GeoClaw_mean,  Compliant')
    for gaugeno in range(gstart,gend+1):  
        ## Now print 
        key = str(gaugeno)
        geoc=geo_amp[gaugeno-1]
        ascee=amp[gaugeno-1]
        ascee_lat=lat[gaugeno-1]
        ascee_long=lon[gaugeno-1]
        ascee_pt8= 0.8*ascee
        ratio3 = .8*ascee/geoc

        ## The north and south gauges, and asce mean for hopefully +-32km of gaugno:
        g_32N = ID_32N[gaugeno-1]
        g_32S = ID_32S[gaugeno-1]
        del_lat = seg32_del_lat[gaugeno-1]
        del_km = seg32_km[gaugeno-1]
        del_mi = seg32_mi[gaugeno-1]
        mean_asce = seg32_amp_mean[gaugeno-1]

        no_seg_gauges = int(g_32S)-int(g_32N)+1
        seg_test = where(geo_amp[int(g_32N-1):int(g_32S-1+1)] >= .8*amp[int(g_32N-1):int(g_32S-1+1)],1,0)

        if (sum(seg_test) == no_seg_gauges):
            satisfied_pt8[gaugeno-1] = 1.0

        ## Now calculate the GeoClaw mean for this segment and save in seg32_geo_amp_mean
        mean_geoc = mean(geo_amp[int(g_32N-1):int(g_32S-1+1)])
        seg32_geo_amp_mean[gaugeno-1]=mean_geoc
        
        #set the compliance for this gauge
        if ((satisfied_pt8[gaugeno-1] >0.0) & (mean_geoc > mean_asce)):
            gauge_is_compliant = True
        else:
            gauge_is_compliant = False 

        print ('%8i, %11.6f, %11.6f, %9.3f, %13.3f, %13.3f, %15.3f, %6i, %5i, %8.3f, %8.3f, %10.3f,        %s' %(gaugeno,ascee_lat,ascee_long,\
                geoc,ascee,ascee_pt8,ratio3,g_32N,g_32S,del_km,mean_asce,mean_geoc,gauge_is_compliant))

    print(' ')
    compliance_condition = (satisfied_pt8 > 0.0) & (seg32_geo_amp_mean >= seg32_amp_mean)
    compliant_indices = where(compliance_condition)[0]
    compliant_gaugenos = ID[compliant_indices] 
    compliant_lon = lon[compliant_indices]
    compliant_lat = lat[compliant_indices]
    num_compliant_gaugenos = len(compliant_gaugenos)
    print(' ')
    print(' The compliant_indices were: ',compliant_indices)
    print(' The number of compliant gaugenos was: ',num_compliant_gaugenos)
    print(' ')
    print('   2028 COMPLIANT LOCATIONS ')
    print(' gauge no     longitude      latitude')
    for i in range(num_compliant_gaugenos):
        print('%8.1f %14.3f %12.3f'  %(compliant_gaugenos[i],compliant_lon[i],compliant_lat[i]) )
