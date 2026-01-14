 
"""
The function report plots special gauge output from GeoClaw and creates 
gauges_report_timesarrival.csv and gauges_report_timesarrival.txt.  For
each gauge, we plot the time series of h and h+B including the pre-quake
(B0) and post-quake (B) bathymetry in the part1 plots, and then we plot
the time series of speed (s), momentum (hs), and momentum flux (hss) in
the part2 plots.  A summary report is given as the .csv file for all
the gauges.  A more detailed .txt file is printed that also includes
the latitude and longitude that GeoClaw used for each gauge.

The function report that takes the following six arguments as below 

plot_gaugereport.report(outdir, plotdir, location, event, dtopofile, run_name)

"""

# for making plots on remote machine:
import matplotlib
matplotlib.use('Agg')  # Use an image backend

from pylab import *
import os, sys
import copy
import clawpack.pyclaw.gauges as gauges
from clawpack.geoclaw import dtopotools

root_dir = os.environ['CHT']
rundir   = os.getcwd()
print ('root_dir printed by plot_gaugesreport was: ',root_dir)
print ('rundir printed by plot_gaugesreport was: ',rundir)
print (' ')

sea_level = 0.0    #MHW runs
import interptools

def report(outdir,plotdir,location,event,dtopofile,run_name):
    other_figures_dir  = plotdir + '/_other_figures'
    if not os.path.isdir(other_figures_dir):
        os.mkdir(other_figures_dir)

    ## Print the directories used
    print('RUN_NAME was: ',run_name)
    print('   location was: ',location)
    print('   event was: ',event)
    print('   dtopofile was: ',dtopofile) 
    print('   outdir is: ',outdir)
    print('   plotdir is: ',plotdir)
    print('   other_figures_dir is: ',other_figures_dir)

    fprint_path = os.path.join(other_figures_dir,'gauges_report_timesarrival.txt')
    print('   Will send detailed .txt report file to \n   ',fprint_path)
    fprint_file = open(fprint_path,'w')
    fprint_csv_path = os.path.join(other_figures_dir,'gauges_report_timesarrival.csv')
    print('   Will send .csv report file to \n   ',fprint_csv_path)
    fprint_csv_file = open(fprint_csv_path,'w')

    def fprint(*args):
        # build a string out of all the arguments passed to fprint:
        line = ''
        for arg in args:
            line = line + str(arg)
        fprint_file.write('%s\n' % line)

    def fprint_csv(*args):
        # build a string out of all the arguments passed to fprint_csv:
        line = ''
        for arg in args:
            line = line + str(arg)
        fprint_csv_file.write('%s\n' % line)

    fprint('\nGAUGE REPORT for RUN_NAME: %s \n' %run_name)
    fprint(' ')

    fprint('sea_level in meters above MHW was: %.3f ' %sea_level)
    fprint(' ')

    figure(400, figsize=(8,8))
    figure(500, figsize=(8,8))

    #####   Pick the gauge numbers for the report. Doing maximums also for
    #####   a smaller set of gauge numbers, so pick that smaller set also.

    ### 15 Aberdeen gauges will be the big list
    gaugeno_list=[11,33,43,56,68,137,138,354,367,368,390,509,558,561,562]

    hmax_orig_dry=0.0; hmax_orig_wet=0.0;
    hmax_area=0.0; zetamax_area = 0.0; etamax_area=0.0; etamax_area_pquake=0.0;
    speedmax_area=0.0; momentummax_area=0.0; mfluxmax_area=0.0;

    hmax_orig_dry_small=0.0; hmax_orig_wet_small=0.0;
    hmax_area_small=0.0; zetamax_area_small = 0.0; etamax_area_small=0.0;
    etamax_area_small_pquake=0.0;
    speedmax_area_small=0.0; momentummax_area_small=0.0; mfluxmax_area_small=0.0;

    # Read in topography and find the average bottom deformation
    dtopo = dtopotools.DTopography()
    dtopo.read(dtopofile, 3)

    gaugeno_dict={}
    for gaugeno in gaugeno_list:
        fprint(' ')
        fprint('Gauge no: ',gaugeno)

        gauge = gauges.GaugeSolution(gaugeno, outdir)
        xlong,ylat = gauge.location

        ## Determine the subsidence or uplift at this xlong, ylat
        try:
            dzi_gauge = interptools.interp(xlong,ylat,dtopo.X,dtopo.Y,dtopo.dZ[-1,:,:])
            fprint( '(long,lat) = (',xlong,',',ylat,')' )
            fprint('Subsidence/uplift at the gauge was:  %.3f meters'  % dzi_gauge )
        except:
            dzi_gauge = 0.
            fprint( '(long,lat) = (',xlong,',',ylat,')' )
            fprint('The gauge does not overlap dtopo, so no subsidence/uplift' )
        fprint(' ')

        t = gauge.t / 60.   # convert to minutes
        q = gauge.q
        h = q[0,:]

        if 1: ##Gets u and v set properly, and does not change h
            hcopy = copy.copy(h)
            hcopy = where(hcopy>0.01, hcopy, 1.e6)
            u = q[1,:] / hcopy
            v = q[2,:] / hcopy
    
        s = sqrt(u**2 + v**2)
        momentum = h*s
        mflux = h*s*s
        eta = q[3,:]

        lent=len(t)
        B  = eta - h     

        #Assuming the first line in the gauges file gives B0
        #That is, the first line reflects the pre-event before any
        #subsidence has happened.  The water has been filled already
        #up to sea-level.

        B0 = eta[0]-h[0]
        Bpost_quake_finest = B0 + dzi_gauge

        #
        # This Bpost_quake_finest should match B[-1].  Actually, if it does that
        # means the interpolation of the dtopo matches pretty well what was used to
        # compute B, and that the gauge was in a fixed region.

        B0_long = B0*ones(lent)

        if 1:
            fprint(' Initial bathymetry at the gauge before quake in meters was: %.3f ' %B0 )
            fprint(' Final bathymetry at the gauge right after quake in meters  was: %.3f ' %Bpost_quake_finest )
            fprint(' Final bathymetry at the gauge at end of job run in meters was: %.3f' %B[-1])
            fprint(' ')

        ## Now, we can calculate zeta
        ## If the initial bathy was in the water, zeta will be B0+h
        ## but if the inital bathy was on land, we want h
        ##
        if (B0 < 0.0):  #zeta = h+B0
           zeta = B0 + h
        else:
           zeta = h
        zetamax=zeta.max()

        hmax=h.max()
        arg_hmax=h.argmax()
        tmax = t[arg_hmax]
        #First wave could be a depression wave, so check abs value for arrival
        hdeep_index = where(abs((h-h[0])) > 0.05)[0]
        if len(hdeep_index) >= 1:
            arr_index = hdeep_index[0]
            tfirst = t[arr_index]
        else:
            tfirst = -9999.

        #Note, these two below are the same when the eta[0]<=etamax.
        #Since the first line in the file gives eta[0] pre-quake. It
        #is possible that the ground drops and floods quite a bit, but
        #not up to the original elevation level eta[0].  eta.max() is
        #includes eta[0] and etamax_pquake does not.
        etamax=eta.max()
        etamax_pquake = hmax + B[-1]

        speedmax=s.max()
        momentummax=momentum.max()
        mfluxmax=mflux.max()

        ## Save the info for this gauge in the dictionary below for later printing
        value_dict={'B0': B0, 'B': B[-1], 'max_h': hmax, 'max_zeta': zetamax, 'max_eta': etamax,\
            'max_eta_pquake': etamax_pquake, 'max_speed': speedmax, 'max_momentum': momentummax,\
            'max_mflux': mfluxmax, 'dzi': dzi_gauge, 'tmax': tmax, 'tfirst': tfirst}
        gaugeno_dict[gaugeno]=value_dict

        if 1: #Change to 0 if already have the figures
            figure(400)
            clf()
            subplot(211)
            plot(t, h, 'b')
            if (B[-1] < 0): 
                Blong = B[-1]*ones(len(t))
                plot(t, -Blong+sea_level,'b--')
            xlabel('')
            ylabel('Flow depth (meters)')
            title('Gauge %i long=%s, lat=%s, B0(m.)=%6.2f, B(m.)=%6.2f \n dzi(m.)=%5.2f, max h(m.)=%5.2f,\
                 max eta post-quake(m.)=%5.2f' %(gaugeno,xlong,ylat,B0,B[-1],dzi_gauge,hmax,etamax_pquake))

            ## eta is always the height of water above the fixed datum called MHW.
            subplot(212)
            plot(t, eta, 'b')
            plot(t, B0_long,'y')
            plot(t, B,'g')
            plot([t[0],t[0]],[B0,B[-1]],'g--')
            plot(t, 0*t+sea_level,'k')
            xlabel('')
            ylabel('Eta (meters): h+B, B0(Y), B(G) (meters)')
            title('')

            tight_layout()

            fdir  = other_figures_dir
            fname_gauge = 'Gauge%s.png' %str(gaugeno).zfill(5)
            fname = run_name + '_' + 'part1_' + fname_gauge
            fwhere = fdir + '/' + fname
            savefig(fwhere, bbox_inches='tight')
            fprint('Created %s' % fname )

        if 1: #Choose 0 if already have the plots
            figure(500)
            clf()
            subplot(311)
            plot(t, s, 'b')
            xlabel('')
            ylabel('speed (m/s)')
            title('Gauge %i long=%s, lat=%s \n s_max(m/sec)=%5.2f, hs_max(m^2/sec)=%5.2f, hss_max(m^3/sec^2)=%5.2f'\
                    %(gaugeno,xlong,ylat,speedmax,momentummax,mfluxmax))
            subplot(312)
            plot(t, momentum, 'b')
            ylabel('momentum (m^2 / s)')
            title('')

            subplot(313)
            plot(t, mflux, 'b')
            xlabel('time (Minutes after earthquake)')
            ylabel('momentum flux (m^3 / s^2)')
            title('')

            tight_layout()
            fdir  = other_figures_dir
            fname_gauge = 'Gauge%s.png' %str(gaugeno).zfill(5)
            fname = run_name + '_'+ 'part2_' + fname_gauge
            fwhere = fdir + '/'+ fname
            savefig(fwhere, bbox_inches='tight')
            fprint('Created %s' % fname )
    
        ## Check to see if the maximums were changed by this gauge
        if ((hmax > hmax_orig_dry) & (B0 > 0)):
            hmax_orig_dry = hmax
        if ((hmax > hmax_orig_wet) & (B0 <= 0)):
            hmax_orig_wet = hmax
        if (hmax > hmax_area):
            hmax_area = hmax
        if (zetamax > zetamax_area):
            zetamax_area = zetamax
        if (etamax > etamax_area):
            etamax_area = etamax
        if (etamax_pquake > etamax_area_pquake):
            etamax_area_pquake = etamax_pquake
        if (speedmax > speedmax_area):
            speedmax_area = speedmax
        if (momentummax > momentummax_area):
            momentummax_area = momentummax
        if (mfluxmax > mfluxmax_area):
            mfluxmax_area = mfluxmax

    # Print the maximums over the area encompassed by all the gauges used
    fprint(' ')
    fprint('   MAXIMUM VALUES OVER ALL GAUGES ' )
    fprint(' ' )
    fprint(' Maximum value of h (flow depth) in m. over all gauges: %.3f ' %hmax_area )
    fprint(' Maximum value of h (flow depth) in m. over all originally dry gauges: %.3f ' \
                                                  %hmax_orig_dry )
    fprint(' Maximum value of h (flow depth) in m. over all originally wet gauges: %.3f ' \
                                                  %hmax_orig_wet )
    fprint(' Maximum value of zeta in m., h (originally land) or h+B0 (originally wet): %.3f '\
             %zetamax_area )
    fprint(' Maximum value of eta (h+B) in m., height above MHW: %.3f' %etamax_area )
    fprint(' Maximum value of eta post quake in m., height above MHW: %.3f' %etamax_area_pquake )
    fprint(' Maximum speed in m/sec: %.3f' %speedmax_area )
    fprint(' Maximum momentum in m^2/sec: %.3f' %momentummax_area )
    fprint(' Maximum momentum flux in m^3/sec^2: %.3f' %mfluxmax_area )
    fprint(' ')

    ####  Write the comma separated file
    fprint_csv('%5s, , , , , , , , , , , , ' %event)
    fprint_csv('       ,      ,         ,       ,   max,     max,   max,    max IE,   max,    max,       max,          ,          ' )
    fprint_csv('  Gauge,    B0,        B,    dzi,    h,     zeta,   eta,   post-eta,   s,     hs,        hss,      tmax,    tfirst' )
    fprint_csv('   No,      (m),      (m),   (m),   (m),     (m),   (m),     (m),    (m/s),  (m*m/s),  (m^3/s^2),  (min),    (min) ')  
    for key in gaugeno_dict:
        value_dict = gaugeno_dict[key]
        fprint_csv('%5i, %8.3f, %8.3f, %5.2f, %6.2f, %6.2f, %6.2f,  %6.2f,  %6.2f,  %6.2f,  %8.2f, %8.1f, %8.1f' %(key,value_dict['B0'],\
           value_dict['B'],value_dict['dzi'],value_dict['max_h'],value_dict['max_zeta'],\
           value_dict['max_eta'],value_dict['max_eta_pquake'],value_dict['max_speed'],\
           value_dict['max_momentum'],value_dict['max_mflux'],value_dict['tmax'],value_dict['tfirst']) )
    fprint_csv_file.close()
    ###  End of comma separated file

    fprint('               SUMMARY FOR EACH OF THESE GAUGES                ' )
    fprint('%5s, , , , , , , , , , , , ' %event)
    fprint('       ,      ,         ,       ,   max,     max,   max,    max IE,   max,    max,       max,          ,          ' )
    fprint('  Gauge,    B0,        B,    dzi,    h,     zeta,   eta,   post-eta,   s,     hs,        hss,      tmax,    tfirst' )
    fprint('   No,      (m),      (m),   (m),   (m),     (m),   (m),     (m),    (m/s),  (m*m/s),  (m^3/s^2),  (min),    (min) ')  
    for key in gaugeno_dict:
        value_dict = gaugeno_dict[key]
        fprint('%5i, %8.3f, %8.3f, %5.2f, %6.2f, %6.2f, %6.2f,  %6.2f,  %6.2f,  %6.2f,  %8.2f, %8.1f, %8.1f' %(key,value_dict['B0'],\
           value_dict['B'],value_dict['dzi'],value_dict['max_h'],value_dict['max_zeta'],\
           value_dict['max_eta'],value_dict['max_eta_pquake'],value_dict['max_speed'],\
           value_dict['max_momentum'],value_dict['max_mflux'],value_dict['tmax'],value_dict['tfirst']) )
    fprint_file.close()
