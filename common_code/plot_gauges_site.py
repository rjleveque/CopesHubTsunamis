
"""
make gauge plots for a particular location and event.
This should now work for any location, sets gaugenos based on '../gauges_B0.csv'
(To use a subset adjust gaugenos in main program.)

make_plot makes plots for a single gaugeno

make_all_plots_and_report make plots for all gauges listed in gaugenos,
    and also makes gauge report and csv file.
    This can be called in parallel for a list of events in a multirun directory.

When run from command line (main program) for a single event, it reads sea_level
from geoclaw.data.  When looping over events, pass the correct value into
make_all_plots_and_report (default is 0.).

"""

import sys,os
if 'matplotlib' not in sys.modules:
    import matplotlib
    matplotlib.use('Agg')  # Use an image backend

from pylab import *
import clawpack.pyclaw.gauges as gauges
from clawpack.clawutil.data import ClawData

# top level directory for this project:
root_dir = os.environ['CHT']   # assuming environment variable set
sys.path.insert(0, os.path.join(root_dir, 'common_code'))
from gauge_B0_tools import read_gauge_B0




def make_plot(gaugeno, location, event, outdir, plotdir, B0, sea_level):
    gauge = gauges.GaugeSolution(gaugeno, outdir)
    x,y = gauge.location
    t = gauge.t / 60.   # convert to minutes
    q = gauge.q
    h = q[0,:]
    #h = where(h>0.01, h, 1.e6)
    u = divide(q[1,:], h, where=h>0.01, out=zeros(h.shape))
    v = divide(q[2,:], h, where=h>0.01, out=zeros(h.shape))
    s = sqrt(u**2 + v**2)
    ###  Loyce fixed mflux and added momentum, mflux=h*s was wrong
    momentum = h*s
    mflux = h*s*s
    eta = q[3,:]
    B = eta - h

    ### Do a test here of the max and min B before the loop below
    Bmin = B.min(); Bmax = B.max();

    #Using only the maximum level info in the gauge file
    #All the gauges here are in a fixed region, so shouldnt have to do this
    maxlevel = gauge.level.max()
    if gauge.level.min() < maxlevel:
        for a in [h,u,v,s,momentum,mflux,eta,B]:
            a = where(gauge.level == maxlevel, a, nan)

    #The quantities h,u,v,s,momentum, mflux,eta,B have been changed to nan
    #if they weren't on finest level for that gauge.
    #The gauge was turned on at time 0.  Don't know when the subsidence
    #was complete for the dynamic rupture.

    hmax = nanmax(h); smax = nanmax(s); mfluxmax = nanmax(mflux); etamax = nanmax(eta);
    momentummax = nanmax(momentum)

    #Say B post-quake is the B post-quake on the finest level for this gauge
    #Note, the gauge won't subside instantaneously, etamax might be achieved before subsidence here.
    #etamax_pquake is the eta corresponding to the maximum hmax plus the post-quake bathymetry B_post.
    ind_hmax = argmax(where(isnan(h), -inf, h))
    ind_B=where(isnan(B) == False)[0].max()
    B_post = B[ind_B]
    etamax_pquake = hmax + B_post

    ### Find h0 from B0 and sea_level of this particular job run
    ### Assuming no etainit applies to this particular gauge.
    ### If etainit does apply, got to know the value of etainit used,
    ### then h0 = etainit - B0

    special_gauge = False    #say we have determined this gauge doesn't get etainit
    no_force_dry = True      #didn't force any gauges to be dry
    if (special_gauge):
        etainit_level = 2.0  #example, never used
    else:
        etainit_level = sea_level
    fill_up_condition = (B0 < etainit_level) & no_force_dry
    if (fill_up_condition):         #fill up to etainit_level initially which is sealevel here.
        h0 = etainit_level - B0
    else:
        h0 = 0.0                    #location stays dry

    ######### Find the time of the maximum hmax, and the first time h0 INCREASES.
    tmax = t[ind_hmax]
    hdeep_index = where((h-h0) > 0.05)[0]
    if len(hdeep_index) >= 1:
        arr_index = hdeep_index[0]
        tfirst = t[arr_index]
    else:
        tfirst = -9999.
    #########

    if 1:
        figure(400, figsize=(8,8))
        clf()

        subplot(311)
        plot(t, h, 'b')
        xlabel('')
        ylabel('Flow depth (m)')
        grid(linewidth=0.5)
        title('Gauge %i at x = %.5f, y = %.5f In %s for Event %s \n \
          h0 = %.2f, hmax = %.2f, smax = %.2f, mfluxmax = %.2f'
          % (gaugeno,x,y,location,event,h0,hmax,smax,mfluxmax))

        subplot(312)
        plot(t, s, 'b')
        xlabel('')
        ylabel('speed (m/s)')
        grid(linewidth=0.5)
        title('')

        subplot(313)
        plot(t, mflux, 'b')
        xlabel('time (Minutes after earthquake)')
        ylabel('momentum flux (m^3 / s^2)')
        grid(linewidth=0.5)
        title('')

        fname = plotdir + '/%s_%s_Gauge%s.png' \
            % (location,event,str(gaugeno).zfill(5))
        savefig(fname, bbox_inches='tight')
        print('Created %s' % fname)

        if 0:
            # plot Bathymetry as sanity check:
            figure(500, figsize=(8,8))
            clf()
            plot(t, B, 'g')
            xlabel('')
            ylabel('Bathymetry B (m)')
            grid(linewidth=0.5)
            title('Gauge %i at x = %.5f, y = %.5f In %s for Event %s \n \
              B_post = %.2f, B0 = %.2f' 
              % (gaugeno,x,y,location,event,B_post,B0))
            fname2 = plotdir + '/%s_%s_B_at_Gauge%s.png' \
                    % (location,event,str(gaugeno).zfill(5))
            savefig(fname2, bbox_inches='tight')
            print('Created %s' % fname2)

        return B_post,h0,hmax,smax,momentummax,mfluxmax,etamax,etamax_pquake,tmax,tfirst

def make_all_plots_and_report(outdir, plotdir, location, event, 
                              gaugenos='all', sea_level=0.):
    gaugeno_dict = {}


    # read dictionary of B0 values indexed by gaugeno:
    fname_gauge_B0 = root_dir + '/geoclaw_runs/sites/%s/' % location \
                     + 'gauges_B0.csv'
    print('Reading B0 at gauges from: ',fname_gauge_B0)
    gauge_B0 = read_gauge_B0(fname_gauge_B0)
    all_gaugenos = gauge_B0.keys()
    print('At %s, found all_gaugenos = %s' % (location,all_gaugenos))
    if gaugenos == 'all':
        gaugenos = all_gaugenos
    else:
        for gaugeno in gaugenos:
            if gaugeno not in all_gaugenos:
                print('*** warning, specified gaugno %i not found for %s' \
                     % (gaugeno, location))

    for gaugeno in gaugenos:
        #OLD
        #B0,B_post,hmax,smax,momentummax,mfluxmax,etamax,etamax_pquake,tmax,tfirst =\
        #      make_plot(gaugeno, location, event, outdir, plotdir)

        # NEW: Passing B0 and sea_level in
        B0 = gauge_B0[gaugeno]
        B_post,h0,hmax,smax,momentummax,mfluxmax,etamax,etamax_pquake,tmax,tfirst =\
              make_plot(gaugeno, location, event, outdir, plotdir, B0, sea_level)

        ## Save the info for this gauge in the dictionary below for later printing
        value_dict={'B0': B0, 'B': B_post, 'h0':h0, 'max_h': hmax, 'max_eta': etamax,\
            'max_eta_pquake': etamax_pquake, 'max_speed': smax, 'max_momentum':momentummax, \
            'max_mflux': mfluxmax, 'tmax': tmax, 'tfirst': tfirst}
        gaugeno_dict[gaugeno]=value_dict

    ## print gauge report as a .txt file
    fprint_path = os.path.join(plotdir, 'gauges_report_timesarrival.txt')
    print('Will send output to \n   ',fprint_path)

    with open(fprint_path,'w') as fprint_file:

        def fprint(*args):
            # build a string out of all the arguments passed to fprint:
            line = ''
            for arg in args:
                line = line + str(arg)
            fprint_file.write('%s\n' % line)


        fprint(' ')
        fprint('\nGAUGE REPORT\n')
        fprint('EVENT: %s\n' %event)
        fprint('LOCATION: %s\n' %location)

        fprint(' ')
        fprint('               SUMMARY FOR EACH OF THESE GAUGES                ' )
        fprint('                                  max     max    max       max    max       max                      ' )
        fprint(' Gauge     B0        B     h0      h      IE    post-IE     s     hs        hss        tmax    tfirst' )
        fprint('   No     (m)       (m)    (m)    (m)     (m)     (m)     (m/s)  (m*m/s)  (m^3/s^2)    (min)   (min) ')
        for key in gaugeno_dict:
            value_dict = gaugeno_dict[key]
            fprint('%5i %8.3f %8.3f %6.2f  %6.2f  %6.2f  %6.2f  %6.2f  %8.2f  %8.2f %8.1f %8.1f' %(key,value_dict['B0'],\
                   value_dict['B'],value_dict['h0'],value_dict['max_h'],\
                   value_dict['max_eta'],value_dict['max_eta_pquake'],value_dict['max_speed'],\
                   value_dict['max_momentum'],value_dict['max_mflux'],value_dict['tmax'],value_dict['tfirst']) )

        fprint(' ')
        fprint(' ')

    ## print gauge report as a .csv file
    fname = os.path.join(plotdir, 'gauges_report_timesarrival.csv')
    with open(fname,'w') as f:
        f.write('# \n')
        f.write('# GAUGE REPORT\n')
        f.write('# EVENT: %s\n' %event)
        f.write('# LOCATION: %s\n' %location)
        f.write('# \n ')
        f.write('      ,  , ,  ,max,max,max,max,max,max,,\n')
        f.write('Gauge ,B0,B,h0,h,IE,post-IE,s,hs,hss,tmax,tfirst\n')
        f.write('Number,(m),(m),(m),(m),(m),(m),(m/s),(m*m/s),(m^3/s^2),(min),(min)\n')

        for key in gaugeno_dict:
            value_dict = gaugeno_dict[key]
            f.write('%5i, %8.3f, %8.3f, %6.2f, %6.2f,  %6.2f,  %6.2f,  %6.2f,  %8.2f,  %8.2f, %8.1f, %8.1f\n'\
                   %(key,value_dict['B0'],value_dict['B'],value_dict['h0'],value_dict['max_h'],\
                   value_dict['max_eta'],value_dict['max_eta_pquake'],value_dict['max_speed'],\
                   value_dict['max_momentum'],value_dict['max_mflux'],value_dict['tmax'],value_dict['tfirst']) )

    print(' Created ',fname)
    close('all')

if __name__ == '__main__':

    sys.path.insert(0,'.')
    from params import event, location

    outdir = os.path.abspath('./_output')
    plotdir = os.path.abspath('./_plots')
    os.system('mkdir -p %s' % plotdir)
    print('Will take output from \n    %s\n   and send plots to \n    %s' \
            % (outdir,plotdir))

    geodata = ClawData()
    fname = outdir + '/geoclaw.data'
    geodata.read(fname,force=True)
    sea_level = geodata.sea_level
    print('+++ sea_level = %.3f' % sea_level)

    gaugenos = 'all'

    if 0:
        # test on small set:
        gaugenos = range(1001,1003)

    make_all_plots_and_report(outdir, plotdir, location=location, event=event,
                              gaugenos=gaugenos, sea_level=sea_level) 

    
