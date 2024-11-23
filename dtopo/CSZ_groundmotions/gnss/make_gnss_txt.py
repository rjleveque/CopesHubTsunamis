from pylab import *
import os
import read_gnss

outdir = 'gnss_txt_files/'
os.system('mkdir -p %s' % outdir)

def make_gnss_txt_file(event,station):
    #event = 'buried-locking-str10-deep'
    #station = 'seat'
    #station = 'ufda'
    event_station = '%s_GNSSstation_%s' % (event,station)

    t,N,E,Z = read_gnss.load_gnss(event, station, make_plot=True)
    fname = outdir + event_station + '.png'
    savefig(fname, bbox_inches='tight')
    print('Created ', fname)

    dd = vstack((t,N,E,Z)).T

    fname = outdir + event_station + '.txt'
    header = 't,N,E,Z for event %s, GNSS station %s' % (event,station)
    savetxt(fname, dd, header=header, fmt='%18.10f')
    print('Created ', fname)
    
    # accelerations for seiching:
    
    Nv = diff(N)/diff(t)
    Na = diff(Nv)/diff(t[1:])
    Ev = diff(E)/diff(t)
    Ea = diff(Ev)/diff(t[1:])
    ta = t[1:-1]
    accel = vstack((ta,Na,Ea)).T
    fname = outdir + event_station + '_accel.txt'
    savetxt(fname, accel, comments='', header='%i' % len(ta), fmt='%14.4f')
    print('Created ', fname)
    
    figure(22, figsize=(12,6))
    clf()
    absa = maximum(abs(Na), abs(Ea))  # maximum abs accel at each time
    maxa = max(abs(Na).max(), abs(Ea).max())
    i1 = where(absa > 0.1*maxa)[0].min()
    i2 = where(absa > 0.1*maxa)[0].max()
    tmin = floor(t[i1]/25) * 25
    tmax = ceil(t[i2]/25) * 25
    plot(ta, Ea, 'g', label='E accel')
    plot(ta, Na, 'r', linewidth=0.5, label='N accel')
    xlim(tmin,tmax)
    grid(True)
    xlabel('time (sec)')
    ylabel('acceleration (m/s^2)')
    title('Station %s acceleration\nScenario %s' % (station,event))
    legend(loc='upper right')
    fname = outdir + event_station + '_accel.png'
    savefig(fname, bbox_inches='tight')
    print('Created ', fname)
