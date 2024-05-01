from pylab import *
import obspy


distances = open('GNSS_distances_to_synthetic_stations').readlines()

print('Omitting stations more than 2 km from a grid point:')
distance = {}
stations = []
for line in distances[1:]:
    tokens = line.split()
    station = tokens[0]
    dist = float(tokens[1])
    if dist>2000:
        print('%s : %.1f km' % (station, dist/1e3))
    else:
        stations.append(station)
        distance[station] = dist

print('Keeping %i stations' % len(distance))


def load_gnss(event, station, make_plot=True):
    print('Loading GNSS data for event %s, station %s' % (event,station))
    fname = 'GNSS_' + event + '.h5'
    waveforms = obspy.read(fname)

    wstation = waveforms.select(station=station)
    if len(wstation) == 0:
        raise ValueError('station %s not found in %s' % (station,fname))
        
    Nstation = wstation.select(channel='CXN')[0]
    Estation = wstation.select(channel='CXE')[0]
    Zstation = wstation.select(channel='CXZ')[0]
    t = Zstation.times()

    if make_plot:
        figure(1,figsize=(10,6))
        clf()
        plot(t, Nstation, 'r', label='N')
        plot(t, Estation, 'g', label='E')
        plot(t, Zstation, 'b', label='Z')
        grid(True)
        legend(framealpha=1)
        xlabel('time (s)')
        ylabel('displacement (m)')
        title('Station %s (distance to grid point = %.1f km) \nScenario %s' \
                % (station,dist/1e3,event))
                
        fname_png = '%s_%s.png' % (event,station)
        savefig(fname_png)
        print('Created ',fname_png)

    return t, Nstation, Estation, Zstation

def make_npy_all_stations(event):

    j = 0
    for station in stations:
        t, Nstation, Estation, Zstation = load_gnss(event, station,
                                                    make_plot=True)
        if j==0:
            d = vstack((t, Nstation, Estation, Zstation)).T
        else:
            d = hstack((d, vstack((Nstation, Estation, Zstation)).T))
        j = j+1
            
    print('d now has shape: ',d.shape)
    fname = 'gnss_%sstations_%s.npy' % (j, event)
    save(fname, d)
    print('Created ', fname)
    
    fname2 = 'gnss_%sstations_list_%s.txt' % (j, event)
    with open(fname2,'w') as f:
        f.write('# list of stations in file %s\n' % fname)
        f.write('# columns are t, then N,E,Z components for each station\n')
        for station in stations:
            f.write('%s\n' % station)

    
def load_sample():
    event = 'buried-locking-str10-middle'
    station = 'BAMF'

    dist = nan
    for line in distances:
        tokens = line.split()
        if tokens[0].lower() == station.lower():
            dist = float(tokens[1])
            print('Distance from station %s to grid point = %.0f meters' \
               % (station,dist))
    
    load_gnss(event, station, make_plot=True)


if __name__=='__main__':

    event = 'buried-locking-str10-middle'
    make_npy_all_stations(event)
