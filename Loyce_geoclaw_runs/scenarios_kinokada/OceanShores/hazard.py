
from pylab import *
import os
import matplotlib
matplotlib.use('Agg')  # Use an image backend

#########  FIX for TACC
#The program should be run from location_dir on the command line: python hazard.py > output_hazard.txt
#On tacc, need to get location_dir set to the directory where geoclaw_plots is (on the scratch directory).
#root_dir = os.environ['CHT']
#location_dir = root_dir +\
#     '/Loyce_geoclaw_runs/scenarios_kinokada/HohRiver'

location_dir = os.getcwd()
########

#########  Set these parameters
### Set instantaneous or Kinokada displacement
#instant = True  
instant = False

### Set the location for printing purposes
location = 'OceanShores'

### Set the names of the earthquake events for analysis
event_stem1=['BL10D','BL10M','BL10S','BL13D',\
   'BL13M','BL13S','BL16D', 'BL16M', 'BL16S']
event_stem2=['BR10D','BR10M','BR10S','BR13D',\
   'BR13M','BR13S','BR16D', 'BR16M', 'BR16S']
event_stem3=['FL10D','FL10M','FL10S','FL13D',\
   'FL13M','FL13S','FL16D', 'FL16M', 'FL16S']
event_stem4=['FR10D','FR10M','FR10S','FR13D',\
   'FR13M','FR13S','FR16D', 'FR16M', 'FR16S']
event_stem = event_stem1 + event_stem2 + event_stem3 + event_stem4

### Set the column in the gauges_report_timesarrival.csv file
### where the qoi (quantity of interest) is for the probability analysis
### Columns are numbered 0, 1, 2, etc.  Column 5 is hmax, for example. 
### Column 6 is hmax-h0
qoi_column = 6
qoi_title = 'Exceedance levels of hmax-h0 in meters'
##########


########### Should not have to change below this line,

location_hazard_plots = location_dir +'/hazard_plots'
location_hazard_csv = location_dir +'/hazard_csv'
if not os.path.isdir(location_hazard_plots):
        os.mkdir(location_hazard_plots)
if not os.path.isdir(location_hazard_csv):
        os.mkdir(location_hazard_csv)

#Read the first report file to find out how many gauges there are.
#There will be the same number of gauges in all report files.
if instant:
    report_dir = location_dir + '/geoclaw_plots/_plots_' + event_stem[0] +\
                 '_instant/_other_figures'
    event_name = event_stem[0]+'_instant'
else:
    report_dir = location_dir + '/geoclaw_plots/_plots_' + event_stem[0] +\
                 '/_other_figures'
    event_name = event_stem[0]
report_file = report_dir + '/gauges_report_timesarrival.csv'
ID=loadtxt(report_file,delimiter=',',skiprows=4,usecols=(0),\
                unpack=True)
nogauges = len(ID)

####
#Set the no_hbar_vals and the hbar_vals, will be the same for all gauges
#Let hbar_vals go from 0 to 15 in increments of .2
no_hbar_vals = 76
hbar_vals=zeros(no_hbar_vals)
for ihbar in range(no_hbar_vals):
    hbar_vals[ihbar]=.2*ihbar

####
# Go over all the events j and store CP[j] for that event.
# Also for each event, access the report_file and save the
# quantity of interest (qoi) for all gauges for that event
# Also create the subdirectories for the event probabilities.

noevents = len(event_stem)
CP = ones(noevents)

##qoi is the quantity of interest for the hazard curve
qoi=zeros((nogauges,noevents))
for j in range(noevents):
    if instant:
        report_dir = location_dir + '/geoclaw_plots/_plots_' + event_stem[j] +\
                     '_instant/_other_figures'
        event_name = event_stem[j]+'_instant'
    else:
        report_dir = location_dir + '/geoclaw_plots/_plots_' + event_stem[j] +\
                     '/_other_figures'
        event_name = event_stem[j]
    report_file = report_dir + '/gauges_report_timesarrival.csv'

    one_third = 1.0/3.0
    CP[j]=.5*one_third
    if 'B' in event_name:
        CP[j]=CP[j]*0.75
    else:
        CP[j]=CP[j]*0.25
    if 'D' in event_name:
        CP[j]=CP[j]*0.3
    elif 'M' in event_name:
        CP[j]=CP[j]*0.5
    else:
        CP[j]=CP[j]*0.2
 
    ### HERE, store the quantity of interest, qoi from the gauge report
    ### for this event j for all gauges (one gauge per row of qoi) by
    ### choosing the correct column using usecols.  Columns start their
    ### numbering at 0 with python.  Here qoi_column=6 corresponded to hmh0
    ### in the gauges report.
    qoi[:,j]=loadtxt(report_file,delimiter=',',skiprows=4,usecols=(qoi_column),\
                    unpack=True)

    #make subdirectory under location_hazard_csv for this event's probs
    event_prob_csv = location_hazard_csv + '/'+ event_name
    if not os.path.isdir(event_prob_csv):
        os.mkdir(event_prob_csv)

    print(' ')
    print('EVENT: ',event_name,' for location: ',location)
    print(' location_dir: ',location_dir)
    print(' location_hazard_plots: /hazard_plots in location_dir')
    print(' location_hazard_csv: /hazard_csv in location_dir')
    print(' event_prob_csv: ' + '/' +event_name+' in location_hazard_csv')
    report_str1 = '/geoclaw_plots/_plots_' + event_name + '/_other_figures'
    report_str2 = ' report_dir: ' + report_str1 + ' in location_dir'
    print(report_str2)
    print(' report_file:  gauges_report_timesarrival.csv in report_dir') 
    print(' conditional probability of occurrence: ',CP[j])
    print(' ')

def prob(qoi,CP,ID,hbar_vals,event_stem,instant):
    no_hbar_vals = len(hbar_vals)
    nogauges = len(ID)
    noevents = len(event_stem)
    #Make the hazard curve probabilities for all gauges for this location.
    for igauge in range(nogauges):
        gauge_no = int(ID[igauge])
        print(' ')
        print('GAUGE NO: ',gauge_no)
        loc_gauge_hazard =zeros(no_hbar_vals)
        for j in range(noevents):
            if instant:
                report_dir = location_dir + '/geoclaw_plots/_plots_' + event_stem[j] +\
                         '_instant/_other_figures'
                event_name = event_stem[j]+'_instant'
            else:
                report_dir = location_dir + '/geoclaw_plots/_plots_' + event_stem[j] +\
                         '/_other_figures'
                event_name = event_stem[j]

            event_prob_csv =  location_hazard_csv + '/' + event_name 

            qoi_thisgauge = qoi[igauge,j]
            pbar_thisgauge = zeros(no_hbar_vals)
            for ihbar in range(no_hbar_vals):
                if (qoi_thisgauge > hbar_vals[ihbar]):
                    pbar_thisgauge[ihbar] = CP[j]
            #Now we have pbar_thisgauge for this gauge for this event. Write 
            #this as a 2-column matrix to event_prob_csv
            event_prob = zeros((no_hbar_vals,2))
            event_prob[:,0]=hbar_vals; event_prob[:,1]=pbar_thisgauge;
            which_gauge = '%s_prob_gauge_%s.csv' %(event_name,gauge_no)
            fname = event_prob_csv + '/' + which_gauge
            savetxt(fname,event_prob,fmt='%.6f',delimiter=',')
            print ('      In event_prob_csv created: ',which_gauge)

            #Update loc_gauge_hazard for this gauge from this event
            #appropriately
            loc_gauge_hazard = loc_gauge_hazard +\
            pbar_thisgauge - pbar_thisgauge*loc_gauge_hazard

        #Now, loc_gauge_hazard is completed over all events for this gauge
        #Now write the location_prob to the directory location_hazard_csv 
        #directory location_hazard_csv for each gauge. Also save a hazard
        #curve as a .png file for each gauge in the directory
        #location_hazard_plots. 

        location_gauge_csv = '%s_prob_gauge_%s.csv' %(location,gauge_no)
        fname = location_hazard_csv + '/' + location_gauge_csv 
        location_prob = zeros((no_hbar_vals,2))
        location_prob[:,0]=hbar_vals; location_prob[:,1]=loc_gauge_hazard;
        savetxt(fname,location_prob,fmt='%.6f',delimiter=',')
        print('      In location_hazard_csv created: ',location_gauge_csv)
        #
        #Could also save a .png file of this hazard curve plot
        figure(100)
        clf()
        plot(location_prob[:,0],location_prob[:,1],'b')
        xlabel(qoi_title)
        ylabel('Conditional Probability of Exceedance')
        title('Location: %s, Gauge %i: ' %(location,gauge_no))
        location_gauge_png = '%s_prob_gauge_%s.png' %(location,gauge_no)
        fname_png = location_hazard_plots + '/' + location_gauge_png
        savefig(fname_png, bbox_inches='tight')
        print('      In location_hazard_plots created: ',location_gauge_png)

prob(qoi,CP,ID,hbar_vals,event_stem,instant)
