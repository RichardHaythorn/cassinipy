import csv
import datetime
import glob
import numpy as np
from cassinipy.misc import toTimestamp

titan_flybydates = {'t17':[2006,9,7],
                    't20':[2006,10,25],'t27':[2007,3,26],
                    't40':[2008,1,5],'t46':[2008,11,3],'t47':[2008,11,19],
                    't55':[2009,5,21],'t56':[2009,6,6],'t57':[2009,6,22],'t58':[2009,7,8],'t59':[2009,7,24],
                    't83':[2012,5,22]}


def read_LP_V1(flyby):
    '''
    Reads RPWS-LP data in V1.TAB format
    '''
    LPdata = {'datetime':[],'RADIAL_DISTANCE':[],'ELECTRON_NUMBER_DENSITY':[],'ELECTRON_TEMPERATURE':[],'SPACECRAFT_POTENTIAL':[]}
   
    if flyby[0] == 't':
        moon = 'titan'
    if flyby[0] == 'e':
        moon = 'enceladus'
   
    print('data/' + moon + '/lp/RPWS_LP_T_2009*' + flyby.upper() + "_V1.TAB")
    with open(glob.glob('data/' + moon + "/lp/RPWS_LP_T_" + str(titan_flybydates[flyby][0]) + "*" + flyby.upper() + "_V1.TAB")[0],'r') as csvfile:
        tempreader = csv.reader(csvfile, delimiter=',', quotechar='|')
        for row in tempreader: 
            if abs(float(row[4])) < 1e32:               
                LPdata['datetime'].append(datetime.datetime.strptime(row[0],"%Y-%m-%dT%H:%M:%S.%fZ"))
                LPdata['RADIAL_DISTANCE'].append(float(row[1]))
                LPdata['ELECTRON_NUMBER_DENSITY'].append(float(row[2]))
                LPdata['ELECTRON_TEMPERATURE'].append(float(row[3]))
                LPdata['SPACECRAFT_POTENTIAL'].append(float(row[4]))

    return LPdata
    
def inst_RPWS_LP(LPdata,datetime):
    '''
    Returns Langmuir probe derived spacecraft potential at a single datetime value
    '''
    counter=0
    while LPdata['datetime'][counter] < datetime:
        counter +=1
        
    lptimestamps = [toTimestamp(d) for d in LPdata['datetime'][counter-1:counter+1]]
    spacecraftpotential = np.interp(toTimestamp(datetime),lptimestamps,LPdata['SPACECRAFT_POTENTIAL'][counter-1:counter+1])
    return spacecraftpotential