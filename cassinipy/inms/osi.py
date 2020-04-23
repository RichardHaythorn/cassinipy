import csv
import datetime
import glob
import numpy as np


def read_INMS_1A(filename):
    '''
    Reads INMS 1A data 
    '''
    INMSdata = {'datetime':[],'alt_t':[],'amu/q':[],'c1':[]}
   
    datacounter = 0 
    
    tempmass = 0
    tempv = 0
    
    with open('data/titan/inms/'+filename+ ".csv",'r') as csvfile:
        tempreader = csv.reader(csvfile, delimiter=',')
        next(tempreader)
        next(tempreader)
        next(tempreader)
        for rowcounter, row in enumerate(tempreader):      
            if row[7] == "osi": 
                #print(row[0],row[26],row[74])
                
                if tempmass == float(row[26]) and tempv < float(row[35]):
                    continue
                
                INMSdata['datetime'].append(datetime.datetime.strptime(row[0],"%Y-%jT%H:%M:%S.%f"))
                INMSdata['alt_t'].append(float(row[42]))
                INMSdata['amu/q'].append(float(row[26]))
                INMSdata['c1'].append(int(row[74]))
                
                tempmass = float(row[26])
                tempv = float(row[35])
                
    return INMSdata

def INMS_massdata(tempdata,mass):
    
    INMS_mass_altitudes = []
    INMS_mass_counts = []
    
    masscounter = 0 
    for counter, i in enumerate(tempdata['amu/q']):
        
        if i == mass:
            #print(counter,i,tempdata['datetime'][counter],tempdata['c1'][counter])
            INMS_mass_altitudes.append(tempdata['alt_t'][counter])
            INMS_mass_counts.append(tempdata['c1'][counter])
            
    
    return INMS_mass_altitudes,INMS_mass_counts
    
def plot_INMS_massdata(flyby,masses,INMSmasstrendax=None):
    
    colordict = {"t55":"C0","t56":"C1","t57":"C2","t58":"C3","t59":"C4"}
    
    tempdata = INMSdatadict[flyby]
    
    INMSmassfig, INMSmassax = plt.subplots()
    INMSmassax.minorticks_on()    
    INMSmassax.grid(b=True, axis='both', which='major', color='k', linestyle='-',alpha=0.5)
    INMSmassax.grid(b=True, axis='both', which='minor', color='k', linestyle='--',alpha=0.25)
    INMSmassax.set_xlabel("Alt")
    INMSmassax.set_ylabel("Counts (high sensitivity)")
    
    if INMSmasstrendax != None:
        temptrendlist = []
    
    for counter,i in enumerate(masses):
        x, y = INMS_massdata(tempdata,i)
        lowalty = np.array(y)[np.array(x) < 1000]
        lowaltx = np.array(x)[np.array(x) < 1000]
        INMSmassax.scatter(lowaltx,lowalty,label=i,color='C'+str(counter))
        
        z = np.polyfit(lowaltx, lowalty, 1)
        p = np.poly1d(z)
        INMSmassax.plot(lowaltx,p(lowaltx),color='C'+str(counter),linestyle='--')
        print(flyby,i,p.c)
        temptrendlist.append(p.c[0]) 
        
    if INMSmasstrendax != None:
        INMSmasstrendax.plot(masses,temptrendlist,color=colordict[flyby],marker='o',label=flyby)

    INMSmassax.set_xlim((950,1000))    
    INMSmassax.set_title("INMS data " + tempdata['datetime'][0].isoformat())
    INMSmassfig.legend()