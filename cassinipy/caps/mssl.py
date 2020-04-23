import datetime
import numpy as np
import pandas as pd
from scipy.io import readsav
from sunpy.timeseries import TimeSeriesMetaData, GenericTimeSeries
from sunpy.util import MetaDict
from sunpy.time import TimeRange


ibscalib = readsav('calib\ibsdisplaycalib.dat')
elscalib = readsav('calib\geometricfactor.dat')
sngcalib = readsav('calib\sngdisplaycalib.dat')

def generate_timeseries_caps_mssl(data,anodefan,elsdatatype='data'):
    '''
    
    Returns a SunPy Timeseries from processed MSSL data
    '''
    #TODO add multi fielded data?
    #TODO fix metadata
    
    if elsdatatype in data.keys():
        dataframe = pd.DataFrame(np.transpose(data[elsdatatype][:,anodefan,:].byteswap().newbyteorder()))
        currentdate = datetime.datetime.strptime(data['sdate'].decode('ascii'),"%d-%b-%Y")
        for key, value in data.items():            
             if isinstance(value,np.ndarray):
                if value.ndim == 1 and value.shape[0] == data[elsdatatype].shape[2]:
                    dataframe[key]= value.byteswap().newbyteorder()
        dataframe = dataframe.set_index(pd.DatetimeIndex([currentdate + datetime.timedelta(seconds=x) for x in dataframe['secofday']]))
        tr = TimeRange(dataframe.index[0],dataframe.index[-1])
        dataframe.drop(columns=['time_ut','secofday','endsec','timehrs','hhmmss'],inplace=True)        
        for i in range(63):
            dataframe.rename(columns={i:str(elscalib['earray'][i])},inplace=True)
        
        caps_metadata = TimeSeriesMetaData(timerange=tr,colnames=['ELS']) 
        #caps_metadata.append(tr,['ELS'],MetaDict([('formatid',dataframe['formatid'])]))
        
        capstimeseries = GenericTimeSeries(dataframe,meta=caps_metadata)
        
    if 'sngdata' in data.keys():
        print("ims")
        instrument="ims"
    if 'ibsdata' in data.keys():
        print("ibs")
        instrument="ibs"
    
    return capstimeseries
        
   # data_metadata = TimeSeriesMetaData(meta=(instrument))   
    
    #return sunpy.timeseries.GenericTimeSeries(meta=data_metadata)
        
        
def CAPS_acutation(data,tempdatetime):
    '''
    Returns the CAPS actuation give a datetime.datetime
    '''
    

    for counter, i in enumerate(data['times_utc']):
        if i >= tempdatetime:
            slicenumber = counter
            break           
    if 'sngact' in data.keys():
        value = data['sngact'][slicenumber]   
    if 'actuator' in data.keys():
        value = data['actuator'][slicenumber]  
        

    return value
    
    
#TODO, make into function for ELS
#TODO, make into function for IBS

def CAPS_ELS_localramangle(tempdatetime,elsdata,anodes=False):
    act = CAPS_acutation(elsdata,tempdatetime)
    ramdir_SC = cassini_ramdirection_SCframe(tempdatetime,output=False)
    ELSvecs = rotate_CAPS_SCframe(act,'els',anodes=anodes)
    
    if anodes == True:
        angle = np.zeros((8))
        for anodenumber, temp in enumerate(ELSvecs):            
            angle[anodenumber] = np.arccos(spice.vdot(temp,ramdir_SC))*spice.dpr()
        
    if anodes == False:
        angle = np.arccos(spice.vdot(ELSvecs,ramdir_SC))*spice.dpr()
    
    return angle

def CAPS_ELS_FOVcentre_azi_elv(tempdatetime,elsdata,anodes=False):
    '''
    Either returns azi and elv for centre of anodes, or centre of full FOV
    '''
    act = CAPS_acutation(elsdata,tempdatetime)
    ELSvecs = rotate_CAPS_SCframe(act,'els',anodes=anodes)
   
   
    if anodes == True:
        ELS_azi = act
        ELS_elv = np.arange(-70,90,20)
            
        
    if anodes == False:
        ELS_azi = act
        ELS_elv = 0
        
    
    
    return ELS_azi, ELS_elv

def CAPS_IBS_FOVcentre_azi_elv(tempdatetime,elsdata):
    '''
    Returns azimuthal and elevation ram angles
    
    '''
    #TODO : Remove dependence on ELS data
    #TODO: Add all fans
     
    act = CAPS_acutation(elsdata,tempdatetime)
    IBSvecs = rotate_CAPS_SCframe(act,'ibs2',anodes=False)
    azimuthvec_norm = [0,-1,0]
    elevationvec_norm = [0,0,-1]

    #IBS_elv = np.arccos(spice.vdot(elevationvec_norm,IBSvecs))*spice.dpr()
    #IBS_azi = np.arcsin(spice.vdot(azimuthvec_norm,IBSvecs))*spice.dpr()
    IBS_elv = 0
    IBS_azi = act
    print(act,IBSvecs,IBS_azi,IBS_elv)

    return IBS_azi, IBS_elv

def CAPS_actuationtimeslice(datetime,elsdata):
    '''
    Returns the start and end slicenumbers while CAPS in actuation in one direction
    '''
    
    #TODO remove ELSdata dependency
    
    slicevalue = CAPS_slicenumber(elsdata,datetime)
    tempslicevalue = slicevalue
    
    if elsdata['actuator'][slicevalue+1] > elsdata['actuator'][slicevalue]:
        direction = "positive"
    else:
        direction = "negative"
        
    if direction == "positive":
        while elsdata['actuator'][tempslicevalue-1] <  elsdata['actuator'][tempslicevalue]:
            tempslicevalue -=1 
        startactslice = tempslicevalue
        tempslicevalue = slicevalue
        while elsdata['actuator'][tempslicevalue] >  elsdata['actuator'][tempslicevalue-1]:
            tempslicevalue +=1 
        endactslice = tempslicevalue
        
    if direction == "negative":
        while elsdata['actuator'][tempslicevalue-1] >  elsdata['actuator'][tempslicevalue]:
            tempslicevalue -=1 
        startactslice = tempslicevalue
        tempslicevalue = slicevalue
        while elsdata['actuator'][tempslicevalue] <  elsdata['actuator'][tempslicevalue-1]:
            tempslicevalue +=1 
        endactslice = tempslicevalue      
        
    return startactslice, endactslice, direction
    
def CAPS_slicenumber(data,tempdatetime):
    
    for counter, i in enumerate(data['times_utc']):
        if i >= tempdatetime:
            slicenumber = counter
            break      
            
    return slicenumber


def CAPS_energyslice(sensor,startenergy,endenergy):
    '''
    Returns the start and end slices of the earray required to cover energy range
    
    '''
    
    if sensor == "els":
        polyearray = elscalib['polyearray']
    if sensor == "ims":
        polyearray = sngcalib['sngpolyearray']
    if sensor == "ibs":
        polyearray = ibscalib['ibspolyearray']
        
    if startenergy < polyearray[0]:
        startcounter = 0
    else:
        for counter, energy in enumerate(polyearray):
            if startenergy >= energy and startenergy < polyearray[counter+1]:
                startcounter = counter
                break
                
    if endenergy > polyearray[-1]:
        endcounter = len(polyearray)-1
    else:                
        for counter, energy in enumerate(polyearray):
            if endenergy >= energy and endenergy < polyearray[counter+1]:
                endcounter = counter
                break
                
    return startcounter,endcounter
    
def ELS_backgroundremoval(data,startslice,endslice):
    
    
    def_backgroundremoved = np.zeros((63,8,endslice-startslice))
    #Background defined as average of 5 loweest count, negative ions unlikely to appear across 3 anodes
    for backgroundcounter, timecounter in enumerate(np.arange(startslice,endslice,1)):
    
        for energycounter in range(63):
            backgroundremoved_temp = np.array(data['def'][energycounter,:8,timecounter]) - np.mean(sorted(data['def'][energycounter,:8,timecounter])[:5])
            backgroundremoved_anodes = [0 if i < 0 else i for i in backgroundremoved_temp] 
            def_backgroundremoved[energycounter,:,backgroundcounter] = backgroundremoved_anodes 
            
    return def_backgroundremoved