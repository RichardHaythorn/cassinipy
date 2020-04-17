from sunpy.timeseries import TimeSeriesMetaData
from sunpy.util import MetaDict

from sunpy.time import TimeRange

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
        
        capstimeseries = sunpy.timeseries.GenericTimeSeries(dataframe,meta=caps_metadata)
        
    if 'sngdata' in data.keys():
        print("ims")
        instrument="ims"
    if 'ibsdata' in data.keys():
        print("ibs")
        instrument="ibs"
    
    return capstimeseries
        
   # data_metadata = TimeSeriesMetaData(meta=(instrument))   
    
    #return sunpy.timeseries.GenericTimeSeries(meta=data_metadata)
        