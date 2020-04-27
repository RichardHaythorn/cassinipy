def caps_hdf5totimeseries(hdfpath,instrumentanode=0):
    '''
    Creates a timeseries for a single anode/fan for a CAPS sensor
    '''   
    #TODO add the metadata
    #TODO add the units 
    hdf5file = h5py.File(hdfpath,'r')
    data_df = pd.DataFrame(hdf5file['DATA'][:,:,instrumentanode])
    times = caps_dateparser([x.decode("ASCII") for x in hdf5file['UTC']])
    data_df = data_df.set_index(pd.DatetimeIndex(times))
    date_timeseries = GenericTimeSeries(data_df)
    hdf5file.close()
    return(date_timeseries)
    
    