def caps_hdf5totimeseries(hdfpath,instrumentanode=0):
    hdf5file = h5py.File('data/ELS_200533900_V01.hdf5','r')
    data_df = pd.DataFrame(hdf5file['DATA'][:,:,instrumentanode])
    times = caps_dateparser([x.decode("ASCII") for x in hdf5file['UTC']])
    data_df = data_df.set_index(pd.DatetimeIndex(times))
    date_timeseries = GenericTimeSeries(data_df)
    
    return(date_timeseries)
    hdf5file.close()
    
caps_hdf5totimeseries('data/ELS_200533900_V01.hdf5')