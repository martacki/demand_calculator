import xarray as xr

def set_lat_lon_attributes(dataset, names=['lat', 'lon']):
    """
    set attributes to lat and lon dimensions of xarray.dataset

    parameters
    ----------
    dataset (xarray.Dataset): dataset for which to update attributes
    parameter (string): option to set for different variables 'demand', 'etc..'  
    
    returns
    -------
    dataset (xarray.DataSet) : dataSet with update attributes
    """
    
    dataset = dataset.rename({names[0]: 'lat', names[1]: 'lon'})
    dataset.lat.attrs.update(
        standard_name = 'latitude',
        long_name = 'latitude',
        units = 'degrees_north',
        axis = 'Y'
    )
    
    dataset.lon.attrs.update(
        standard_name = 'longitude',
        long_name = 'longitude',
        units = 'degrees_east',
        axis = 'X'
    )
    
    return dataset


year = snakemake.wildcards.yr

rename_dict = {'x': 'lon', 'y': 'lat'}

temp = xr.open_dataset(snakemake.input[0])
temp = temp.drop(['height', 'wnd100m', 'wnd_azimuth', 'roughness',
                  'influx_toa', 'influx_direct', 'influx_diffuse',
                  'albedo', 'solar_altitude', 'solar_azimuth',
                  'soil temperature', 'runoff', 'lat', 'lon'])
temp = temp.rename_dims(rename_dict).rename(rename_dict)
temp = temp.resample(time="1D").mean(dim='time')
temp = set_lat_lon_attributes(temp)
temp.to_netcdf(snakemake.output[0])
