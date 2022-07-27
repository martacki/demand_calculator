import pandas as pd
import xarray as xr
import datetime
import numpy as np
import matplotlib.pyplot as plt
import random as rd

alpha3_to_alpha2 = {
    'AUT': 'AT', 'UKR': 'UK', 'BIH': 'BA', 'BEL': 'BE',
    'BGR': 'BG', 'CHE': 'CH', 'CYP': 'CY', 'CZE': 'CZ',
    'DEU': 'DE', 'DNK': 'DK', 'EST': 'EE', 'GRC': 'GR',
    'ESP': 'ES', 'FIN': 'FI', 'FRA': 'FR', 'HRV': 'HR',
    'HUN': 'HU', 'IRL': 'IR', 'ISL': 'IS', 'ITA': 'IT',
    'LTU': 'LT', 'LUX': 'LU', 'LVA': 'LV', 'MDA': 'MD',
    'MNE': 'ME', 'MKD': 'MK', 'NLD': 'NL', 'NOR': 'NO',
    'POL': 'PL', 'PRT': 'PT', 'ROU': 'RO', 'SRB': 'RS',
    'SWE': 'SE', 'SVN': 'SI', 'SVK': 'SK', 'GBR': 'GB',
}

reference_files = snakemake.input.reference

demand_ref = pd.read_csv(reference_files[0], index_col=0)
for year_fn in reference_files[1:]:
    demand_ref = pd.concat([demand_ref, pd.read_csv(year_fn, index_col=0)])
demand_ref['month'] = (
    demand_ref.reset_index().utc_timestamp.apply(lambda b: b.split(' ')[0].split('-')[1]).values
)

# renamde index: date -> weekday
new_idx = [datetime.datetime(int(x.split(' ')[0].split('-')[0]),
                             int(x.split(' ')[0].split('-')[1]),
                             int(x.split(' ')[0].split('-')[2]))
           .strftime('%A') for x in demand_ref.index]
demand_ref.index = new_idx

res_year = snakemake.wildcards.yr
result = xr.open_dataset(snakemake.input.demand)
result = result.sortby('time')

daterange = pd.date_range(
    start=f"{res_year[0:4]}-01-01", end=f"{int(res_year[0:4])+1}-01-01", inclusive="left", freq="1H"
)

demand_timeseries = pd.DataFrame(index=daterange, columns=demand_ref.columns[:-1], dtype='object')
for doy in range(result.dims['time']): #doy = day of year
    result_month = str(result.time[doy].values).split('T')[0].split('-')[1]
    result_day = datetime.datetime(int(str(result.time[doy].values).split('T')[0].split('-')[0]),
                                   int(str(result.time[doy].values).split('T')[0].split('-')[1]),
                                   int(str(result.time[doy].values).split('T')[0].split('-')[2]),
                                  ).strftime('%A')

    demand_ref_day = demand_ref.query("month == @result_month").loc[result_day]
    demand_ref_day.month = demand_ref_day.month.astype(int)

    random_sel = rd.randint(0, len(demand_ref_day)/24-1)
    ref_idx = [(random_sel)*24, (random_sel+1)*24]
    demand_ref_day = demand_ref_day.iloc[ref_idx[0]:ref_idx[1]]
    demand_ref_day_norm = demand_ref_day / demand_ref_day.max()
    demand_ref_day_norm /= demand_ref_day_norm.sum()
    
    daily_demand = pd.Series(index=result.country, data=result.demand[doy].values)*1e3
    daily_demand.index = (
        pd.DataFrame(daily_demand).reset_index()['index'].map(alpha3_to_alpha2).values
    )
    demand_timeseries.iloc[doy*24:(doy+1)*24] = (
        (daily_demand * demand_ref_day_norm)[demand_timeseries.columns]
    )

demand_timeseries = demand_timeseries.fillna(0)
demand_timeseries.to_csv(snakemake.output[0])
