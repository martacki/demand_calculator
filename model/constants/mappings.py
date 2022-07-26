import pandas as pd
import numpy as np
import config

### europe area 
l = config.l
r = config.r
t = config.t
b = config.b
area = dict(lon=slice(l,r), lat=slice(t, b))

l0=['ALB','AUT','UKR','BIH','BEL','BGR','CHE','CYP','CZE','DEU','DNK','EST','GRC','ESP','FIN',
    'FRA','HRV','HUN','IRL','ISL','ITA','LIE','LTU','LUX','LVA','MDA','MNE','MKD','MLT','NLD',
    'NOR','POL','PRT','ROU','SRB','SWE','SVN','SVK','TUR','GBR', 'GBR','XKX'
   ]
l1=[ 'AL', 'AT',   '', 'BA', 'BE', 'BG', 'CH', 'CY', 'CZ', 'DE', 'DK', 'EE', 'EL', 'ES', 'FI',
    'FR', 'HR', 'HU', 'IE', 'IS', 'IT', 'LI', 'LT', 'LU', 'LV',   '', 'ME', 'MK', 'MT', 'NL', 
    'NO', 'PL', 'PT', 'RO', 'RS', 'SE', 'SI', 'SK', 'TR', 'UK',    '','XK'
   ]
l2=[   '', 'AT',   '', 'BA', 'BE', 'BG', 'CH', 'CY', 'CZ', 'DE', 'DK', 'EE', 'GR', 'ES', 'FI',
    'FR', 'HR', 'HU', 'IE', 'IS', 'IT', 'LI', 'LT', 'LU', 'LV',   '', 'ME', 'MK', 'MT', 'NL', 
    'NO', 'PL', 'PT', 'RO',   '', 'SE', 'SI', 'SK', 'TR', 'GB',  'NI',  ''
   ]
l3=[ 'AL', 'AT', 'UA', 'BA', 'BE', 'BG', 'CH', 'CY', 'CZ', 'DE', 'DK', 'EE', 'GR', 'ES', 'FI',
    'FR', 'HR', 'HU', 'IE', 'IS', 'IT',   '', 'LT', 'LU', 'LV', 'MD', 'ME', 'MK', 'MT', 'NL',
    'NO', 'PL', 'PT', 'RO', 'RS', 'SE', 'SI', 'SK', 'TR', 'GB',    '',  ''
   ]
l4 = ['ALB','AUT','UKR','BIH','BEL','BGR','CHE','CYP','CZE','DEU','DNK','EST','GRC','ESP','FIN',
      'FRA','HRV','HUN','IRL','ISL','ITA','LIE','LTU','LUX','LVA','MDA','MNE','MKD','MLT','NLD',
      'NOR','POL','PRT','ROU','SRB','SWE','SVN','SVK','TUR','GBR',    '', 'KOS'
     ]
names = ['Albania', 'Austria', 'Ukraine', 'Bosnia and Herzegovina', 'Belgium', 'Bulgaria', 
         'Switzerland', 'Cyprus', 'Czechia', 'Germany', 'Denmark', 'Estonia', 'Greece', 'Spain',
         'Finland', 'France', 'Croatia', 'Hungary', 'Ireland', 'Iceland', 'Italy', 
         'Liechtenstein', 'Lithuania', 'Luxembourg', 'Latvia', 'Moldova, Republic of', 
         'Montenegro', 'Macedonia, the former Yugoslav Republic of', 'Malta', 'Netherlands', 
         'Norway', 'Poland', 'Portugal', 'Romania', 'Serbia', 'Sweden', 'Slovenia', 'Slovakia', 
         'Turkey', 'United Kingdom', 'United Kingdom', 'Kosovo'
        ]


### country mapping dataframe
df_countries = pd.DataFrame(
    data = dict(
        ISO31661A3 = l0,
        EU_map = l4,
        nuts_id = l1,
        entsoe = l2,
        entsoe_transparency = l3,
        name=names,
    )
).rename_axis('index_nr').reset_index().set_index('ISO31661A3')

# def save_country_dataset(ofile, dataset = df_countries):
#     pd.to_csv(ofile)


if config.selected_countries == 'europe':
    country_selection = df_countries.EU_map.dropna()
else: 
    country_selection = selected_countries
df_countries_select = df_countries.loc[country_selection.index.unique()]
df_countries_select.replace("", np.nan, inplace=True)
#------------------------------------------------------------------------------
## Hydropower database mapping
#------------------------------------------------------------------------------

data_vars = ['name', 'id', 'fuel_type', 'technology_type', 'country', 'country_code',
             'lat', 'lon', 'plant_capacity_MW', 'dam_height_m', 'pump_capacity_MW', 
             'storage_capacity_MWh', 'volume_capacity_Mm3', 'avg_annual_generation_GWh', 
             'efficiency', 'geo_id', 'wri_id', 'jrc_id', 'pypsa_id', 'gpd_id', 
             'entsoe_id', 'opsd_id', 'carma_id', 'set', 'date_in', 'date_retrofit', 'duration'
            ]

jrc_mapping = ['name', 'id', 'Fueltype', 'type', 'country',  'country_code','lat',
               'lon', 'installed_capacity_MW', 'dam_height_m', 'pumping_MW', 
               'storage_capacity_MWh', 'Volume_Mm3', 'avg_annual_generation_GWh', '', 
               'GEO', 'WRI', '', 'pypsa_id', '', '', '', '', '', '', '', ''
              ]

ppm_mapping = ['Name', '', 'Fueltype', 'Technology', 'Country', 'country_code','lat', 'lon',
               'Capacity', 'DamHeight_m', '', '', 'Volume_Mm3', '', 'Efficiency', 'GEO_id', 
               '', 'JRC_id', '', 'GPD_id', 'ENTSOE_id', 'OPSD_id', 'CARMA_id', 'Set', 'DateIn', 
               'DateRetrofit', 'Duration'
              ]

df_hydro = pd.DataFrame(data= {'powerplantmatching' : ppm_mapping, 'jrc':jrc_mapping},
                  index = data_vars)
df_hydro.replace('', np.nan, inplace=True)