import pandas as pd
import numpy as np
import config

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

if config.selected_countries == 'europe':
    country_selection = df_countries.EU_map.dropna()
else: 
    country_selection = config.selected_countries
df_countries_select = df_countries.loc[country_selection.index.unique()]
df_countries_select.replace("", np.nan, inplace=True)
