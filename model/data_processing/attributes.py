from datetime import datetime
import config

def set_global_attributes(dataset, source, grid='gaussian n80', area='EU13+2'):
    """
    set global attributes to xarray.dataset

    parameters
    ----------
    dataset (xarray.Dataset): dataset for which to update attributes
    parameter (string): option to set for different variables 'demand', 'etc..'  
    
    returns
    -------
    dataset (xarray.DataSet) : dataSet with update attributes
    """

    dataset.attrs.update(
        author = config.author_name,
        project = config.project_name,
        source = source,
        history = f'Computed {datetime.now().strftime("%d-%b-%Y (%H:%M)")}',
        area = area,
        grid = grid,
    )
    
    return dataset


def set_demand_attributes(dataset):
    dataset.demand.attrs.update(
        standard_name = 'demand',
        long_name = 'demand computed from weighted T',
        units = 'GWh',
    )
    return dataset
