import logging
from pathlib import Path
from typing import List

import rioxarray
import xarray as xr
from pyproj import Transformer
import geopandas as gpd

_log = logging.getLogger(__name__)


def _get_dimension(dims: dict, candidates: List[str]):
    for name in candidates:
        if name in dims:
            return name
    error = f'Dimension matching one of the candidates {candidates} not found! The available ones are {dims}. Please rename the dimension accordingly and try again. This local collection will be skipped.'
    raise Exception(error)


def _get_netcdf_zarr_metadata(file_path):
    if '.zarr' in file_path.suffixes:
        data = xr.open_dataset(file_path.as_posix(),chunks={},engine='zarr')
    else:
        data = xr.open_dataset(file_path.as_posix(),chunks={}) # Add decode_coords='all' if the crs as a band gives some issues
    file_path = file_path.as_posix()
    try:
        try:
            t_dim = _get_dimension(data.dims, ['t', 'time', 'temporal', 'DATE'])
        except Exception as e:
            t_dim = None
            _log.warning(e)
        try:
            x_dim = _get_dimension(data.dims, ['x', 'X', 'lon', 'longitude'])
        except Exception as e:
            x_dim = None
            _log.warning(e)
        try:
            y_dim = _get_dimension(data.dims, ['y', 'Y', 'lat', 'latitude'])
        except Exception as e:
            y_dim = None
            _log.warning(e)
        try:
            geometry_dim = _get_dimension(data.dims, ['cells', 'cell_id', 'cell_ids', 'geometry'])
        except Exception as e:
            geometry_dim = None
            _log.warning(e)
    except Exception as e:
        _log.warning(e)
        raise Exception(f'Error creating metadata for {file_path}') from e
    metadata = {}
    metadata['stac_version'] = '1.0.0'
    metadata['type'] = 'Collection'
    metadata['id'] = file_path
    data_attrs_lowercase = [x.lower() for x in data.attrs]
    data_attrs_original  = [x for x in data.attrs]
    data_attrs = dict(zip(data_attrs_lowercase,data_attrs_original))
    if 'title' in data_attrs_lowercase:
        metadata['title'] = data.attrs[data_attrs['title']]
    else:
        metadata['title'] = file_path
    if 'description' in data_attrs_lowercase:
        metadata['description'] = data.attrs[data_attrs['description']]
    else:
        metadata['description'] = ''
    if 'license' in data_attrs_lowercase:
        metadata['license'] = data.attrs[data_attrs['license']]
    else:
        metadata['license'] = ''
    providers = [{'name':'',
                 'roles':['producer'],
                 'url':''}]
    if 'providers' in data_attrs_lowercase:
        providers[0]['name'] = data.attrs[data_attrs['providers']]
        metadata['providers'] = providers
    elif 'institution' in data_attrs_lowercase:
        providers[0]['name'] = data.attrs[data_attrs['institution']]
        metadata['providers'] = providers
    else:
        metadata['providers'] = providers
    if 'links' in data_attrs_lowercase:
        metadata['links'] = data.attrs[data_attrs['links']]
    else:
        metadata['links'] = ''
    
    if x_dim is not None and y_dim is not None:
        x_min = data[x_dim].min().item(0)
        x_max = data[x_dim].max().item(0)
        y_min = data[y_dim].min().item(0)
        y_max = data[y_dim].max().item(0)
    elif geometry_dim is not None:
        try:
            import xdggs
            data = data.pipe(xdggs.decode)
            # cell_boundaries = data.dggs.cell_boundaries()
            # bounds = gpd.GeoSeries(cell_boundaries.values).total_bounds
            # x_min = bounds[0]
            # x_max = bounds[2]
            # y_min = bounds[1]
            # y_max = bounds[3]
            
            cell_centers = data.dggs.cell_centers()

            x_min = cell_centers.latitude.min().values.item()
            x_max = cell_centers.latitude.max().values.item()
            y_min = cell_centers.longitude.min().values.item() - 360
            y_max = cell_centers.longitude.max().values.item() - 360

            data.rio.write_crs(4326,inplace=True)
        except Exception as e:
            # raise Exception("Could not extract the spatial extent of the data!")
            raise Exception(e)

    crs_present = False
    bands = list(data.data_vars)
    if 'crs' in bands:
        crs = data.crs.attrs["crs_wkt"]
        bands.remove('crs')
        crs_present = True
    elif data.rio.crs:
        crs = data.rio.crs
        crs_present = True
    extent = {}
    if crs_present:
            transformer = Transformer.from_crs(crs, "epsg:4326")
            lat_min, lon_min = transformer.transform(x_min, y_min)
            lat_max, lon_max = transformer.transform(x_max, y_max)
            extent["spatial"] = {"bbox": [[lon_min, lat_min, lon_max, lat_max]]}
    if t_dim is not None:
        t_min = str(data[t_dim].min().values)
        t_max = str(data[t_dim].max().values)
        extent['temporal'] = {'interval': [[t_min,t_max]]}

    metadata['extent'] = extent

    t_dimension = {}
    if t_dim is not None:
        t_dimension = {t_dim: {'type': 'temporal', 'extent':[t_min,t_max]}}

    x_dimension = {}
    if x_dim is not None:
        x_dimension = {x_dim: {'type': 'spatial','axis':'x','extent':[x_min,x_max]}}
        if crs_present:
            x_dimension[x_dim]['reference_system'] = str(crs)

    y_dimension = {}
    if y_dim is not None:
        y_dimension = {y_dim: {'type': 'spatial','axis':'y','extent':[y_min,y_max]}}
        if crs_present:
            y_dimension[y_dim]['reference_system'] = str(crs)

    geometry_dimension = {}
    if geometry_dim is not None:
        geometry_dimension = {geometry_dim: {'type': 'geometry','bbox': [[y_min, x_min, y_max, x_max]]}}

    b_dimension = {}
    if len(bands)>0:
        b_dimension = {'bands': {'type': 'bands', 'values':bands}}

    metadata['cube:dimensions'] = {**t_dimension,**x_dimension,**y_dimension,**geometry_dimension,**b_dimension}

    return metadata


def _get_geotiff_metadata(file_path):
    data = rioxarray.open_rasterio(file_path.as_posix(),chunks={},band_as_variable=True)
    file_path = file_path.as_posix()
    try:
        t_dim = _get_dimension(data.dims, ['t', 'time', 'temporal', 'DATE'])
    except Exception:
        t_dim = None
    try:
        x_dim = _get_dimension(data.dims, ['x', 'X', 'lon', 'longitude'])
        y_dim = _get_dimension(data.dims, ['y', 'Y', 'lat', 'latitude'])
    except Exception as e:
        _log.warning(e)
        raise Exception(f'Error creating metadata for {file_path}') from e

    metadata = {}
    metadata['stac_version'] = '1.0.0-rc.2'
    metadata['type'] = 'Collection'
    metadata['id'] = file_path
    data_attrs_lowercase = [x.lower() for x in data.attrs]
    data_attrs_original  = [x for x in data.attrs]
    data_attrs = dict(zip(data_attrs_lowercase,data_attrs_original))
    if 'title' in data_attrs_lowercase:
        metadata['title'] = data.attrs[data_attrs['title']]
    else:
        metadata['title'] = file_path
    if 'description' in data_attrs_lowercase:
        metadata['description'] = data.attrs[data_attrs['description']]
    else:
        metadata['description'] = ''
    if 'license' in data_attrs_lowercase:
        metadata['license'] = data.attrs[data_attrs['license']]
    else:
        metadata['license'] = ''
    providers = [{'name':'',
                 'roles':['producer'],
                 'url':''}]
    if 'providers' in data_attrs_lowercase:
        providers[0]['name'] = data.attrs[data_attrs['providers']]
        metadata['providers'] = providers
    elif 'institution' in data_attrs_lowercase:
        providers[0]['name'] = data.attrs[data_attrs['institution']]
        metadata['providers'] = providers
    else:
        metadata['providers'] = providers
    if 'links' in data_attrs_lowercase:
        metadata['links'] = data.attrs[data_attrs['links']]
    else:
        metadata['links'] = ''
    x_min = data[x_dim].min().item(0)
    x_max = data[x_dim].max().item(0)
    y_min = data[y_dim].min().item(0)
    y_max = data[y_dim].max().item(0)

    crs_present = False
    coords = list(data.coords)
    if 'spatial_ref' in coords:
        # bands.remove('crs')
        crs_present = True
    bands = []
    for d in data.data_vars:
        data_attrs_lowercase = [x.lower() for x in data[d].attrs]
        data_attrs_original  = [x for x in data[d].attrs]
        data_attrs = dict(zip(data_attrs_lowercase,data_attrs_original))
        if 'description' in data_attrs_lowercase:
            bands.append(data[d].attrs[data_attrs['description']])
        else:
            bands.append(d)
    extent = {}
    if crs_present:
        if 'crs_wkt' in data.spatial_ref.attrs:
            transformer = Transformer.from_crs(data.spatial_ref.attrs['crs_wkt'], 'epsg:4326')
            lat_min,lon_min = transformer.transform(x_min,y_min)
            lat_max,lon_max = transformer.transform(x_max,y_max)
            extent['spatial'] = {'bbox': [[lon_min, lat_min, lon_max, lat_max]]}

    if t_dim is not None:
        t_min = str(data[t_dim].min().values)
        t_max = str(data[t_dim].max().values)
        extent['temporal'] = {'interval': [[t_min,t_max]]}

    metadata['extent'] = extent

    t_dimension = {}
    if t_dim is not None:
        t_dimension = {t_dim: {'type': 'temporal', 'extent':[t_min,t_max]}}

    x_dimension = {x_dim: {'type': 'spatial','axis':'x','extent':[x_min,x_max]}}
    y_dimension = {y_dim: {'type': 'spatial','axis':'y','extent':[y_min,y_max]}}
    if crs_present:
        if 'crs_wkt' in data.spatial_ref.attrs:
            x_dimension[x_dim]['reference_system'] = data.spatial_ref.attrs['crs_wkt']
            y_dimension[y_dim]['reference_system'] = data.spatial_ref.attrs['crs_wkt']

    b_dimension = {}
    if len(bands)>0:
        b_dimension = {'bands': {'type': 'bands', 'values':bands}}

    metadata['cube:dimensions'] = {**t_dimension,**x_dimension,**y_dimension,**b_dimension}

    return metadata


def _get_local_collections(local_collections_path):
    if isinstance(local_collections_path,str):
        local_collections_path = [local_collections_path]
    local_collections_list = []
    for flds in local_collections_path:
        local_collections_netcdf_zarr = [p for p in Path(flds).rglob('*') if p.suffix in  ['.nc','.zarr']]
        for local_file in local_collections_netcdf_zarr:
            try:
                metadata = _get_netcdf_zarr_metadata(local_file)
                local_collections_list.append(metadata)
            except Exception as e:
                _log.error(e)
                continue
        local_collections_geotiffs = [p for p in Path(flds).rglob('*') if p.suffix in  ['.tif','.tiff']]
        for local_file in local_collections_geotiffs:
            try:
                metadata = _get_geotiff_metadata(local_file)
                local_collections_list.append(metadata)
            except Exception as e:
                _log.error(e)
                continue
    local_collections_dict = {'collections':local_collections_list}

    return local_collections_dict
