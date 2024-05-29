# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown] user_expressions=[]
# # Access Sentinel 2 Data from AWS
#
# https://odc-stac.readthedocs.io/en/latest/notebooks/stac-load-e84-aws.html
#
# https://registry.opendata.aws/sentinel-2-l2a-cogs/
#
# Guides:
#
# https://www.matecdev.com/posts/landsat-sentinel-aws-s3-python.html
# https://docs.digitalearthafrica.org/en/latest/platform_tools/odc_stac.html#Connect-to-the-Digital-Earth-Africa-STAC-Catalog
#
# Maybe i need to follow this guide for authentication??
# https://pystac-client.readthedocs.io/en/stable/tutorials/authentication.html
#
# https://buildmedia.readthedocs.org/media/pdf/rasterio/stable/rasterio.pdf
#
# https://boto3.amazonaws.com/v1/documentation/api/latest/guide/configuration.html
#
# https://odc-stac.readthedocs.io/en/latest/_api/odc.stac.load.html
#
# https://snowex-2021.hackweek.io/tutorials/geospatial/raster.html
#

# %%

# %%
# %pip install odc-stac

# %%
import dask.distributed
import folium
import folium.plugins
import geopandas as gpd
import shapely.geometry
from IPython.display import HTML, display
from pystac_client import Client
import os
import boto3
import rasterio as rio
import stackstac

from odc.stac import configure_rio, stac_load


# %%
def convert_bounds(bbox, invert_y=False):
    """
    Helper method for changing bounding box representation to leaflet notation

    ``(lon1, lat1, lon2, lat2) -> ((lat1, lon1), (lat2, lon2))``
    """
    x1, y1, x2, y2 = bbox
    if invert_y:
        y1, y2 = y2, y1
    return ((y1, x1), (y2, x2))

cfg = {
    "sentinel-s2-l1c": {
        "assets": {
            "*": {"data_type": "uint16", "nodata": 0},
            "SCL": {"data_type": "uint8", "nodata": 0},
            "visual": {"data_type": "uint8", "nodata": 0},
        },
        "aliases": {"red": "B04", "green": "B03", "blue": "B02"},
    },
    "*": {"warnings": "ignore"},
}

# %%
os.environ["AWS_REQUEST_PAYER"] = "requester"
os.environ["AWS_PROFILE"] = "akpetty"

# %%


# https://stackoverflow.com/questions/67812512/rasterio-does-not-exist-in-the-file-system-and-is-not-recognized-as-a-support

#os.environ['CURL_CA_BUNDLE'] = '/etc/ssl/certs/ca-certificates.crt'

#print("Creating AWS Session")
#aws_session = rio.session.AWSSession(boto3.Session(), requester_pays=True)

# %%
#configure_rio(
#    cloud_defaults=True,
#    aws={"aws_unsigned": True},
#    AWS_S3_ENDPOINT="s3.af-south-1.amazonaws.com",
#)

#configure_rio(cloud_defaults=True, aws={"aws_unsigned": True, "requester_pays":True})

# %% [markdown]
# ## Start Dask Client
#
# This step is optional, but it does improve load speed significantly. You
# don't have to use Dask, as you can load data directly into memory of the
# notebook.

# %%
#client = dask.distributed.Client()
#configure_rio(cloud_defaults=True, aws={"aws_unsigned": True, "requester_pays":True}, client=client)
#display(client)

# %% [markdown]
# ## Find STAC Items to Load

# %%
km2deg = 1.0 / 111
x, y = (113.887, -25.843)  # Center point of a query
r = 100 * km2deg
bbox = (x - r, y - r, x + r, y + r)

catalog = Client.open("https://earth-search.aws.element84.com/v1")

query = catalog.search(
    collections=["sentinel-2-l1c"], datetime="2021-09-16", limit=10, bbox=bbox
)

stac_items = query.items()

#items = list(query.get_items())
#print(f"Found: {len(items):d} datasets")



# %%

# %%

# %%
# Convert STAC items into a GeoJSON FeatureCollection
stac_json = query.get_all_items_as_dict()

# %%
gdf = gpd.GeoDataFrame.from_features(stac_json, "epsg:4326")
gdf

# %%
#gdf.iloc[0]['earthsearch:s2:datastrip_id']

# %% [markdown] user_expressions=[]
# ## Review Query Result
#
# We'll use GeoPandas DataFrame object to make plotting easier.

# %%

# Compute granule id from components
#gdf["granule"] = (
#    gdf["sentinel:utm_zone"].apply(lambda x: f"{x:02d}")
#    + gdf["sentinel:latitude_band"]
#    + gdf["sentinel:grid_square"]
#)

fig = gdf.plot("s2:granule_id",
    edgecolor="black",
    categorical=True,
    aspect="equal",
    alpha=0.5,
    figsize=(6, 12),
    legend=True,
    legend_kwds={"loc": "upper left", "frameon": False, "ncol": 1},
)
_ = fig.set_title("STAC Query Results")

# %% [markdown]
# ## Plot STAC Items on a Map

# %%
# https://github.com/python-visualization/folium/issues/1501
from branca.element import Figure

fig = Figure(width="400px", height="500px")
map1 = folium.Map()
fig.add_child(map1)

folium.GeoJson(
    shapely.geometry.box(*bbox),
    style_function=lambda x: dict(fill=False, weight=1, opacity=0.7, color="olive"),
    name="Query",
).add_to(map1)

gdf.explore(
    "s2:granule_id",
    categorical=True,
    tooltip=[
        "s2:granule_id",
        "eo:cloud_cover",
    ],
    popup=True,
    style_kwds=dict(fillOpacity=0.1, width=2),
    name="STAC",
    m=map1,
)

map1.fit_bounds(bounds=convert_bounds(gdf.unary_union.bounds))
display(fig)

# %% [markdown]
#
# Note that even though there are 9 STAC Items on input, there is only one
# timeslice on output. This is because of `groupby="solar_day"`. With that
# setting `stac_load` will place all items that occured on the same day (as
# adjusted for the timezone) into one image plane.

# %%
#for stac_asset in items.assets.values():
#    stac_asset.href = stac_asset.href.replace(
#        "s3://sentinel-s2-l2a/", "s3://sentinel-s2-l1c/"
#    )

stac_item = next(stac_items)

for stac_asset in stac_item.assets.values():
    stac_asset.href = stac_asset.href.replace(
        "s3://sentinel-s2-l2a/", "s3://sentinel-s2-l1c/"
    )
    
    
dataarray = stackstac.stack(items=stac_item, dtype="float16", resolution=10)
print(dataarray)

# %%
da_rgb = dataarray.sel(band=["red", "green", "blue"]).squeeze()[:]  # get subset
da_rgb.astype("int").plot.imshow(rgb="band", robust=True)

# %%
dataarray.band

# %%
da_rgb

# %%
_ = (
    da_rgb.plot.imshow(
        col="band",
        size=4,
        vmin=0,
        vmax=4000,
    )
)

# %%
#odc.stac.configure_s3_access(profile=None, region_name='auto', aws_unsigned=None, requester_pays=True, cloud_defaults=True)
#configure_rio(profile=None, region_name='auto', aws_unsigned=None, requester_pays=False, cloud_defaults=True, **gdal_opts)## Construct Dask Dataset

#configure_rio(
#    cloud_defaults=True,
#    aws={"aws_unsigned": True},
#    AWS_S3_ENDPOINT="s3.af-south-1.amazonaws.com",
#)


# %%
# Since we will plot it on a map we need to use `EPSG:3857` projection
crs = "epsg:3857"
zoom = 2**5  # overview level 5

xx = stac_load(
    stac_item,
    bands=("red", "green", "blue"),
    crs=crs,
    resolution=10 * zoom,
    groupby="solar_day",
    chunks={},
    stac_cfg=cfg,
)
display(xx)

# %%
stac_item

# %%
xx

# %% [markdown]
# Note that data is not loaded yet. But we can review memory requirement. We can also check data footprint.

# %%
xx.odc.geobox

# %%
#items[0]

# %%
#gdf["earthsearch:s3_path"][0]

# %%
#import rasterio as rio
#from rasterio.session import AWSSession
#import rioxarray
#import boto3
#session = boto3.Session()
#rio_env = rio.Env(AWSSession(session, requester_pays=True), 
#                  AWS_NO_SIGN_REQUEST='NO',
#                  GDAL_DISABLE_READDIR_ON_OPEN='TRUE')
#rio_env.__enter__()

# %%
#s3_url = 's3://usgs-landsat/collection02/level-1/standard/oli-tirs/2021/016/042/LC08_L1TP_016042_20211027_20211104_02_T1/LC08_L1TP_016042_20211027_20211104_02_T1_B1.TIF'
#da = rioxarray.open_rasterio(s3_url, chunks='auto').squeeze('band', drop=True)
#da

# %%
#with rio.Env(rio_env):
#    da = rioxarray.open_rasterio(gdf["earthsearch:s3_path"][0], chunks='auto').squeeze('band', drop=True)


# %%
#gdf["earthsearch:s3_path"][0]

# %%
#da = rioxarray.open_rasterio(gdf["earthsearch:s3_path"][0], chunks='auto').squeeze('band', drop=True)


# %%
#xx.load()

# %%

# %% [markdown]
# ## Load data into local memory

# %%
# %%time
xx2 = xx.compute()

# %%
_ = (
    xx.isel(time=0)
    .to_array("band")
    .plot.imshow(
        col="band",
        size=4,
        vmin=0,
        vmax=4000,
    )
)

# %% [markdown]
# ## Load with bounding box
#
# As you can see `stac_load` returned all the data covered by STAC items
# returned from the query. This happens by default as `stac_load` has no way of
# knowing what your query was. But it is possible to control what region is
# loaded. There are several mechanisms available, but probably simplest one is
# to use `bbox=` parameter (compatible with `stac_client`).
#
# Let's load a small region at native resolution to demonstrate.

# %%
r = 6.5 * km2deg
small_bbox = (x - r, y - r, x + r, y + r)

yy = stac_load(
    items,
    bands=("red", "green", "blue"),
    crs=crs,
    resolution=10,
    chunks={},  # <-- use Dask
    groupby="solar_day",
    stac_cfg=cfg,
    bbox=small_bbox,
)
display(yy.odc.geobox)

# %%
yy = yy.compute()

# %%
_ = (
    yy.isel(time=0)
    .to_array("band")
    .plot.imshow(
        col="band",
        size=4,
        vmin=0,
        vmax=4000,
    )
)

# %% [markdown]
# --------------------------------------------------------------
