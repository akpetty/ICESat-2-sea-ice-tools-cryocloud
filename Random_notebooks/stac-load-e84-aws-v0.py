# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown] user_expressions=[]
# # Access Sentinel 2 Data from AWS
#
# [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/opendatacube/odc-stac/develop?labpath=notebooks%2Fstac-load-e84-aws.ipynb)
#
# https://registry.opendata.aws/sentinel-2-l2a-cogs/
#
# https://odc-stac.readthedocs.io/en/latest/notebooks/stac-load-e84-aws.html
#

# %% tags=[]
import dask.distributed
import folium
import folium.plugins
import geopandas as gpd
import shapely.geometry
from IPython.display import HTML, display
from pystac_client import Client

from odc.stac import configure_rio, stac_load


def convert_bounds(bbox, invert_y=False):
    """
    Helper method for changing bounding box representation to leaflet notation

    ``(lon1, lat1, lon2, lat2) -> ((lat1, lon1), (lat2, lon2))``
    """
    x1, y1, x2, y2 = bbox
    if invert_y:
        y1, y2 = y2, y1
    return ((y1, x1), (y2, x2))


# %% tags=[]
cfg = {
    "sentinel-s2-l2a-cogs": {
        "assets": {
            "*": {"data_type": "uint16", "nodata": 0},
            "SCL": {"data_type": "uint8", "nodata": 0},
            "visual": {"data_type": "uint8", "nodata": 0},
        },
        "aliases": {"red": "B04", "green": "B03", "blue": "B02"},
    },
    "*": {"warnings": "ignore"},
}

# %% [markdown] user_expressions=[]
# ## Start Dask Client
#
# This step is optional, but it does improve load speed significantly. You
# don't have to use Dask, as you can load data directly into memory of the
# notebook.

# %%
#client = dask.distributed.Client()
#configure_rio(cloud_defaults=True, aws={"aws_unsigned": True}, client=client)
#display(client)

# %% [markdown] user_expressions=[]
# ## Find STAC Items to Load

# %% tags=[]
km2deg = 1.0 / 111
x, y = (113.887, -25.843)  # Center point of a query
r = 100 * km2deg
bbox = (x - r, y - r, x + r, y + r)

catalog = Client.open("https://earth-search.aws.element84.com/v0")

query = catalog.search(
    collections=["sentinel-s2-l2a-cogs"], datetime="2021-09-16", limit=100, bbox=bbox
)

items = list(query.get_items())
print(f"Found: {len(items):d} datasets")

# Convert STAC items into a GeoJSON FeatureCollection
stac_json = query.get_all_items_as_dict()

# %% tags=[]
query = catalog.search(datetime="2021-09-16", limit=100, bbox=bbox)
items = list(query.get_items())
print(f"Found: {len(items):d} datasets")
items

# %% [markdown] user_expressions=[] tags=[]
# ## Review Query Result
#
# We'll use GeoPandas DataFrame object to make plotting easier.

# %% tags=[]
gdf = gpd.GeoDataFrame.from_features(stac_json, "epsg:4326")

# Compute granule id from components
gdf["granule"] = (
    gdf["sentinel:utm_zone"].apply(lambda x: f"{x:02d}")
    + gdf["sentinel:latitude_band"]
    + gdf["sentinel:grid_square"]
)

fig = gdf.plot(
    "granule",
    edgecolor="black",
    categorical=True,
    aspect="equal",
    alpha=0.5,
    figsize=(6, 12),
    legend=True,
    legend_kwds={"loc": "upper left", "frameon": False, "ncol": 1},
)
_ = fig.set_title("STAC Query Results")

# %% [markdown] user_expressions=[]
# ## Plot STAC Items on a Map

# %% tags=[]
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
    "granule",
    categorical=True,
    tooltip=[
        "granule",
        "datetime",
        "sentinel:data_coverage",
        "eo:cloud_cover",
    ],
    popup=True,
    style_kwds=dict(fillOpacity=0.1, width=2),
    name="STAC",
    m=map1,
)

map1.fit_bounds(bounds=convert_bounds(gdf.unary_union.bounds))
display(fig)

# %% [markdown] user_expressions=[]
# ## Construct Dask Dataset
#
# Note that even though there are 9 STAC Items on input, there is only one
# timeslice on output. This is because of `groupby="solar_day"`. With that
# setting `stac_load` will place all items that occured on the same day (as
# adjusted for the timezone) into one image plane.

# %% tags=[]
# Since we will plot it on a map we need to use `EPSG:3857` projection
crs = "epsg:3857"
zoom = 2**5  # overview level 5

xx = stac_load(
    items,
    bands=("red", "green", "blue"),
    crs=crs,
    resolution=10 * zoom,
    chunks={},  # <-- use Dask
    groupby="solar_day",
    stac_cfg=cfg,
)
display(xx)

# %% [markdown] user_expressions=[]
# Note that data is not loaded yet. But we can review memory requirement. We can also check data footprint.

# %% tags=[]
xx.odc.geobox

# %% [markdown] user_expressions=[]
# ## Load data into local memory

# %% tags=[]
# %%time
xx = xx.compute()

# %% tags=[]
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

# %% [markdown] user_expressions=[]
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
km2deg = 1.0 / 111
x, y = (113.887, -25.843)  # Center point of a query
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

# %% [markdown] user_expressions=[] tags=[]
# --------------------------------------------------------------
