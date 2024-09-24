# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.16.2
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown] user_expressions=[]
# # Primer on using odc-stac/odc-geo
#
# With the new releases of `odc-geo` and `odc-stac` there are a lot of new features that make working with geospatial data in a notebook much more pleasant and productive. This notebook walks through some of the features
#
# ```bash
# pip install odc-geo==0.2.0 odc-stac==0.3.0
# # OR
# mambba install odc-geo==0.2.0 odc-stac==0.3.0
# ```
#
# This notebook was tested on [Planetary Computer Hub](https://planetarycomputer.microsoft.com/). It loads 16-days worth of Landsat8 scenes in a UTM zone 33 (covers parts of Europe and Africa). Imagery is loaded at 1km per pixel resolution in Web Mercator projection. We convert pixels to `RGBA<uint8>` by clipping `6_000, 30_000 -> 0, 255` and display result on the map.
#
# You can download [this notebook](odc-stac-demo.ipynb) to try it out on your own.

# %%
# %pip install planetary-computer

# %% [markdown] user_expressions=[]
# ## Query STAC Items
#
# Query 16 days at the start of June 2021. Landsat 8 revisit is 16 days so we should get one observation from each scene in the queried region.

# %%
import pystac_client

catalog = pystac_client.Client.open("https://planetarycomputer.microsoft.com/api/stac/v1")

items_all = catalog.search(
    collections=["landsat-8-c2-l2"],
    datetime="2021-06-01T00:00:00Z/2021-06-16T23:59:59Z",
    bbox=[10, -71, 21, 71],
).get_all_items()

# only keep items in utm=33, to select one tall column of data
items = [item for item in items_all if item.properties["proj:epsg"] in [32633, 32733]]
print(f"Using {len(items)} of {len(items_all)} items (utm 33)")

# %% [markdown] user_expressions=[]
# ## Construct Dask Graph
#
# Request to load data with `odc.stac.load`.
#
# - Only load `red`,`green`,`blue` bands
# - 1km sized pixels in Web Mercator projection (`EPSG:3857`)
# - Specify `uint16` pixels `0` marking nodata pixels (this can also be extracted from STAC extension data if configured)
# - Use `average` resampling for converting to desired projection
# - Put all items into one pixel plane, using custom groupby (`groupby=`)
#   - Most of the time you'll use `groupby="time"` or `groupby="solar_day"`
# - Configure Dask chunking parameters and enable delayed computation (default is to load directly) `chunks=`
#   - This data is tall an narrow so we only chunk tall dimension
# - Supply url signing function `patch_url=planetary_computer.sign`
#   - Data is publicly accessible, but requests need to be signed
# - Convert to RGBa, replacing missing data with transparent pixels
#
# No pixel loading happens yet, this just builds a recipe that we will execute later on.

# %%
import planetary_computer
from dask.distributed import Client
from odc.stac import configure_rio
from odc.stac import load as stac_load

configure_rio(cloud_defaults=True)  # No need to look for side-car files

# fmt: off
xx = stac_load(
    items,
    ["red", "green", "blue"],  # Can be omited, then you'll get all the bancs
    #################################################################
    # All options below have reasonable defaults derived
    # from supplied metadata.
    #################################################################
    
    crs="epsg:3857",       # Defaults to most common CRS across Items
    resolution=1000,       # Defaults to native resolution

    dtype="uint16",        # Defaults to native pixel type
    nodata=0,              # Defaults to native
    
    resampling="average",  # Default is `nearest`
    
    # Default is to group by time
    # groupby="solar_day", # is often what one wants
    groupby=lambda item, parsed, idx: 1, # all items in the same "group"
    
    # Enable and configure Dask
    # since each pixel is 1x1km we opt for smaller chunk sizes
    # Image will be Tall and narrow, so only chunk along Y
    chunks={"x": -1, "y": 1000},
    
    # Needed for planetary computer data sources
    patch_url=planetary_computer.sign,
)
# fmt: on

thumb = xx.odc.to_rgba(vmin=6_000, vmax=30_000)

display(xx)
display(thumb)

# %%
# Start local Dask client to run across CPUs
#  one can also connect to remote cluster
dask_client = Client()

# %% [markdown] user_expressions=[]
# ### Actually Load and Convert to RGBA
#
# Should take about a minute, longer if running from "home".

# %%
# %%time
thumb = thumb.compute()

# %% [markdown] user_expressions=[]
# ## Plot Color Image on a Map
#
# We use `folium` here, but same works with `ipyleaflet` too.
#
# **note**: we are using `webp` compression here to reduce size of rendered notebook. If you don't see an image on a map below it could be that your browser doesn't support webp.

# %%
import folium

_map = folium.Map()
# fmt: off
thumb.odc.add_to(
    _map,
    fmt="webp",        # Default is png, jpeg is also supported
    quality=60,        # webp/jpeg quality configuration
    max_size=10_000,  # disable further shrinking, keep 1km resolution
    opacity=0.90,      # folium layer options are passed-through
)
# fmt: on
_map.fit_bounds(thumb.odc.map_bounds())
display(_map)

# %% [markdown] user_expressions=[]
# ## Add Metadata to the Map
#
# We can also use `geopandas` to plot stac item metadata on top of imagery. Handy for narrowing down test scenes when looking for experiment sites.

# %%
import geopandas
import pystac

items = pystac.ItemCollection(sorted(items, key=lambda x: x.datetime), clone_items=False)
df = geopandas.GeoDataFrame.from_features(items.to_dict(), crs="epsg:4326")
# https://github.com/geopandas/geopandas/issues/1208
df["id"] = [x.id for x in items]
df[["geometry", "id", "datetime", "landsat:wrs_path", "landsat:wrs_row", "proj:epsg"]].explore(
    "landsat:wrs_path",
    style_kwds=dict(fillOpacity=0, width=0.3, opacity=0.3),
    popup=True,
    categorical=True,
    m=_map,
)
display(_map)

# %% [markdown] user_expressions=[]
# ## Review Single Scene
#
# Let's have a closer look at one scene in Africa.

# %%
items_179_073 = [
    item
    for item in items
    if (item.properties["landsat:wrs_path"], item.properties["landsat:wrs_row"]) == ("179", "073")
]
items_179_073

# %% [markdown] user_expressions=[]
# ### Load without using Dask
#
# If `chunks=` parameter is not supplied data loading happens locally without using Dask. Concurrency is still suported though, we request to use 4 CPU cores with `pool=4`. Progress can be observed by passing `progress=tqdm` option.

# %%
from IPython.display import Image
from tqdm.auto import tqdm

yy = stac_load(
    items_179_073,
    ["red", "green", "blue", "nir08"],
    dtype="uint16",
    nodata=0,
    patch_url=planetary_computer.sign,
    pool=4,
    progress=tqdm,
)

# display raster GeoBox
# - Loaded in Native projection/resolution
yy.odc.geobox

# %% [markdown] user_expressions=[]
# ### Generate Preview
#
# - Convert to RGBA
# - Shrink image to have no more than 1600 pixels per side using `cubic` resampling
# - Compress with `jpeg` and display
#   - `jpeg` does not support transparency so we ask for missing pixels to be white
#   - Default is `png` that does support transparency and is lossless, but is more heavy in terms of bytes
#   - Use `webp` when size matters but transparency is also needed (some browsers might not support it though)

# %%
#Image(
#    data=yy.odc.to_rgba(vmin=6_000, vmax=30_000)
#    .isel(time=0)
#    .odc.reproject(yy.odc.geobox.zoom_to(1600), resampling="cubic")
#    .odc.compress("jpeg", 90, transparent=(255, 255, 255))
#)

# %% [markdown] user_expressions=[]
# ### Visualize NIR band
#
# Individual bands can be visualized in false colour with `.odc.colorize` method (supports matplotlib colormaps).

# %%
nir = yy.isel(time=0).nir08
Image(
    data=nir.odc.reproject(yy.odc.geobox.zoom_to(1600), resampling="cubic")
    .odc.colorize("magma")
    .odc.compress("webp", 80)
)

# %% [markdown] user_expressions=[]
# ### Save Images to GeoTIFF
#
# Cloud Optimized GeoTIFF (COG) is a popular and well supported format. We can easily save xarrays to COG format.

# %%
fname_rgb = f"{items_179_073[0].id}-rgba-cog.tif"
fname_nir = f"{items_179_073[0].id}-nir-cog.tif"

print(f"Saving {fname_rgb}")
yy.isel(time=0).odc.to_rgba(vmin=6_000, vmax=30_000).odc.write_cog(
    fname_rgb,
    overwrite=True,
    overview_resampling="average",
    # there are other options, but we'll leave them as defaults
)

print(f"Saving {fname_nir}")
nir.odc.write_cog(
    fname_nir,
    overwrite=True,
    overview_resampling="average",
    # change compression to ZSTD (default is DEFLATE)
    compress="zstd",
    zstd_level=9,
)
print("Done")

# %%
# !ls -lh {fname_rgb} {fname_nir}
# !gdalinfo {fname_rgb}
# !gdalinfo {fname_nir}

# %% [markdown] user_expressions=[]
# ---------------------------------

# %% [markdown] user_expressions=[]
# ## Dump Versions

# %%
import odc.geo
import odc.stac

print(f"odc.stac=={odc.stac.__version__}")
print(f"odc.geo=={odc.geo.__version__}")
