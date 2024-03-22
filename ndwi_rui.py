import ee
import geopandas as gpd

# Initialize Earth Engine
ee.Initialize()

# AOI
# Read AOI shapefile --------
aoi_file = gpd.read_file("AOI/luanda_catchment_level4.shp").to_crs(epsg = 4326)

# Convert shapefile to ee.Geometry ------------
jsonDict = eval(aoi_file['geometry'].to_json())

if len(jsonDict['features']) > 1:
    print('Need to convert polygons into a multipolygon')
    print('or do something else, like creating individual raster for each polygon and then merge')
    exit()

jsonDict['features'][0]['geometry']['coordinates'][0] = [x[:-1] for x in jsonDict['features'][0]['geometry']['coordinates'][0]]
AOI = ee.Geometry.MultiPolygon(jsonDict['features'][0]['geometry']['coordinates'])

id = "luanda"

# Create Water Mask
waterMask = ee.Image("JRC/GSW1_0/GlobalSurfaceWater").select('transition')
blank = ee.Image(0)
nonWater = blank.addBands(waterMask).unmask().select('transition').eq(0).rename('non_water')

# Create Slope Mask
srtm = ee.Image("USGS/SRTMGL1_003")
slope = ee.Terrain.slope(srtm)
Grade15 = slope.gt(15)
gtGrade15 = Grade15.updateMask(Grade15.neq(0))
slopeMask = gtGrade15.clip(AOI).unmask(0).subtract(1).multiply(-1)

# Load Sentinel-1 C-band SAR Ground Range collection (log scaling, VV co-polar)
collection = ee.ImageCollection('COPERNICUS/S1_GRD') \
    .filterBounds(AOI) \
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV')) \
    .filter(ee.Filter.eq('orbitProperties_pass', 'DESCENDING')) \
    .select('VV')

# Filter by date
# Flood date: 04/02/2023 to 04/19/2023
during_flood_collection = collection.filterDate('2023-04-03', '2023-04-20')
# Iterate over the during_flood_collection and print image dates
image_dates = during_flood_collection.aggregate_array('system:time_start').getInfo()

print("Dates of images in during_flood_collection:")
for date in image_dates:
    print(ee.Date(date).format('YYYY-MM-dd').getInfo())
exit()

during_flood = during_flood_collection.first().clip(AOI).updateMask(nonWater).updateMask(slopeMask)

# Before Flood & in Dry Season: May - Sept
before_flood_collection = collection.filterDate('2022-08-01', '2022-09-30')
before_flood = before_flood_collection.median().clip(AOI).updateMask(nonWater).updateMask(slopeMask)

# Set smoothing box (30m is 3x size of S1 pixel)
smoothing_box = 30

# Create difference images (during flood - before flood)
difference_smoothed = during_flood \
    .focal_median(smoothing_box, 'square', 'meters') \
    .subtract(before_flood.focal_median(smoothing_box, 'square', 'meters'))
difference_raw = during_flood.subtract(before_flood)
print(difference_smoothed.projection().getInfo())
exit()
# Resample difference raster
difference_smoothed_resampled = difference_smoothed.resample('bilinear').reproject(
    crs=difference_smoothed.projection().getInfo()['crs'],
    scale=10.0
)

# Calculate Lower 10% value for setting upper (Darker) threshold
UpperThreshold = difference_smoothed_resampled.reduceRegion(
    reducer=ee.Reducer.percentile([6]),
    geometry=AOI,
    scale=10,
    bestEffort=True
).get('VV')

# Calculate Top 10% value for setting lower (Brighter) threshold
LowerThreshold = difference_smoothed_resampled.reduceRegion(
    reducer=ee.Reducer.percentile([94]),
    geometry=AOI,
    scale=10,
    bestEffort=True
).get('VV')

print(UpperThreshold.getInfo())
print(LowerThreshold.getInfo())

# Make Difference Darker Rasters
difference_darker_smoothed_thresholded = difference_smoothed_resampled.lt(UpperThreshold)
difference_darker_thresholded = difference_raw.lt(UpperThreshold)

# Make Difference Brighter Rasters
difference_brighter_smoothed_thresholded = difference_smoothed_resampled.gt(LowerThreshold)
difference_brighter_thresholded = difference_raw.gt(LowerThreshold)

# DARKER
# Compute connected pixel counts; stop searching for connected pixels
# once the size of the connected neighborhood reaches 30 pixels, and
# use 8-connected rules.
conn_D = difference_darker_smoothed_thresholded.connectedPixelCount(400, True)
smallClusters_D = conn_D.lt(200)
onesSurroundedByZeros_D = difference_darker_smoothed_thresholded.And(smallClusters_D)
cleaned1Islands_D = difference_darker_smoothed_thresholded.where(onesSurroundedByZeros_D, 0)
zerosSurroundedByOnes_D = smallClusters_D.And(difference_darker_smoothed_thresholded.Not())
cleaned1and0Islands_D = cleaned1Islands_D.where(zerosSurroundedByOnes_D, 1)
floodedAreas_filtered_D = cleaned1and0Islands_D.mask(cleaned1and0Islands_D)

# BRIGHTER
# Compute connected pixel counts; stop searching for connected pixels
# once the size of the connected neighborhood reaches 30 pixels, and
# use 8-connected rules.
conn_B = difference_brighter_smoothed_thresholded.connectedPixelCount(400, True)
smallClusters_B = conn_B.lt(200)
onesSurroundedByZeros_B = difference_brighter_smoothed_thresholded.And(smallClusters_B)
cleaned1Islands_B = difference_brighter_smoothed_thresholded.where(onesSurroundedByZeros_B, 0)
zerosSurroundedByOnes_B = smallClusters_B.And(difference_brighter_smoothed_thresholded.Not())
cleaned1and0Islands_B = cleaned1Islands_B.where(zerosSurroundedByOnes_B, 1)
floodedAreas_filtered_B = cleaned1and0Islands_B.mask(cleaned1and0Islands_B)


# Save Resampled & Smoothed Difference Raster to Cloud Storage
task_smoothed_resampled = ee.batch.Export.image.toCloudStorage(
    image=difference_smoothed_resampled,
    description=id + "_Difference_Smoothed_S1",
    bucket='gee-test-20231019',
    region=AOI,
    scale=10,
    maxPixels=1e9
)
task_smoothed_resampled.start()

# Save Flood Extent (Darker) to Cloud Storage
task_floodedAreas_filtered_D = ee.batch.Export.image.toCloudStorage(
    image=floodedAreas_filtered_D,
    description=id + "_FloodExtent_Darker_S1",
    bucket='gee-test-20231019',
    region=AOI,
    scale=10,
    maxPixels=1e9
)
task_floodedAreas_filtered_D.start()

# Save Flood Extent (Brighter) to Cloud Storage
task_floodedAreas_filtered_B = ee.batch.Export.image.toCloudStorage(
    image=floodedAreas_filtered_B,
    description=id + "_FloodExtent_Brighter_S1",
    bucket='gee-test-20231019',
    region=AOI,
    scale=10,
    maxPixels=1e9
)
task_floodedAreas_filtered_B.start()

# Save Dry Season Scene to Cloud Storage
task_before_flood = ee.batch.Export.image.toCloudStorage(
    image=before_flood,
    description=id + "_DrySeason_S1",
    bucket='gee-test-20231019',
    region=AOI,
    scale=10,
    maxPixels=1e9
)
task_before_flood.start()

# Save Wet Scene to Cloud Storage
task_during_flood = ee.batch.Export.image.toCloudStorage(
    image=during_flood,
    description=id + "_DuringFlood_S1",
    bucket='gee-test-20231019',
    region=AOI,
    scale=10,
    maxPixels=1e9
)
task_during_flood.start()
