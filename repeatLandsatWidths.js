//// repeatLandsatWidths
// 2016 George H Allen, see license file
// georgehenryallen@gmail.com

// Masks water, clouds, & shadows on a collection of Landsat 5, 7, and 8 scenes.
// Measures river width using GRWL-derived cross sections
// exports width measurements (with date and location) as a CSV
 
////////////////////////////////////////////////////////////////////////
// Load GRWL cross sections:
// feature with 40 cross sections per seg 
var fc = ee.FeatureCollection('ft:16S4bhTENwMkDVbp3jNOmCUUgwiy-8MsTCGgvke9I');
//print('Cross Section FC feature properties:', fcProperties);

// buffer feature and edit attributes: 
var buffered = fc.map(function(f) {
  f.set({xSecLength: f.geometry().length().toLong()}); // add feature length attribute
  f.buffer(15); // buffer feature geometry
  return ee.Feature(f.select(["segmentID", "segmentInd", "xSecLength"], null, true)); // extract needed attributes
});

print("buffered feature collection:", buffered);

////////////////////////////////////////////////////////////////////////
// Load Landsat collections:
// fixme: consider pansharpening for 15 m resolution:

// Match common names to the sensor-specific bands:
var LT5_BANDS = ['B1',   'B2',    'B3',  'B4',  'B5',    'B7',    'B6'];
var LE7_BANDS = ['B1',   'B2',    'B3',  'B4',  'B5',    'B7',    'B6_VCID_1'];
var LC8_BANDS = ['B2',   'B3',    'B4',  'B5',  'B6',    'B7',    'B10'];
var STD_NAMES = ['blue', 'green', 'red', 'nir', 'swir1', 'swir2', 'temp'];

// load Landsat 5,7,8 collections:
// TODO(GHA): combine 5, 7, and 8 collections:
//var LT5 = ee.ImageCollection('LANDSAT/LT5_L1T_TOA')
//    .select(LT5_BANDS, STD_NAMES);
var LT5 = ee.ImageCollection('LANDSAT/LT5_L1T')
    .select(LT5_BANDS, STD_NAMES); 
// var LE7 = ee.ImageCollection('LANDSAT/LE7_L1T_TOA')
//   .select(LE7_BANDS, STD_NAMES);
// var LC8 = ee.ImageCollection('LANDSAT/LC8_L1T_TOA')
//   .select(LC8_BANDS, STD_NAMES);

var collection = LT5;//.merge(LE7).merge(LC8);

//var collection = ee.ImageCollection('LANDSAT/LC8_L1T_TOA');
var filtered = collection
    //.filterMetadata('WRS_PATH', 'equals', 34)
    //.filterMetadata('WRS_ROW', 'equals', 27)
    .filterDate('2012-11-01', '2012-12-31')
    .sort('system:time_start', false)
    .filterMetadata('CLOUD_COVER', 'less_than', 20)
    .filterBounds(fc);


////////////////////////////////////////////////////////////////////////
// classify filtered collection as water or cloud+shadow, & add cloud attribute:
var maskCollection = filtered.map(function(image){
  var cm = ee.Image(0); // cloud mask
  // water
  var mndwi = image.normalizedDifference(['green', 'nir']); 
  cm = cm.where(mndwi.gt(-0.14), 1);
  // shadows
  cm = cm.where(image.select('green').lt(0.07), 2); 
  // clouds:
  // A helper to apply an expression and linearly rescale the output.
  var rs = function(image, exp, thresholds) {
    return image.expression(exp, {image:image})
        .subtract(thresholds[0]).divide(thresholds[1]-thresholds[0]);
  };
  
  var s = ee.Image(1.0);
  s = s.min(rs(image, 'image.blue', [0.1, 0.3]));
  s = s.min(rs(image, 'image.red + image.green + image.blue', [0.2, 0.8])); 
  s = s.min(rs(image, 'image.nir + image.swir1 + image.swir2', [0.3, 0.8]));
  s = s.min(rs(image, 'image.temp', [300, 290])); 
  var ndsi = image.normalizedDifference(['green', 'swir1']);
  var cloudiness = s.min(rs(ndsi, 'image', [0.8, 0.6]));
  cm = cm.where(cloudiness.gt(0.15), 3);
  
  var wm = ee.Image(0); // water mask
  wm = wm.where(cm.eq(1), 1);
  cm = cm.where(cm.eq(1), 0);

  // http://stackoverflow.com/questions/2666112/convert-miliseconds-to-date-in-excel
  var date = ee.Number(image.get('system:time_start'));
  
  // Add bands to the wm
  wm = wm.addBands([date, cm])
            .rename(["water", "date", "cloud"]);
  
  return(wm);
  
});

////////////////////////////////////////////////////////////////////////
// take average of masks in buffered cross sections:
var maskTab = maskCollection.map(function(image){
  var reduced = image.reduceRegions({
                                      collection: buffered,
                                      reducer: ee.Reducer.mean(),
                                      scale: 30
                                    });
  return reduced;
});

maskTab = maskTab.flatten();

////////////////////////////////////////////////////////////////////////
// export table to drive:

// keep only essential properties in table: 
var maskTab = maskTab.map(function(f) {
                            return ee.Feature(f.select(
                              ["segmentID", "segmentInd", "date", 
                              'water', "cloud", "xSecLength"]
                              , null, false));
                          });


// TODO(GHA): multiply properties to get width property:
// var maskTab = maskTab.map(function(f) {
//                             return ee.Feature(f.select(['water']);
//                           });
// var water = ee.List(maskTab.aggregate_array('water'));
// var xSecLength = ee.List(maskTab.aggregate_array('xSecLength'));
// var width = water.toArray().divide(xSecLength.toArray());
// var width = water.zip(xSecLength).map(function(f) { return ee.Number(f.get(0)).divide(f.get(1))});

Export.table(maskTab, 'AMHG_LT5_1', 
            {fileFormat: 'CSV',
              driveFolder: 'GEE',
              driveFileNamePrefix: 'AMHG_LT5_1'
            });

////////////////////////////////////////////////////////////////////////
//Export mask video to drive:
// var vid = wmMn
//     .select(['clpd', 'nir', 'green'])
//     .map(function(image) {
//       return  image.multiply(512).uint8();
//     });



// Export.video(vid, 'sinTile_masks_vid01', 
//             {dimensions: '720',
//             framesPerSecond: 4,
//             region: ee.Geometry.Rectangle([-104.27, 48.07, -104.08, 47.98]).toGeoJSON(),
//             driveFolder: 'GEE',
//             driveFileNamePrefix: 'waterVideo_v01'
//             });


////////////////////////////////////////////////////////////////////////
// Display:
Map.addLayer(filtered.reduce(ee.Reducer.first()), {bands: ['nir_first', 'green_first', 'green_first'], 
                      gamma: 1.6}, 'first orig image');

Map.addLayer(maskCollection.select('water').reduce(ee.Reducer.first()), 
              {palette: ['f2f2f2', '0033cc'],
              opacity: 1}, 
              "first water mask"
              );

Map.addLayer(maskCollection.select('cloud').reduce(ee.Reducer.first()),
            {palette: ['000000', 'ff0000', 'ffffff'],
            opacity: 1}, 
            "first cloud mask"
            );

Map.addLayer(buffered.draw({
            color: 'ffbf00', 
            strokeWidth: 1, 
            pointRadius: 4}), 
            {},
            'GRWL buffer'
            );

Map.centerObject(fc.geometry(), 9);


