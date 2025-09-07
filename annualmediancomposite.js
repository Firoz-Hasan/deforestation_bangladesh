// ┌─────────────────────────────────────────────────────────────────────┐
// │ Cloud-Free Annual Median Composite for Tropical Forest Monitoring   │
// │ Specifically designed for regions like Sundarbans with monsoon clouds│
// └─────────────────────────────────────────────────────────────────────┘

// Define region of interest — Example: Sundarbans (adjust as needed)
//var sundarbans = ee.Geometry.Rectangle([88.5, 21.5, 90.0, 22.7]); // Approximate Sundarbans bounds
var roi = table; //  use 'geometry'

// Define year of interest
var year = 2020; // Change as needed
var startDate = ee.Date.fromYMD(year, 1, 1);
var endDate = ee.Date.fromYMD(year, 12, 31);

print('Processing year:', year);

// Load Landsat 7 Surface Reflectance (C02) — suitable for forest monitoring
var l7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
  .filterBounds(roi)
  .filterDate(startDate, endDate);


// " Google Earth Engine's (GEE) CLOUD_COVER metadata to filter out scenes with cloud cover exceeding less than 20%."

var l7_filtered = l7.filter(ee.Filter.lte('CLOUD_COVER', 20));

print('Scenes after filtering (≤20% cloud cover):', l7_filtered.size());

//Apply pixel-level cloud and saturation masking
function maskL7sr(image) {
  var qaMask = image.select('QA_PIXEL').bitwiseAnd(parseInt('11', 2)).eq(0); // clear or water
  var saturationMask = image.select('QA_RADSAT').eq(0); // no band saturation
  // Scale optical surface reflectance bands
  var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
  // Scale surface temperature band
  var thermalBand = image.select('ST_B6').multiply(0.00341802).add(149.0);
  return image
    .addBands(opticalBands, null, true)
    .addBands(thermalBand, null, true)
    .updateMask(qaMask)
    .updateMask(saturationMask);
}

var l7_masked = l7_filtered.map(maskL7sr);


// "Additionally, we generated a temporal median composite from multiple cloud-free observations taken within the same year."

var annualMedianComposite = l7_masked.median(); // Temporal median from multiple cloud-free obs in same year

print('annual median composite generated');


// "This approach effectively suppressed residual cloud artefacts and provided a cleaner, cloud-free image that preserved the underlying spectral information of the forest."

//  Visualization (Natural Color for Forest)
var visParams = {
  bands: ['SR_B3', 'SR_B2', 'SR_B1'], // Red, Green, Blue
  min: 0.0,
  max: 0.3,
  gamma: 1.3 // Slight contrast boost for forest structure
};

Map.centerObject(roi, 10);
Map.addLayer(annualMedianComposite, visParams, 'Annual Median Composite (Cloud-Free Forest)', true);

//Export for further analysis or classification
Export.image.toDrive({
  image: annualMedianComposite,
  description: 'L7_Annual_Median_Composite_' + year + '_Sundarbans',
  scale: 30,
  region: roi,
  maxPixels: 1e10,
  fileFormat: 'GeoTIFF',
  formatOptions: {
    cloudOptimized: true
  }
});