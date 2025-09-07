// ┌─────────────────────────────────────────────────────────────────────┐
// │ Generate Multi-Temporal Seasonal Composite Using Median             │
// │ To reduce noise from seasonal lighting/vegetation changes           │
// │ Ensures consistent spectral info for classification                 │
// └─────────────────────────────────────────────────────────────────────┘

// Define area and season
var roi = table2; //  study area (import from assets folder : clip for sundarban and modhupur)


var year = 2022; // Adjust as needed
var startDate = ee.Date.fromYMD(year, 6, 1);   // June 1
var endDate = ee.Date.fromYMD(year, 8, 31);    // August 31

print('Season:', startDate.format('YYYY-MM-dd'), 'to', endDate.format('YYYY-MM-dd'));

//  Load Landsat 7 Surface Reflectance (C02)
var l7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
  .filterBounds(roi)
  .filterDate(startDate, endDate)
  .filter(ee.Filter.lte('CLOUD_COVER', 20)); // Optional: reduce cloud impact

print(' Images in season:', l7.size());

// Cloud masking function (Landsat 7 Level 2)
function maskL7sr(image) {
  var qaMask = image.select('QA_PIXEL').bitwiseAnd(parseInt('11', 2)).eq(0); // clear or water
  var saturationMask = image.select('QA_RADSAT').eq(0); // no saturation
  // Scale optical bands
  var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
  // Scale thermal band
  var thermalBand = image.select('ST_B6').multiply(0.00341802).add(149.0);
  return image
    .addBands(opticalBands, null, true)
    .addBands(thermalBand, null, true)
    .updateMask(qaMask)
    .updateMask(saturationMask);
}

//  Apply mask
var l7masked = l7.map(maskL7sr);

// 
// "Generated multi-temporal composites by combining images from different dates within the same season, using median values"

var seasonalMedianComposite = l7masked.median(); // Median reduces noise, smooths short-term variations

print('seasonal median composite generated');

// Visualization (Natural Color)
var visParams = {
  bands: ['SR_B3', 'SR_B2', 'SR_B1'], // Red, Green, Blue
  min: 0.0,
  max: 0.3,
  gamma: 1.2
};

Map.centerObject(roi, 10);
Map.addLayer(seasonalMedianComposite, visParams, ' Seasonal Median Composite (Reduced Seasonal Noise)', true);

//  Optional: Export for classification
Export.image.toDrive({
  image: seasonalMedianComposite,
  description: 'L7_Seasonal_Median_Composite_' + year + '_Summer',
  scale: 30,
  region: roi,
  maxPixels: 1e10,
  fileFormat: 'GeoTIFF',
  formatOptions: {
    cloudOptimized: true
  }
});