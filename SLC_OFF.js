// ┌─────────────────────────────────────────────────────────────────────┐
// │ Hybrid Gap-Filling for Landsat 7 SLC-off Images — FIXED + CLOUD FILTER │
// │ Method: Spatial NN (reduced kernel) + Temporal Mean + ≤20% Cloud Cover │
// └─────────────────────────────────────────────────────────────────────┘

var targetDate = '2020-07-15';
var dateRange = 365;
var roi = geometry;

// ─────────────────────────────────────────────────────────────────────
// 1. Load and preprocess —  Filter by cloud cover ≤ 20%
// ─────────────────────────────────────────────────────────────────────

var l7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2')
  .filterBounds(roi)
  .filterDate(
    ee.Date(targetDate).advance(-dateRange, 'day'),
    ee.Date(targetDate).advance(dateRange, 'day')
  )
  .filter(ee.Filter.lte('CLOUD_COVER', 20)); //  Only use scenes with ≤20% cloud cover

print(' Total images after cloud filtering:', l7.size());

function maskL7sr(image) {
  var qaMask = image.select('QA_PIXEL').bitwiseAnd(parseInt('11', 2)).eq(0); // clear or water
  var saturationMask = image.select('QA_RADSAT').eq(0); // no saturation
  var opticalBands = image.select('SR_B.').multiply(0.0000275).add(-0.2);
  var thermalBand = image.select('ST_B6').multiply(0.00341802).add(149.0);
  return image.addBands(opticalBands, null, true)
              .addBands(thermalBand, null, true)
              .updateMask(qaMask)
              .updateMask(saturationMask);
}

var l7masked = l7.map(maskL7sr);

// ─────────────────────────────────────────────────────────────────────
// 2. Select target image (closest to target date)
// ─────────────────────────────────────────────────────────────────────

var targetImage = l7masked
  .map(function(img) {
    var dateDiff = ee.Number(img.date().difference(ee.Date(targetDate), 'day')).abs();
    return img.set('dateDiff', dateDiff);
  })
  .sort('dateDiff')
  .first();

print('Target image selected (cloud_cover =', targetImage.get('CLOUD_COVER'), '%)');

// ─────────────────────────────────────────────────────────────────────
// 3. Build temporal mean composite (exclude target image)
// ─────────────────────────────────────────────────────────────────────

var referenceImages = l7masked.filter(ee.Filter.neq('system:index', targetImage.get('system:index')));
print('Reference images for temporal composite:', referenceImages.size());

var meanComposite = referenceImages.mean();

// ─────────────────────────────────────────────────────────────────────
// 4. Spatial gap filling — reduced kernel + clip
// ─────────────────────────────────────────────────────────────────────

var bandNames = targetImage.bandNames();

function fillBand(image, bandName) {
  var band = image.select(bandName);
  var filled = band
    .clip(roi)
    .focal_max({radius: 30, units: 'pixels', iterations: 10}); // Reduced kernel
  return band.unmask(filled).rename(bandName);
}

var firstBandName = ee.String(bandNames.get(0));
var filledImage = fillBand(targetImage, firstBandName);

var filledImageFinal = ee.Image(
  bandNames.slice(1).iterate(function(bandName, img) {
    img = ee.Image(img);
    var filledBand = fillBand(targetImage, ee.String(bandName));
    return img.addBands(filledBand);
  }, filledImage)
);

filledImageFinal = ee.Image(filledImageFinal).select(bandNames);

// ─────────────────────────────────────────────────────────────────────
// 5. Final fill with temporal mean
// ─────────────────────────────────────────────────────────────────────

var finalImage = filledImageFinal.unmask(meanComposite)
  .copyProperties(targetImage, targetImage.propertyNames());

finalImage = ee.Image(finalImage); // explicit cast

// ─────────────────────────────────────────────────────────────────────
// 6. Visualization
// ─────────────────────────────────────────────────────────────────────

var visParams = {
  bands: ['SR_B3', 'SR_B2', 'SR_B1'],
  min: 0.0,
  max: 0.3,
  gamma: 1.2
};

Map.centerObject(roi, 10);
Map.addLayer(targetImage, visParams, '1. Original (with gaps)', false);
Map.addLayer(filledImageFinal, visParams, '2. After Spatial NN', false);
Map.addLayer(finalImage, visParams, '3. Final: Spatial + Temporal', true);

// ─────────────────────────────────────────────────────────────────────
// 7. Export
// ─────────────────────────────────────────────────────────────────────

Export.image.toDrive({
  image: finalImage.clip(roi),
  description: 'L7_SLCoff_CloudFiltered_GapFilled',
  scale: 30,
  region: roi,
  maxPixels: 1e10,
  fileFormat: 'GeoTIFF',
  formatOptions: { cloudOptimized: true }
});