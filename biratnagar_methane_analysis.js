// =====================================================================
// BIRATNAGAR ATMOSPHERIC METHANE STUDY
// Sentinel-5P XCH4 and environmental drivers, 2019–2025

// =====================================================================
// 1. CONFIGURATION
// =====================================================================

var START_DATE = ee.Date('2019-03-01');
var END_DATE   = ee.Date('2026-01-01');

var START_YEAR = 2019;
var END_YEAR   = 2025;

var CH4_BAND =
  'CH4_column_volume_mixing_ratio_dry_air_bias_corrected';

// Approximate effective analysis scale for gridded methane data.
var CH4_SCALE = 7000;

// Exploratory Sentinel-1 water threshold.
// Sensitivity testing should also use -15, -19 and -21 dB.
var WATER_THRESHOLD_DB = -17;

// Minimum requirements for valid monthly methane data.
var MIN_CH4_IMAGES = 5;
var MIN_COVERAGE_PERCENT = 25;


// =====================================================================
// 2. STUDY AREA
// =====================================================================

var brt = ee.FeatureCollection(
  'projects/graphite-bliss-486713-v5/assets/Biratnagar_submetro'
);

var brtGeometry = brt.geometry();


// Urban core
var urbanCore = brtGeometry;


// Peri-urban zone: 0–10 km outside Biratnagar
var buffer10 = brtGeometry.buffer(10000);

var periUrbanRing = buffer10.difference(
  brtGeometry,
  ee.ErrorMargin(100)
);


// Regional-background zone: 25–50 km
var buffer25 = brtGeometry.buffer(25000);
var buffer50 = brtGeometry.buffer(50000);

var regionalBackground = buffer50.difference(
  buffer25,
  ee.ErrorMargin(100)
);


// Overall analytical domain
var studyDomain = buffer50;


// Analytical zones
var zones = ee.FeatureCollection([

  ee.Feature(urbanCore, {
    Zone: 'Urban_Core'
  }),

  ee.Feature(periUrbanRing, {
    Zone: 'Peri_Urban_0_10km'
  }),

  ee.Feature(regionalBackground, {
    Zone: 'Regional_Background_25_50km'
  })
]);


Map.centerObject(brt, 9);

Map.addLayer(
  brt,
  {color: 'red'},
  'Biratnagar boundary'
);

Map.addLayer(
  periUrbanRing,
  {color: 'yellow'},
  'Peri-urban ring: 0–10 km',
  false
);

Map.addLayer(
  regionalBackground,
  {color: 'blue'},
  'Regional background: 25–50 km',
  false
);


// =====================================================================
// 3. SENTINEL-5P METHANE
// =====================================================================

// Do not map the mask over the complete 2019–2025 collection.
var methaneRaw = ee.ImageCollection(
  'COPERNICUS/S5P/OFFL/L3_CH4'
)
.filterBounds(studyDomain)
.select(CH4_BAND);


// Physical-range screening
function maskMethane(image) {

  var ch4 = image.select(CH4_BAND);

  var valid = ch4.gt(1500)
                 .and(ch4.lt(2500));

  return image
    .updateMask(valid)
    .copyProperties(
      image,
      ['system:time_start']
    );
}


// Date-filter first and then apply masking.
function getMethaneCollection(start, end) {

  return methaneRaw
    .filterDate(start, end)
    .map(maskMethane);
}


// Return a monthly mean methane image.
function getMonthlyMethaneMean(start) {

  start = ee.Date(start);

  var end = start.advance(1, 'month');

  return getMethaneCollection(start, end)
    .mean()
    .select(CH4_BAND)
    .set({
      'system:time_start': start.millis(),
      Date: start.format('YYYY-MM')
    });
}


// Print only a limited example.
print(
  'Example methane image',
  getMethaneCollection(
    '2025-01-01',
    '2025-02-01'
  ).first()
);

print(
  'January 2025 source-image count',
  getMethaneCollection(
    '2025-01-01',
    '2025-02-01'
  ).size()
);


// =====================================================================
// 4. ENVIRONMENTAL DATASETS
// =====================================================================

// Daily CHIRPS rainfall
var chirps = ee.ImageCollection(
  'UCSB-CHG/CHIRPS/DAILY'
)
.filterDate(START_DATE, END_DATE)
.filterBounds(studyDomain)
.select('precipitation');


// ERA5-Land monthly data
var era5 = ee.ImageCollection(
  'ECMWF/ERA5_LAND/MONTHLY_AGGR'
)
.filterDate(START_DATE, END_DATE)
.filterBounds(studyDomain);


// MODIS NDVI
var modisNDVI = ee.ImageCollection(
  'MODIS/061/MOD13Q1'
)
.filterDate(START_DATE, END_DATE)
.filterBounds(studyDomain)
.select('NDVI');


// Sentinel-1 VV data
var sentinel1 = ee.ImageCollection(
  'COPERNICUS/S1_GRD'
)
.filterDate(START_DATE, END_DATE)
.filterBounds(studyDomain)
.filter(
  ee.Filter.eq(
    'instrumentMode',
    'IW'
  )
)
.filter(
  ee.Filter.listContains(
    'transmitterReceiverPolarisation',
    'VV'
  )
)
.select('VV');


// ESA WorldCover
var worldCover = ee.ImageCollection(
  'ESA/WorldCover/v200'
)
.first()
.select('Map');


// Cropland class
var croplandMask = worldCover
  .eq(40)
  .rename('cropland');


// Built-up class
var builtMask = worldCover
  .eq(50)
  .rename('built_up');


// Permanent water
var permanentWater = ee.Image(
  'JRC/GSW1_4/GlobalSurfaceWater'
)
.select('occurrence')
.gte(90);


// =====================================================================
// 5. REDUCTION FUNCTIONS
// =====================================================================

function regionMean(image, band, geometry, scale) {

  return image.select(band).reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: geometry,
    scale: scale,
    bestEffort: true,
    maxPixels: 1e13,
    tileScale: 4
  }).get(band);
}


function regionMedian(image, band, geometry, scale) {

  return image.select(band).reduceRegion({
    reducer: ee.Reducer.median(),
    geometry: geometry,
    scale: scale,
    bestEffort: true,
    maxPixels: 1e13,
    tileScale: 4
  }).get(band);
}


function regionSD(image, band, geometry, scale) {

  return image.select(band).reduceRegion({
    reducer: ee.Reducer.stdDev(),
    geometry: geometry,
    scale: scale,
    bestEffort: true,
    maxPixels: 1e13,
    tileScale: 4
  }).get(band);
}


// Calculate percentage of geometry represented by a binary mask.
function areaPercentage(maskImage, geometry, scale) {

  var pixelArea = ee.Image.pixelArea();

  var selectedAreaDictionary = pixelArea
    .updateMask(maskImage)
    .reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: geometry,
      scale: scale,
      bestEffort: true,
      maxPixels: 1e13,
      tileScale: 4
    });

  var totalAreaDictionary = pixelArea.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: geometry,
    scale: scale,
    bestEffort: true,
    maxPixels: 1e13,
    tileScale: 4
  });

  var selectedArea = ee.Number(
    ee.Algorithms.If(
      selectedAreaDictionary.contains('area'),
      selectedAreaDictionary.get('area'),
      0
    )
  );

  var totalArea = ee.Number(
    totalAreaDictionary.get('area')
  );

  return selectedArea
    .divide(totalArea)
    .multiply(100);
}


// Safely subtract a background value.
function calculateEnhancement(localValue, backgroundValue) {

  return ee.Algorithms.If(
    ee.Algorithms.IsEqual(
      localValue,
      null
    ),
    null,

    ee.Algorithms.If(
      ee.Algorithms.IsEqual(
        backgroundValue,
        null
      ),
      null,

      ee.Number(localValue).subtract(
        ee.Number(backgroundValue)
      )
    )
  );
}


// =====================================================================
// 6. STATIC LAND-COVER STATISTICS
// =====================================================================

var urbanCroplandPercent = areaPercentage(
  croplandMask,
  urbanCore,
  10
);

var periUrbanCroplandPercent = areaPercentage(
  croplandMask,
  periUrbanRing,
  10
);

var backgroundCroplandPercent = areaPercentage(
  croplandMask,
  regionalBackground,
  10
);

var urbanBuiltPercent = areaPercentage(
  builtMask,
  urbanCore,
  10
);

var periUrbanBuiltPercent = areaPercentage(
  builtMask,
  periUrbanRing,
  10
);

var backgroundBuiltPercent = areaPercentage(
  builtMask,
  regionalBackground,
  10
);


print(
  'Urban cropland percentage',
  urbanCroplandPercent
);

print(
  'Peri-urban cropland percentage',
  periUrbanCroplandPercent
);

print(
  'Urban built-up percentage',
  urbanBuiltPercent
);


// =====================================================================
// 7. MONTHLY DATE SEQUENCE
// =====================================================================

var numberOfMonths = END_DATE.difference(
  START_DATE,
  'month'
);

var monthlyDates = ee.List.sequence(
  0,
  numberOfMonths.subtract(1)
).map(function(number) {

  return START_DATE.advance(
    number,
    'month'
  );
});


// =====================================================================
// 8. MONTHLY ANALYSIS TABLE
// =====================================================================

var monthlyAnalysis = ee.FeatureCollection(

  monthlyDates.map(function(dateValue) {

    var start = ee.Date(dateValue);
    var end = start.advance(1, 'month');

    var year = ee.Number(
      start.get('year')
    );

    var month = ee.Number(
      start.get('month')
    );


    // ---------------------------------------------------------------
    // 8.1 Monthly methane
    // ---------------------------------------------------------------

    var monthlyCH4Collection =
      getMethaneCollection(start, end);

    var monthlyCH4Mean =
      monthlyCH4Collection
        .mean()
        .rename('XCH4');

    var monthlyCH4Median =
      monthlyCH4Collection
        .median()
        .rename('XCH4_median');

    var monthlyTemporalSD =
      monthlyCH4Collection
        .reduce(ee.Reducer.stdDev())
        .rename('XCH4_temporal_sd');

    var observationCount =
      monthlyCH4Collection
        .count()
        .rename('valid_count');


    var meanValidObservations = regionMean(
      observationCount,
      'valid_count',
      urbanCore,
      CH4_SCALE
    );


    var validCoveragePercent = areaPercentage(
      observationCount.gt(0),
      urbanCore,
      CH4_SCALE
    );


    var urbanCH4 = regionMean(
      monthlyCH4Mean,
      'XCH4',
      urbanCore,
      CH4_SCALE
    );

    var periurbanCH4 = regionMean(
      monthlyCH4Mean,
      'XCH4',
      periUrbanRing,
      CH4_SCALE
    );

    var backgroundCH4 = regionMean(
      monthlyCH4Mean,
      'XCH4',
      regionalBackground,
      CH4_SCALE
    );

    var urbanMedianCH4 = regionMedian(
      monthlyCH4Median,
      'XCH4_median',
      urbanCore,
      CH4_SCALE
    );

    var urbanTemporalSD = regionMean(
      monthlyTemporalSD,
      'XCH4_temporal_sd',
      urbanCore,
      CH4_SCALE
    );


    var urbanEnhancement =
      calculateEnhancement(
        urbanCH4,
        backgroundCH4
      );

    var periurbanEnhancement =
      calculateEnhancement(
        periurbanCH4,
        backgroundCH4
      );


    // ---------------------------------------------------------------
    // 8.2 Monthly rainfall
    // ---------------------------------------------------------------

    var monthlyRainfall = chirps
      .filterDate(start, end)
      .sum()
      .rename('rainfall_mm');

    var rainfallMM = regionMean(
      monthlyRainfall,
      'rainfall_mm',
      urbanCore,
      5566
    );


    // ---------------------------------------------------------------
    // 8.3 ERA5-Land drivers
    // ---------------------------------------------------------------

    var monthlyERA5 = ee.Image(
      era5.filterDate(start, end).first()
    );


    var temperatureC = monthlyERA5
      .select('temperature_2m')
      .subtract(273.15)
      .rename('temperature_c');


    var soilMoisture = monthlyERA5
      .select(
        'volumetric_soil_water_layer_1'
      )
      .rename('soil_moisture');


    var windU = monthlyERA5
      .select(
        'u_component_of_wind_10m'
      );


    var windV = monthlyERA5
      .select(
        'v_component_of_wind_10m'
      );


    var windSpeed = windU
      .pow(2)
      .add(windV.pow(2))
      .sqrt()
      .rename('wind_speed');


    var windDirection = windU
      .atan2(windV)
      .multiply(180 / Math.PI)
      .add(180)
      .mod(360)
      .rename('wind_direction');


    var meanTemperatureC = regionMean(
      temperatureC,
      'temperature_c',
      urbanCore,
      11132
    );

    var meanSoilMoisture = regionMean(
      soilMoisture,
      'soil_moisture',
      urbanCore,
      11132
    );

    var meanWindSpeed = regionMean(
      windSpeed,
      'wind_speed',
      urbanCore,
      11132
    );

    var meanWindDirection = regionMean(
      windDirection,
      'wind_direction',
      urbanCore,
      11132
    );


    // ---------------------------------------------------------------
    // 8.4 Monthly NDVI
    // ---------------------------------------------------------------

    var monthlyNDVI = modisNDVI
      .filterDate(start, end)
      .mean()
      .multiply(0.0001)
      .rename('NDVI');


    var meanNDVI = regionMean(
      monthlyNDVI,
      'NDVI',
      urbanCore,
      250
    );


    // ---------------------------------------------------------------
    // 8.5 Monthly inundation
    // ---------------------------------------------------------------

    var monthlyS1 = sentinel1.filterDate(
      start,
      end
    );


    var monthlyWaterCollection = monthlyS1.map(
      function(image) {

        var temporaryWater = image
          .lt(WATER_THRESHOLD_DB)
          .and(permanentWater.not())
          .rename('temporary_water');

        return temporaryWater.copyProperties(
          image,
          ['system:time_start']
        );
      }
    );


    var inundationFrequency =
      monthlyWaterCollection
        .mean()
        .rename('inundation_frequency');


    var inundationMask =
      inundationFrequency.gte(0.5);


    var urbanInundatedPercent = areaPercentage(
      inundationMask,
      urbanCore,
      30
    );

    var periurbanInundatedPercent = areaPercentage(
      inundationMask,
      periUrbanRing,
      30
    );


    // ---------------------------------------------------------------
    // 8.6 Nepal season
    // ---------------------------------------------------------------

    var season = ee.String(
      ee.Algorithms.If(
        month.eq(12).or(month.lte(2)),
        'Winter',

        ee.Algorithms.If(
          month.gte(3).and(month.lte(5)),
          'Pre-Monsoon',

          ee.Algorithms.If(
            month.gte(6).and(month.lte(9)),
            'Monsoon',
            'Post-Monsoon'
          )
        )
      )
    );


    // ---------------------------------------------------------------
    // 8.7 Return monthly record
    // ---------------------------------------------------------------

    return ee.Feature(null, {

      Date:
        start.format('YYYY-MM'),

      Year:
        year,

      Month:
        month,

      Season:
        season,

      CH4_image_count:
        monthlyCH4Collection.size(),

      Mean_valid_observations:
        meanValidObservations,

      Valid_coverage_percent:
        validCoveragePercent,

      Urban_XCH4_mean_ppb:
        urbanCH4,

      Urban_XCH4_median_ppb:
        urbanMedianCH4,

      Urban_XCH4_temporal_SD_ppb:
        urbanTemporalSD,

      Periurban_XCH4_mean_ppb:
        periurbanCH4,

      Background_XCH4_mean_ppb:
        backgroundCH4,

      Urban_CH4_enhancement_ppb:
        urbanEnhancement,

      Periurban_CH4_enhancement_ppb:
        periurbanEnhancement,

      Rainfall_mm:
        rainfallMM,

      Temperature_C:
        meanTemperatureC,

      Soil_moisture_fraction:
        meanSoilMoisture,

      Wind_speed_ms:
        meanWindSpeed,

      Wind_direction_deg:
        meanWindDirection,

      NDVI:
        meanNDVI,

      Urban_inundated_percent:
        urbanInundatedPercent,

      Periurban_inundated_percent:
        periurbanInundatedPercent,

      Urban_cropland_percent:
        urbanCroplandPercent,

      Periurban_cropland_percent:
        periUrbanCroplandPercent,

      Background_cropland_percent:
        backgroundCroplandPercent,

      Urban_built_percent:
        urbanBuiltPercent,

      Periurban_built_percent:
        periUrbanBuiltPercent,

      Background_built_percent:
        backgroundBuiltPercent,

      S1_image_count:
        monthlyS1.size()
    });
  })
);


print(
  'Monthly methane and driver table',
  monthlyAnalysis
);


// =====================================================================
// 9. VALID MONTHLY OBSERVATIONS
// =====================================================================

var monthlyValid = monthlyAnalysis

  .filter(
    ee.Filter.gte(
      'CH4_image_count',
      MIN_CH4_IMAGES
    )
  )

  .filter(
    ee.Filter.gte(
      'Valid_coverage_percent',
      MIN_COVERAGE_PERCENT
    )
  )

  .filter(
    ee.Filter.notNull([
      'Urban_XCH4_mean_ppb',
      'Background_XCH4_mean_ppb',
      'Urban_CH4_enhancement_ppb'
    ])
  );


print(
  'Months passing coverage criteria',
  monthlyValid
);

print(
  'Number of valid months',
  monthlyValid.size()
);


// =====================================================================
// 10. MONTHLY EXPORTS
// =====================================================================

Export.table.toDrive({

  collection:
    monthlyAnalysis,

  description:
    'Biratnagar_Monthly_CH4_Drivers_All_2019_2025',

  fileNamePrefix:
    'Biratnagar_Monthly_CH4_Drivers_All_2019_2025',

  fileFormat:
    'CSV'
});


Export.table.toDrive({

  collection:
    monthlyValid,

  description:
    'Biratnagar_Monthly_CH4_Drivers_Valid_2019_2025',

  fileNamePrefix:
    'Biratnagar_Monthly_CH4_Drivers_Valid_2019_2025',

  fileFormat:
    'CSV'
});


// =====================================================================
// 11. SEASONAL ANALYSIS
// Winter 2019 is removed.
// =====================================================================

var years = ee.List.sequence(
  START_YEAR,
  END_YEAR
);


var seasonalAnalysis = ee.FeatureCollection(

  years.map(function(yearValue) {

    var year = ee.Number(yearValue);


    var seasonDefinitions = ee.List([

      ee.Dictionary({
        Season: 'Winter',
        Season_Order: 1,

        Start: ee.Date.fromYMD(
          year.subtract(1),
          12,
          1
        ),

        End: ee.Date.fromYMD(
          year,
          3,
          1
        )
      }),


      ee.Dictionary({
        Season: 'Pre-Monsoon',
        Season_Order: 2,

        Start: ee.Date.fromYMD(
          year,
          3,
          1
        ),

        End: ee.Date.fromYMD(
          year,
          6,
          1
        )
      }),


      ee.Dictionary({
        Season: 'Monsoon',
        Season_Order: 3,

        Start: ee.Date.fromYMD(
          year,
          6,
          1
        ),

        End: ee.Date.fromYMD(
          year,
          10,
          1
        )
      }),


      ee.Dictionary({
        Season: 'Post-Monsoon',
        Season_Order: 4,

        Start: ee.Date.fromYMD(
          year,
          10,
          1
        ),

        End: ee.Date.fromYMD(
          year,
          12,
          1
        )
      })
    ]);


    return seasonDefinitions.map(
      function(seasonValue) {

        var season = ee.Dictionary(
          seasonValue
        );

        var seasonName = ee.String(
          season.get('Season')
        );

        var start = ee.Date(
          season.get('Start')
        );

        var end = ee.Date(
          season.get('End')
        );


        // Maximum period is four months, preventing collection overflow.
        var seasonalCollection =
          getMethaneCollection(start, end);


        var seasonalMean = seasonalCollection
          .mean()
          .rename('XCH4');


        var seasonalTemporalSD = seasonalCollection
          .reduce(ee.Reducer.stdDev())
          .rename('XCH4_temporal_sd');


        var seasonalCount = seasonalCollection
          .count()
          .rename('valid_count');


        var coveragePercent = areaPercentage(
          seasonalCount.gt(0),
          urbanCore,
          CH4_SCALE
        );


        var urbanMean = regionMean(
          seasonalMean,
          'XCH4',
          urbanCore,
          CH4_SCALE
        );

        var periurbanMean = regionMean(
          seasonalMean,
          'XCH4',
          periUrbanRing,
          CH4_SCALE
        );

        var backgroundMean = regionMean(
          seasonalMean,
          'XCH4',
          regionalBackground,
          CH4_SCALE
        );

        var urbanSD = regionMean(
          seasonalTemporalSD,
          'XCH4_temporal_sd',
          urbanCore,
          CH4_SCALE
        );


        var urbanEnhancement =
          calculateEnhancement(
            urbanMean,
            backgroundMean
          );


        var periurbanEnhancement =
          calculateEnhancement(
            periurbanMean,
            backgroundMean
          );


        return ee.Feature(null, {

          Year:
            year,

          Season:
            seasonName,

          Season_Order:
            season.get('Season_Order'),

          Start_Date:
            start.format('YYYY-MM-dd'),

          End_Date:
            end.advance(-1, 'day')
               .format('YYYY-MM-dd'),

          Image_Count:
            seasonalCollection.size(),

          Valid_coverage_percent:
            coveragePercent,

          Urban_XCH4_ppb:
            urbanMean,

          Urban_XCH4_temporal_SD_ppb:
            urbanSD,

          Periurban_XCH4_ppb:
            periurbanMean,

          Background_XCH4_ppb:
            backgroundMean,

          Urban_CH4_enhancement_ppb:
            urbanEnhancement,

          Periurban_CH4_enhancement_ppb:
            periurbanEnhancement
        });
      }
    );
  }).flatten()
);


// Exclude Winter 2019.
// Retain records if year is not 2019 OR season is not Winter.
seasonalAnalysis = seasonalAnalysis.filter(

  ee.Filter.or(

    ee.Filter.neq(
      'Year',
      2019
    ),

    ee.Filter.neq(
      'Season',
      'Winter'
    )
  )
);


print(
  'Seasonal analysis — Winter 2019 excluded',
  seasonalAnalysis
);

print(
  'Seasonal record count; expected = 27',
  seasonalAnalysis.size()
);


Export.table.toDrive({

  collection:
    seasonalAnalysis,

  description:
    'Biratnagar_Seasonal_CH4_Enhancement_2019_2025',

  fileNamePrefix:
    'Biratnagar_Seasonal_CH4_Enhancement_2019_2025',

  fileFormat:
    'CSV'
});


// =====================================================================
// 12. ANNUAL ANALYSIS USING MONTHLY COMPOSITES
//
// Annual collections may exceed 5,000 images. Therefore, each annual
// value is calculated from monthly mean images.
// =====================================================================

function getAnnualMonthlyComposite(yearValue) {

  var year = ee.Number(yearValue);

  var startMonth = ee.Number(
    ee.Algorithms.If(
      year.eq(2019),
      3,
      1
    )
  );

  var monthNumbers = ee.List.sequence(
    startMonth,
    12
  );


  var monthlyImages = monthNumbers.map(
    function(monthValue) {

      var month = ee.Number(monthValue);

      var start = ee.Date.fromYMD(
        year,
        month,
        1
      );

      return getMonthlyMethaneMean(start);
    }
  );


  return ee.ImageCollection.fromImages(
    monthlyImages
  );
}


var annualAnalysis = ee.FeatureCollection(

  years.map(function(yearValue) {

    var year = ee.Number(yearValue);

    var monthlyCompositeCollection =
      getAnnualMonthlyComposite(year);


    // Each available month receives equal weight.
    var annualMean =
      monthlyCompositeCollection
        .mean()
        .rename('XCH4');


    var annualSD =
      monthlyCompositeCollection
        .reduce(ee.Reducer.stdDev())
        .rename('XCH4_monthly_sd');


    var urbanMean = regionMean(
      annualMean,
      'XCH4',
      urbanCore,
      CH4_SCALE
    );

    var periurbanMean = regionMean(
      annualMean,
      'XCH4',
      periUrbanRing,
      CH4_SCALE
    );

    var backgroundMean = regionMean(
      annualMean,
      'XCH4',
      regionalBackground,
      CH4_SCALE
    );

    var urbanSD = regionMean(
      annualSD,
      'XCH4_monthly_sd',
      urbanCore,
      CH4_SCALE
    );


    var urbanEnhancement =
      calculateEnhancement(
        urbanMean,
        backgroundMean
      );


    return ee.Feature(null, {

      Year:
        year,

      Months_Included:
        monthlyCompositeCollection.size(),

      Urban_XCH4_ppb:
        urbanMean,

      Urban_monthly_SD_ppb:
        urbanSD,

      Periurban_XCH4_ppb:
        periurbanMean,

      Background_XCH4_ppb:
        backgroundMean,

      Urban_CH4_enhancement_ppb:
        urbanEnhancement,

      Note:
        ee.Algorithms.If(
          year.eq(2019),
          'March-December only',
          'January-December'
        )
    });
  })
);


print(
  'Annual methane analysis',
  annualAnalysis
);


Export.table.toDrive({

  collection:
    annualAnalysis,

  description:
    'Biratnagar_Annual_CH4_Enhancement_2019_2025',

  fileNamePrefix:
    'Biratnagar_Annual_CH4_Enhancement_2019_2025',

  fileFormat:
    'CSV'
});


// =====================================================================
// 13. MONTHLY CH4 CHART
// =====================================================================

var monthlyChart = ui.Chart.feature.byFeature({

  features:
    monthlyValid,

  xProperty:
    'Date',

  yProperties: [
    'Urban_XCH4_mean_ppb',
    'Background_XCH4_mean_ppb'
  ]
})
.setChartType('LineChart')
.setOptions({

  title:
    'Monthly XCH4: Biratnagar and Regional Background',

  hAxis: {
    title: 'Month',
    slantedText: true,
    slantedTextAngle: 45
  },

  vAxis: {
    title: 'XCH4 concentration (ppb)'
  },

  lineWidth:
    2,

  pointSize:
    2,

  colors: [
    '#d73027',
    '#4575b4'
  ],

  legend: {
    position: 'right'
  }
});


print(monthlyChart);


// =====================================================================
// 14. MONTHLY METHANE-ENHANCEMENT CHART
// =====================================================================

var enhancementChart = ui.Chart.feature.byFeature({

  features:
    monthlyValid,

  xProperty:
    'Date',

  yProperties: [
    'Urban_CH4_enhancement_ppb'
  ]
})
.setChartType('LineChart')
.setOptions({

  title:
    'Background-Normalized Methane Enhancement over Biratnagar',

  hAxis: {
    title: 'Month',
    slantedText: true,
    slantedTextAngle: 45
  },

  vAxis: {
    title: 'Delta XCH4 (ppb)'
  },

  lineWidth:
    2,

  pointSize:
    3,

  colors: [
    '#7b3294'
  ],

  legend: {
    position: 'none'
  }
});


print(enhancementChart);


// =====================================================================
// 15. SEASONAL CHART
// =====================================================================

var seasonalChart = ui.Chart.feature.groups({

  features:
    seasonalAnalysis,

  xProperty:
    'Year',

  yProperty:
    'Urban_XCH4_ppb',

  seriesProperty:
    'Season'
})
.setChartType('LineChart')
.setOptions({

  title:
    'Seasonal Mean Methane Concentration over Biratnagar',

  hAxis: {
    title: 'Year',
    format: '####'
  },

  vAxis: {
    title: 'XCH4 concentration (ppb)'
  },

  lineWidth:
    2,

  pointSize:
    5,

  colors: [
    '#4575b4',
    '#d73027',
    '#1a9850',
    '#fdae61'
  ],

  legend: {
    position: 'right'
  }
});


print(seasonalChart);


// =====================================================================
// 16. ANNUAL CHART
// =====================================================================

var annualChart = ui.Chart.feature.byFeature({

  features:
    annualAnalysis,

  xProperty:
    'Year',

  yProperties: [
    'Urban_XCH4_ppb',
    'Background_XCH4_ppb'
  ]
})
.setChartType('LineChart')
.setOptions({

  title:
    'Annual Mean Methane over Biratnagar and Regional Background',

  hAxis: {
    title: 'Year',
    format: '####'
  },

  vAxis: {
    title: 'XCH4 concentration (ppb)'
  },

  lineWidth:
    3,

  pointSize:
    6,

  colors: [
    '#d73027',
    '#4575b4'
  ],

  legend: {
    position: 'right'
  }
});


print(annualChart);


// =====================================================================
// 17. ANNUAL MAPS
// =====================================================================

var methaneVis = {

  min:
    1850,

  max:
    2025,

  palette: [
    '313695',
    '4575b4',
    '74add1',
    'abd9e9',
    'e0f3f8',
    'ffffbf',
    'fee090',
    'fdae61',
    'f46d43',
    'd73027',
    'a50026'
  ]
};


[2019, 2021, 2023, 2025].forEach(
  function(year) {

    var annualImage =
      getAnnualMonthlyComposite(year)
        .mean()
        .select(CH4_BAND);


    Map.addLayer(

      annualImage.clip(studyDomain),

      methaneVis,

      'Annual XCH4 ' + year,

      false
    );
  }
);


// =====================================================================
// 18. MEAN 2019–2025 METHANE MAP
// =====================================================================

var annualMeanImageCollection =
  ee.ImageCollection.fromImages(

    years.map(function(yearValue) {

      var year = ee.Number(yearValue);

      return getAnnualMonthlyComposite(year)
        .mean()
        .select(CH4_BAND)
        .set({
          Year: year,
          'system:time_start':
            ee.Date.fromYMD(
              year,
              1,
              1
            ).millis()
        });
    })
  );


var meanMethane2019to2025 =
  annualMeanImageCollection
    .mean()
    .rename(CH4_BAND);


Map.addLayer(

  meanMethane2019to2025.clip(
    studyDomain
  ),

  methaneVis,

  'Mean XCH4, 2019–2025',

  true
);


// =====================================================================
// 19. INUNDATION-FREQUENCY MAP
// =====================================================================

var waterFrequency = sentinel1.map(
  function(image) {

    return image
      .lt(WATER_THRESHOLD_DB)
      .and(permanentWater.not())
      .rename('temporary_water')
      .copyProperties(
        image,
        ['system:time_start']
      );
  }
).mean();


Map.addLayer(

  waterFrequency.clip(studyDomain),

  {
    min: 0,
    max: 0.7,

    palette: [
      'ffffff',
      'ffffcc',
      'a1dab4',
      '41b6c4',
      '225ea8'
    ]
  },

  'Temporary-water frequency',

  false
);


// =====================================================================
// 20. LAND-COVER MAPS
// =====================================================================

Map.addLayer(

  croplandMask
    .selfMask()
    .clip(studyDomain),

  {
    palette: ['e6b800']
  },

  'Cropland',

  false
);


Map.addLayer(

  builtMask
    .selfMask()
    .clip(studyDomain),

  {
    palette: ['d73027']
  },

  'Built-up land',

  false
);


// =====================================================================
// 21. EXPORT MEAN METHANE MAP
// =====================================================================

Export.image.toDrive({

  image:
    meanMethane2019to2025.clip(
      studyDomain
    ),

  description:
    'Biratnagar_Mean_XCH4_2019_2025',

  fileNamePrefix:
    'Biratnagar_Mean_XCH4_2019_2025',

  region:
    studyDomain,

  scale:
    CH4_SCALE,

  maxPixels:
    1e13,

  fileFormat:
    'GeoTIFF'
});


// =====================================================================
// 22. EXPORT ANALYTICAL ZONES
// =====================================================================

Export.table.toDrive({

  collection:
    zones,

  description:
    'Biratnagar_CH4_Analytical_Zones',

  fileNamePrefix:
    'Biratnagar_CH4_Analytical_Zones',

  fileFormat:
    'SHP'
});


// =====================================================================
// 23. FINAL CONSOLE CHECKS
// =====================================================================

print(
  'Expected seasonal records = 27',
  seasonalAnalysis.size()
);

print(
  'Expected annual records = 7',
  annualAnalysis.size()
);

print(
  'Total monthly records before screening',
  monthlyAnalysis.size()
);

print(
  'Monthly records after coverage screening',
  monthlyValid.size()
);

print(
  'Analysis completed. Start exports from the Tasks tab.'
);
