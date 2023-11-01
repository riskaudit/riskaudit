-- Tutorial
-- https://github.com/OvertureMaps/data/blob/main/webinar_queries.md

-- INSTALL spatial;
-- INSTALL httpfs;

-- Initial requirements
duckdb
LOAD spatial;
LOAD httpfs;
SET s3_region='us-west-2';

-- Getting started
SELECT COUNT(id) FROM places
SELECT * FROM places LIMIT 10

-- What is the distribution of confidence values in the places theme?
SELECT ROUND(confidence * 10)/10 AS _confidence,
	ROUND(CAST(COUNT(id) AS double) * 100 / (select count(id) from places),2) as _percent
FROM places
GROUP BY ROUND(confidence * 10)
ORDER BY _confidence DESC

-- Spatial Distribution of Places above 0.8 confidence.
WITH places_with_quadkey AS (
	SELECT bing_tile_quadkey(
			BING_TILE_AT(
				ST_Y(ST_GeomFromBinary(geometry)),
				ST_X(ST_GeomFromBinary(geometry)),
				8
			)
		) AS q8,
		id
	FROM places
	WHERE confidence > 0.8
)
SELECT BING_TILE_POLYGON(BING_TILE(q8)),
count(id) as num_places
FROM places_with_quadkey
GROUP BY q8

-- Now download some data
SELECT TRY(
		FILTER(
			names [ 'common' ],
			name->name [ 'language' ] = 'local'
		) [ 1 ] [ 'value' ]
	) AS name,
	categories.main AS category,
	confidence,
	ST_GeomFromBinary(geometry) as wkt
FROM places
WHERE confidence > 0.8
	AND bbox.minX > -126.7952
	AND bbox.maxX < -118.5453
	AND bbox.minY > 43.5453
	AND bbox.maxY < 50.4344

ogr2ogr places.shp places.csv

-- Explore global building coverage at zoom level 8 tiles.
WITH buildings_with_quadkey AS (
  SELECT
    bing_tile_quadkey(
      BING_TILE_AT(
        (bbox.maxY + bbox.minY)/2,
        (bbox.maxX + bbox.minX)/2,
        8
      )
    ) AS q8,
    id,
    CARDINALITY(
FILTER(sources, x -> x['dataset'] = 'OpenStreetMap')
    )>0 AS osm_building
  FROM buildings
)
SELECT
  BING_TILE_POLYGON(BING_TILE(q8)),
  COUNT(id) as num_buildings,
  COUNT_IF(osm_building) AS osm_building,
  COUNT_IF(osm_building) / CAST(COUNT(id) AS double) AS percent_osm
FROM buildings_with_quadkey
GROUP BY q8

-- Try downloading a small subset of buildings around Seattle.
SELECT class, height,
	ST_GeomFromBinary(geometry) as wkt
FROM buildings
WHERE ST_CONTAINS(
		ST_GeometryFromText(
			'POLYGON((-122.36719284258956 47.618321237733284,-122.33594394153602 47.632404470851924,-122.2775808079059 47.61236859966664,-122.34462990362489 47.58012171471199,-122.36719284258956 47.618321237733284))'
		),
		ST_GeomFromBinary(geometry)
	)

-- Save the following file as a buildings.csv and then convert it to a GeoJSON:
ogr2ogr -select height,class buildings.geojson buildings.csv

-- Then try using tippecanoe to create tiles
tippecanoe -z13 -Z13 -fo buildings.pmtiles buildings.geojson





LOAD httpfs;
LOAD spatial;
SET s3_region='us-west-2';

COPY (
    SELECT
        type,
        version,
        CAST(updatetime as varchar) as updateTime,
        height,
        numfloors as numFloors,
        level,
        class,
        JSON(names) as names,
        JSON(sources) as sources,
        ST_GeomFromWKB(geometry) as geometry
    FROM
        read_parquet('s3://overturemaps-us-west-2/release/2023-10-19-alpha.0/theme=buildings/type=*/*', hive_partitioning=1)
    LIMIT
        100
    ) TO 'buildings_sample.geojsonseq'
WITH (FORMAT GDAL, DRIVER 'GeoJSONSeq');