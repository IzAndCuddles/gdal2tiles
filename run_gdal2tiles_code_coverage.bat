@echo off
set OLD_PATH=%PATH%
set PATH=C:\Program Files (x86)\GDAL;%PATH%
set GDAL_DATA=C:\Program Files (x86)\GDAL\gdal-data
set GDAL_DRIVER_PATH=C:\Program Files (x86)\GDAL\gdalplugins

set INPUT_IMAGE=D:\Stage\gdal2tiles\image_test_miniature\image.vrt

set chemin=%cd%

cd D:\Stage\gdal2tiles\test_image\


pushd mercator_average
echo mercator average
%chemin%\gdal2tiles.py -p mercator %INPUT_IMAGE%
popd
echo -------------------------------------------------------------
pushd mercator_near
echo mercator near
%chemin%\gdal2tiles.py -p mercator -r near %INPUT_IMAGE%
popd
echo -------------------------------------------------------------
pushd mercator_antialias
echo mercator antialias
%chemin%\gdal2tiles.py -p mercator -r antialias %INPUT_IMAGE%
popd
echo =============================================================
pushd geodetic_average
echo geodetic average
%chemin%\gdal2tiles.py -p geodetic %INPUT_IMAGE%
popd
echo -------------------------------------------------------------
pushd geodetic_near
echo geodetic near
%chemin%\gdal2tiles.py -p geodetic -r near %INPUT_IMAGE%
popd
echo -------------------------------------------------------------
pushd geodetic_antialias
echo geodetic antialias
%chemin%\gdal2tiles.py -p geodetic -r antialias %INPUT_IMAGE%
popd
echo =============================================================
pushd raster_average
echo raster average
%chemin%\gdal2tiles.py -p raster %INPUT_IMAGE%
popd
echo -------------------------------------------------------------
pushd raster_near
echo raster near
%chemin%\gdal2tiles.py -p raster -r near %INPUT_IMAGE%
popd
rem echo -------------------------------------------------------------
rem : config raster+antialias > fail de generate_base_tiles
rem pushd raster_antialias
rem %chemin%\gdal2tiles.py -p raster -r antialias %INPUT_IMAGE%
rem popd

set PATH=%OLD_PATH%
pause
