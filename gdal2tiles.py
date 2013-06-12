#!/usr/bin/env python
#******************************************************************************
#  $Id: gdal2tiles.py 19288 2010-04-02 18:36:17Z rouault $
# 
# Project:  Google Summer of Code 2007, 2008 (http://code.google.com/soc/)
# Support:  BRGM (http://www.brgm.fr)
# Purpose:  Convert a raster into TMS (Tile Map Service) tiles in a directory.
#           - generate Google Earth metadata (KML SuperOverlay)
#           - generate simple HTML viewer based on Google Maps and OpenLayers
#           - support of global tiles (Spherical Mercator) for compatibility
#               with interactive web maps a la Google Maps
# Author:   Klokan Petr Pridal, klokan at klokan dot cz
# Web:      http://www.klokan.cz/projects/gdal2tiles/
# GUI:      http://www.maptiler.org/
#
###############################################################################
# Copyright (c) 2008, Klokan Petr Pridal
# 
#  Permission is hereby granted, free of charge, to any person obtaining a
#  copy of this software and associated documentation files (the "Software"),
#  to deal in the Software without restriction, including without limitation
#  the rights to use, copy, modify, merge, publish, distribute, sublicense,
#  and/or sell copies of the Software, and to permit persons to whom the
#  Software is furnished to do so, subject to the following conditions:
# 
#  The above copyright notice and this permission notice shall be included
#  in all copies or substantial portions of the Software.
# 
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
#  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
#  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#  DEALINGS IN THE SOFTWARE.
#******************************************************************************

# TODO essayer de mettre des docstring pour les fonctions / methodes / classes
# TODO faire essais sur comportement en cas d'erreur pendant traitement d'un/plusieurs jobs
# TODO supprimer tous nos TODO en francais

import sys
import multiprocessing

try:
    from osgeo import gdal
    from osgeo import osr
except:
    import gdal
    print('You are using "old gen" bindings. gdal2tiles needs "new gen" bindings.')
    sys.exit(1)

import os
import math

try:
    from PIL import Image
    import numpy
    import osgeo.gdal_array as gdalarray
except:
    # 'antialias' resampling is not available
    pass

__version__ = "$Id: gdal2tiles.py 19288 2010-04-02 18:36:17Z rouault $"

resampling_list = ('average','near','bilinear','cubic','cubicspline','lanczos','antialias')
profile_list = ('mercator','geodetic','raster') #,'zoomify')
webviewer_list = ('all','google','openlayers','none')

# =============================================================================
# =============================================================================
# =============================================================================

__doc__globalmaptiles = """
globalmaptiles.py

Global Map Tiles as defined in Tile Map Service (TMS) Profiles
==============================================================

Functions necessary for generation of global tiles used on the web.
It contains classes implementing coordinate conversions for:

  - GlobalMercator (based on EPSG:900913 = EPSG:3785)
       for Google Maps, Yahoo Maps, Microsoft Maps compatible tiles
  - GlobalGeodetic (based on EPSG:4326)
       for OpenLayers Base Map and Google Earth compatible tiles

More info at:

http://wiki.osgeo.org/wiki/Tile_Map_Service_Specification
http://wiki.osgeo.org/wiki/WMS_Tiling_Client_Recommendation
http://msdn.microsoft.com/en-us/library/bb259689.aspx
http://code.google.com/apis/maps/documentation/overlays.html#Google_Maps_Coordinates

Created by Klokan Petr Pridal on 2008-07-03.
Google Summer of Code 2008, project GDAL2Tiles for OSGEO.

In case you use this class in your product, translate it to another language
or find it usefull for your project please let me know.
My email: klokan at klokan dot cz.
I would like to know where it was used.

Class is available under the open-source GDAL license (www.gdal.org).
"""

MAXZOOMLEVEL = 32
TILESIZE = 256


class GlobalMercator(object):
    """
    TMS Global Mercator Profile
    ---------------------------

    Functions necessary for generation of tiles in Spherical Mercator projection,
    EPSG:900913 (EPSG:gOOglE, Google Maps Global Mercator), EPSG:3785, OSGEO:41001.

    Such tiles are compatible with Google Maps, Microsoft Virtual Earth, Yahoo Maps,
    UK Ordnance Survey OpenSpace API, ...
    and you can overlay them on top of base maps of those web mapping applications.
    
    Pixel and tile coordinates are in TMS notation (origin [0,0] in bottom-left).

    What coordinate conversions do we need for TMS Global Mercator tiles::

         LatLon      <->       Meters      <->     Pixels    <->       Tile     

     WGS84 coordinates   Spherical Mercator  Pixels in pyramid  Tiles in pyramid
         lat/lon            XY in metres     XY pixels Z zoom      XYZ from TMS 
        EPSG:4326           EPSG:900913                                         
         .----.              ---------               --                TMS      
        /      \     <->     |       |     <->     /----/    <->      Google    
        \      /             |       |           /--------/          QuadTree   
         -----               ---------         /------------/                   
       KML, public         WebMapService         Web Clients      TileMapService

    What is the coordinate extent of Earth in EPSG:900913?

      [-20037508.342789244, -20037508.342789244, 20037508.342789244, 20037508.342789244]
      Constant 20037508.342789244 comes from the circumference of the Earth in meters,
      which is 40 thousand kilometers, the coordinate origin is in the middle of extent.
      In fact you can calculate the constant as: 2 * math.pi * 6378137 / 2.0
      $ echo 180 85 | gdaltransform -s_srs EPSG:4326 -t_srs EPSG:900913
      Polar areas with abs(latitude) bigger then 85.05112878 are clipped off.

    What are zoom level constants (pixels/meter) for pyramid with EPSG:900913?

      whole region is on top of pyramid (zoom=0) covered by 256x256 pixels tile,
      every lower zoom level resolution is always divided by two
      initialResolution = 20037508.342789244 * 2 / 256 = 156543.03392804062

    What is the difference between TMS and Google Maps/QuadTree tile name convention?

      The tile raster itself is the same (equal extent, projection, pixel size),
      there is just different identification of the same raster tile.
      Tiles in TMS are counted from [0,0] in the bottom-left corner, id is XYZ.
      Google placed the origin [0,0] to the top-left corner, reference is XYZ.
      Microsoft is referencing tiles by a QuadTree name, defined on the website:
      http://msdn2.microsoft.com/en-us/library/bb259689.aspx

    The lat/lon coordinates are using WGS84 datum, yeh?

      Yes, all lat/lon we are mentioning should use WGS84 Geodetic Datum.
      Well, the web clients like Google Maps are projecting those coordinates by
      Spherical Mercator, so in fact lat/lon coordinates on sphere are treated as if
      the were on the WGS84 ellipsoid.
     
      From MSDN documentation:
      To simplify the calculations, we use the spherical form of projection, not
      the ellipsoidal form. Since the projection is used only for map display,
      and not for displaying numeric coordinates, we don't need the extra precision
      of an ellipsoidal projection. The spherical projection causes approximately
      0.33 percent scale distortion in the Y direction, which is not visually noticable.

    How do I create a raster in EPSG:900913 and convert coordinates with PROJ.4?

      You can use standard GIS tools like gdalwarp, cs2cs or gdaltransform.
      All of the tools supports -t_srs 'epsg:900913'.

      For other GIS programs check the exact definition of the projection:
      More info at http://spatialreference.org/ref/user/google-projection/
      The same projection is degined as EPSG:3785. WKT definition is in the official
      EPSG database.

      Proj4 Text:
        +proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0
        +k=1.0 +units=m +nadgrids=@null +no_defs

      Human readable WKT format of EPGS:900913:
         PROJCS["Google Maps Global Mercator",
             GEOGCS["WGS 84",
                 DATUM["WGS_1984",
                     SPHEROID["WGS 84",6378137,298.257223563,
                         AUTHORITY["EPSG","7030"]],
                     AUTHORITY["EPSG","6326"]],
                 PRIMEM["Greenwich",0],
                 UNIT["degree",0.0174532925199433],
                 AUTHORITY["EPSG","4326"]],
             PROJECTION["Mercator_1SP"],
             PARAMETER["central_meridian",0],
             PARAMETER["scale_factor",1],
             PARAMETER["false_easting",0],
             PARAMETER["false_northing",0],
             UNIT["metre",1,
                 AUTHORITY["EPSG","9001"]]]
    """

    def __init__(self, tileSize=256):
        "Initialize the TMS Global Mercator pyramid"
        self.tileSize = tileSize
        self.initialResolution = 2 * math.pi * 6378137 / self.tileSize
        # 156543.03392804062 for tileSize 256 pixels
        self.originShift = 2 * math.pi * 6378137 / 2.0
        # 20037508.342789244

    def LatLonToMeters(self, lat, lon ):
        "Converts given lat/lon in WGS84 Datum to XY in Spherical Mercator EPSG:900913"

        mx = lon * self.originShift / 180.0
        my = math.log( math.tan((90 + lat) * math.pi / 360.0 )) / (math.pi / 180.0)

        my = my * self.originShift / 180.0
        return mx, my

    def MetersToLatLon(self, mx, my ):
        "Converts XY point from Spherical Mercator EPSG:900913 to lat/lon in WGS84 Datum"

        lon = (mx / self.originShift) * 180.0
        lat = (my / self.originShift) * 180.0

        lat = 180 / math.pi * (2 * math.atan( math.exp( lat * math.pi / 180.0)) - math.pi / 2.0)
        return lat, lon

    def PixelsToMeters(self, px, py, zoom):
        "Converts pixel coordinates in given zoom level of pyramid to EPSG:900913"

        res = self.Resolution( zoom )
        mx = px * res - self.originShift
        my = py * res - self.originShift
        return mx, my
        
    def MetersToPixels(self, mx, my, zoom):
        "Converts EPSG:900913 to pyramid pixel coordinates in given zoom level"
                
        res = self.Resolution( zoom )
        px = (mx + self.originShift) / res
        py = (my + self.originShift) / res
        return px, py
    
    def PixelsToTile(self, px, py):
        "Returns a tile covering region in given pixel coordinates"

        tx = int( math.ceil( px / float(self.tileSize) ) - 1 )
        ty = int( math.ceil( py / float(self.tileSize) ) - 1 )
        return tx, ty

    def PixelsToRaster(self, px, py, zoom):
        "Move the origin of pixel coordinates to top-left corner"
        
        mapSize = self.tileSize << zoom
        return px, mapSize - py
        
    def MetersToTile(self, mx, my, zoom):
        "Returns tile for given mercator coordinates"
        
        px, py = self.MetersToPixels( mx, my, zoom)
        return self.PixelsToTile( px, py)

    def TileBounds(self, tx, ty, zoom):
        "Returns bounds of the given tile in EPSG:900913 coordinates"
        
        minx, miny = self.PixelsToMeters( tx*self.tileSize, ty*self.tileSize, zoom )
        maxx, maxy = self.PixelsToMeters( (tx+1)*self.tileSize, (ty+1)*self.tileSize, zoom )
        return ( minx, miny, maxx, maxy )

    def TileLatLonBounds(self, tx, ty, zoom ):
        "Returns bounds of the given tile in latutude/longitude using WGS84 datum"

        bounds = self.TileBounds( tx, ty, zoom)
        minLat, minLon = self.MetersToLatLon(bounds[0], bounds[1])
        maxLat, maxLon = self.MetersToLatLon(bounds[2], bounds[3])
         
        return ( minLat, minLon, maxLat, maxLon )
        
    def Resolution(self, zoom ):
        "Resolution (meters/pixel) for given zoom level (measured at Equator)"
        
        # return (2 * math.pi * 6378137) / (self.tileSize * 2**zoom)
        return self.initialResolution / (2**zoom)
        
    def ZoomForPixelSize(self, pixelSize ):
        "Maximal scaledown zoom of the pyramid closest to the pixelSize."
        
        for i in range(MAXZOOMLEVEL):
            if pixelSize > self.Resolution(i):
                if i!=0:
                    return i-1
                else:
                    return 0 # We don't want to scale up
        
    def GoogleTile(self, tx, ty, zoom):
        "Converts TMS tile coordinates to Google Tile coordinates"
        
        # coordinate origin is moved from bottom-left to top-left corner of the extent
        return tx, (2**zoom - 1) - ty

    def QuadTree(self, tx, ty, zoom ):
        "Converts TMS tile coordinates to Microsoft QuadTree"
        
        quadKey = ""
        ty = (2**zoom - 1) - ty
        for i in range(zoom, 0, -1):
            digit = 0
            mask = 1 << (i-1)
            if (tx & mask) != 0:
                digit += 1
            if (ty & mask) != 0:
                digit += 2
            quadKey += str(digit)
            
        return quadKey

#---------------------

class GlobalGeodetic(object):
    """
    TMS Global Geodetic Profile
    ---------------------------

    Functions necessary for generation of global tiles in Plate Carre projection,
    EPSG:4326, "unprojected profile".

    Such tiles are compatible with Google Earth (as any other EPSG:4326 rasters)
    and you can overlay the tiles on top of OpenLayers base map.
    
    Pixel and tile coordinates are in TMS notation (origin [0,0] in bottom-left).

    What coordinate conversions do we need for TMS Global Geodetic tiles?

      Global Geodetic tiles are using geodetic coordinates (latitude,longitude)
      directly as planar coordinates XY (it is also called Unprojected or Plate
      Carre). We need only scaling to pixel pyramid and cutting to tiles.
      Pyramid has on top level two tiles, so it is not square but rectangle.
      Area [-180,-90,180,90] is scaled to 512x256 pixels.
      TMS has coordinate origin (for pixels and tiles) in bottom-left corner.
      Rasters are in EPSG:4326 and therefore are compatible with Google Earth.

         LatLon      <->      Pixels      <->     Tiles     

     WGS84 coordinates   Pixels in pyramid  Tiles in pyramid
         lat/lon         XY pixels Z zoom      XYZ from TMS 
        EPSG:4326                                           
         .----.                ----                         
        /      \     <->    /--------/    <->      TMS      
        \      /         /--------------/                   
         -----        /--------------------/                
       WMS, KML    Web Clients, Google Earth  TileMapService
    """

    def __init__(self, tileSize = 256):
        self.tileSize = tileSize

    def LatLonToPixels(self, lat, lon, zoom):
        "Converts lat/lon to pixel coordinates in given zoom of the EPSG:4326 pyramid"

        res = 180.0 / self.tileSize / 2**zoom
        px = (180 + lat) / res
        py = (90 + lon) / res
        return px, py

    def PixelsToTile(self, px, py):
        "Returns coordinates of the tile covering region in pixel coordinates"

        tx = int( math.ceil( px / float(self.tileSize) ) - 1 )
        ty = int( math.ceil( py / float(self.tileSize) ) - 1 )
        return tx, ty
    
    def LatLonToTile(self, lat, lon, zoom):
        "Returns the tile for zoom which covers given lat/lon coordinates"
        
        px, py = self.LatLonToPixels( lat, lon, zoom)
        return self.PixelsToTile(px,py)

    def Resolution(self, zoom ):
        "Resolution (arc/pixel) for given zoom level (measured at Equator)"
        
        return 180.0 / self.tileSize / 2**zoom
        #return 180 / float( 1 << (8+zoom) )
        
    def ZoomForPixelSize(self, pixelSize ):
        "Maximal scaledown zoom of the pyramid closest to the pixelSize."

        for i in range(MAXZOOMLEVEL):
            if pixelSize > self.Resolution(i):
                if i!=0:
                    return i-1
                else:
                    return 0 # We don't want to scale up

    def TileBounds(self, tx, ty, zoom):
        "Returns bounds of the given tile"
        res = 180.0 / self.tileSize / 2**zoom
        return (
            tx*self.tileSize*res - 180,
            ty*self.tileSize*res - 90,
            (tx+1)*self.tileSize*res - 180,
            (ty+1)*self.tileSize*res - 90
        )
        
    def TileLatLonBounds(self, tx, ty, zoom):
        "Returns bounds of the given tile in the SWNE form"
        b = self.TileBounds(tx, ty, zoom)
        return (b[1],b[0],b[3],b[2])


def processBaseTileJobs(q, s_srs, input, profile, verbose, in_nodata, n):
    "Create an only 'out_data' for the process that targets the function, then get jobs from the queue and create base tiles"
    out_data = init_out_data(s_srs, input, profile, verbose, in_nodata)
    job = q.get()
    
    while job !='TERMINATE':
        tiledriver, options, resampling, tilebands, querysize, tz, ty, lx, tile,output,tileext,tmaxx,tmaxy,mercator,geodetic = job
        for tx in lx:
            process_base_tile_job(out_data,tiledriver, options, resampling, tilebands, querysize, tz, ty, tx, tile,output,tileext,tmaxx,tmaxy,mercator,geodetic)
            n.value += 1        # _not_ synchronized / not accurate but faster than proper lock
            
        job = q.get()
        
def processOverviewTileJobs(q,n):
    "Get jobs from the queue and create base tiles. Handle case were lx is a list, and when it is a value"
    job = q.get()
    
    while job !='TERMINATE':
        output,options,tiledriver,resampling,tileext, tile, tilebands, tz, ty, lx = job
        try:
            iterator = iter(lx)
        except TypeError:
            lx=[lx]
        for tx in lx:
            process_overview_tile_job(output,options,tiledriver,resampling,tileext, tile, tilebands, tz, ty, tx)
            n.value += 1        # _not_ synchronized / not accurate but faster than proper lock
        job = q.get()

class MultiProcess(object):
    "The MultiProcess class share a queue between a list of processes"
    def __init__(self,num_process):
        self.q=multiprocessing.Queue()
        self.process=[]
        self.num_process = num_process
    
    def _processTile(self,target,args):
        "Create processes and start them"
        for i in xrange(self.num_process):
            p = multiprocessing.Process(target=target, args=args)
            p.daemon=True
            p.start()
            self.process.append(p)
    
    def processBase(self, s_srs, input, profile, verbose, in_nodata, n):
        "Create processes which will create base tiles"
        self._processTile(processBaseTileJobs,(self.q,s_srs,input,profile,verbose,in_nodata,n))
            
    def processOverview(self, n):
        "Create processes which will create overview tiles"
        self._processTile(processOverviewTileJobs,(self.q,n))
            
    def sendJob(self, job):
        "Put a job to the queue"
        self.q.put(job)

    def finish(self,n,ncount):
        "Finish all processes, and print the progressbar while waiting for all processes to join"
        for p in self.process:
            self.q.put('TERMINATE')
        # In order to keep the progressbar as it was in the original version, we can't
        # just join(0) the processes, or use a JoinableQueue : a timeout is needed.
        # While all the processes aren't done, we join each process for 1 sec, than 
        # call the progressbar.
        done = False
        while not done:
            alive_count = len(self.process)
            for p in self.process:
                p.join(1.0)
                alive = p.is_alive()
                if not alive:
                    alive_count -= 1
                progressbar(n.value/float(ncount))
            done = (alive_count == 0)


class Profile(object):
    "The Profile class gathers the profile-related variables"
    def __init__(self):
        self.mercator=GlobalMercator()
        self.geodetic=GlobalGeodetic()
        self.tileswne=None

class Tile(object):
    "The Tile class gathers the tile-related variables"
    def __init__(self):
        self.tminz=None
        self.tmaxz=None
        self.tminmax=None
        self.tsize=None
        self.nativezoom=None
        self.ominx=None
        self.omaxx=None
        self.ominy=None
        self.omaxy=None
        
        
class OutData(object):
    "The OutData class gathers the variables related to the output dataset"
    def __init__(self):
        self.out_ds=None
        self.alphaband=None
        self.dataBandsCount=None
        self.out_gt=None

# The following functions come directly from the GDAL2Tiles.__init__ and 
# GDAL2Tiles.open_input methods of the original version of GDAL2tiles.
# To make the code easier to understand, and because of some unpickleable
# elements, we created these new functions. They are called in
# Configuration.__init__ and Configuration.open_input, and during the
# multiprocessing of the tiles, to create unpickleable elements needed to
# create the tiles, from pickleable arguments.

def def_in_srs(s_srs,input):
    in_ds=gdal.Open(input, gdal.GA_ReadOnly)
    in_srs = None
    if s_srs:
        in_srs = osr.SpatialReference()
        in_srs.SetFromUserInput(s_srs)
        in_srs_wkt = in_srs.ExportToWkt()
    else:
        in_srs_wkt = in_ds.GetProjection()
        if not in_srs_wkt and in_ds.GetGCPCount() != 0:
            in_srs_wkt = in_ds.GetGCPProjection()
        if in_srs_wkt:
            in_srs = osr.SpatialReference()
            in_srs.ImportFromWkt(in_srs_wkt)
    return in_srs,in_srs_wkt

def def_out_srs(profile, in_srs):
    out_srs = osr.SpatialReference()
    if profile == 'mercator':
        out_srs.ImportFromEPSG(900913)
    elif profile == 'geodetic':
        out_srs.ImportFromEPSG(4326)
    else:
        out_srs = in_srs
    return out_srs

def init_drv(tiledriver):
    out_drv = gdal.GetDriverByName(tiledriver)
    mem_drv = gdal.GetDriverByName('MEM')
    return out_drv,mem_drv

def init_out_data(s_srs,input,profile,verbose,in_nodata):
    in_srs,in_srs_wkt=def_in_srs(s_srs,input)
    out_srs=def_out_srs(profile,in_srs)
    in_ds=gdal.Open(input, gdal.GA_ReadOnly)
    out_data = OutData()
    out_data.out_ds = in_ds
    if profile in ('mercator', 'geodetic'):
        out_ds = None
        if (in_ds.GetGeoTransform() == (0.0, 1.0, 0.0, 0.0, 0.0, 1.0)) and (in_ds.GetGCPCount() == 0):
            error("There is no georeference - neither affine transformation (worldfile) nor GCPs. You can generate only 'raster' profile tiles.", "Either gdal2tiles with parameter -p 'raster' or use another GIS software for georeference e.g. gdal_transform -gcp / -a_ullr / -a_srs")
        if in_srs:
            if (in_srs.ExportToProj4() != out_srs.ExportToProj4()) or (in_ds.GetGCPCount() != 0):
                out_ds = image_warping(in_ds,in_srs_wkt,out_srs,verbose,in_nodata)
        else:
            error("Input file has unknown SRS.", "Use --s_srs ESPG:xyz (or similar) to provide source reference system.")
        if out_ds is not None:
            out_data.out_ds = out_ds
            if verbose:
                print "Projected file:", "tiles.vrt", "( %sP x %sL - %s bands)" % (out_ds.RasterXSize, out_ds.RasterYSize, out_ds.RasterCount)
    #
    # Here we should have a raster (out_ds) in the correct Spatial Reference system
    #
    # Get alpha band (either directly or from NODATA value)
    out_data.alphaband = out_data.out_ds.GetRasterBand(1).GetMaskBand()
    if (out_data.alphaband.GetMaskFlags() & gdal.GMF_ALPHA) or out_data.out_ds.RasterCount == 4 or out_data.out_ds.RasterCount == 2: # TODO: Better test for alpha band in the dataset
        out_data.dataBandsCount = out_data.out_ds.RasterCount - 1
    else:
        out_data.dataBandsCount = out_data.out_ds.RasterCount
    # Read the georeference
    out_data.out_gt = out_data.out_ds.GetGeoTransform()
    return out_data


def correct_ACW_nodata(in_nodata,verbose, out_ds):
    import tempfile
    tempfilename = tempfile.mktemp('-gdal2tiles.vrt')
    out_ds.GetDriver().CreateCopy(tempfilename, out_ds) # open as a text file
    s = open(tempfilename).read() # Add the warping options
    s = s.replace("""<GDALWarpOptions>""", """<GDALWarpOptions>
          <Option name="INIT_DEST">NO_DATA</Option>
          <Option name="UNIFIED_SRC_NODATA">YES</Option>""")
    # replace BandMapping tag for NODATA bands....
    for i in range(len(in_nodata)):
        s = s.replace("""<BandMapping src="%i" dst="%i"/>""" % ((i + 1), (i + 1)), 
            """<BandMapping src="%i" dst="%i">
  <SrcNoDataReal>%i</SrcNoDataReal>
  <SrcNoDataImag>0</SrcNoDataImag>
  <DstNoDataReal>%i</DstNoDataReal>
  <DstNoDataImag>0</DstNoDataImag>
</BandMapping>""" % ((i + 1), (i + 1), in_nodata[i], in_nodata[i])) # Or rewrite to white by: , 255 ))
    
# save the corrected VRT
    open(tempfilename, "w").write(s) # open by GDAL as out_ds
    out_ds = gdal.Open(tempfilename) #, gdal.GA_ReadOnly)
# delete the temporary file
    os.unlink(tempfilename) # set NODATA_VALUE metadata
    out_ds.SetMetadataItem('NODATA_VALUES', '%i %i %i' % (in_nodata[0], in_nodata[1], in_nodata[2]))
    if verbose:
        print "Modified warping result saved into 'tiles1.vrt'"
        open("tiles1.vrt", "w").write(s)
    return out_ds



def correct_ACW_mono_rgb(verbose, out_ds):
    import tempfile
    tempfilename = tempfile.mktemp('-gdal2tiles.vrt')
    out_ds.GetDriver().CreateCopy(tempfilename, out_ds) # open as a text file
    s = open(tempfilename).read() # Add the warping options
    s = s.replace("""<BlockXSize>""", 
        """<VRTRasterBand dataType="Byte" band="%i" subClass="VRTWarpedRasterBand">
        <ColorInterp>Alpha</ColorInterp>
      </VRTRasterBand>
      <BlockXSize>""" % (out_ds.RasterCount + 1))
    s = s.replace("""</GDALWarpOptions>""", 
        """<DstAlphaBand>%i</DstAlphaBand>
      </GDALWarpOptions>""" % (out_ds.RasterCount + 1))
    s = s.replace("""</WorkingDataType>""", """</WorkingDataType>
        <Option name="INIT_DEST">0</Option>""")
# save the corrected VRT
    open(tempfilename, "w").write(s) # open by GDAL as out_ds
    out_ds = gdal.Open(tempfilename) #, gdal.GA_ReadOnly)
# delete the temporary file
    os.unlink(tempfilename)
    if verbose:
        print "Modified -dstalpha warping result saved into 'tiles1.vrt'"
        open("tiles1.vrt", "w").write(s)
    return out_ds


def image_warping(in_ds,in_srs_wkt,out_srs,verbose,in_nodata):
    # Generation of VRT dataset in tile projection, default 'nearest neighbour' warping
    out_ds = gdal.AutoCreateWarpedVRT(in_ds, in_srs_wkt, out_srs.ExportToWkt())
# TODO: HIGH PRIORITY: Correction of AutoCreateWarpedVRT according the max zoomlevel for correct direct warping!!!
    if verbose:
        print "Warping of the raster by AutoCreateWarpedVRT (result saved into 'tiles.vrt')"
        out_ds.GetDriver().CreateCopy("tiles.vrt", out_ds)
# Note: self.in_srs and self.in_srs_wkt contain still the non-warped reference system!!!
# Correction of AutoCreateWarpedVRT for NODATA values
    if in_nodata != []:
        out_ds = correct_ACW_nodata(in_nodata,verbose, out_ds)

# Correction of AutoCreateWarpedVRT for Mono (1 band) and RGB (3 bands) files without NODATA:
# equivalent of gdalwarp -dstalpha
    if in_nodata == [] and out_ds.RasterCount in [1, 3]:
        out_ds = correct_ACW_mono_rgb(verbose,out_ds)
    s = '''
'''
    return out_ds


class Configuration (object):
    "The Configuration class gathers the option-related variables, and the methods initializing variables needed to create the tiles"
    
    def create_tile(self):
        tile = Tile()
    # User specified zoom levels
        if self.options.zoom:
            minmax = self.options.zoom.split('-', 1)
            minmax.extend([''])
            min_, max_ = minmax[:2]
            tile.tminz = int(min_)
            if max:
                tile.tmaxz = int(max_)
            else:
                tile.tmaxz = int(min_)
        return tile

    def __init__(self,arguments):
        
        self.input = None
        self.output = None
        
        self.querysize_c = 4 * TILESIZE
        
        # RUN THE ARGUMENT PARSER:
        parser=self.optparse_init()
        self.options, args = parser.parse_args(arguments)
        
        if not args:
            error("No input file specified")
            
        # POSTPROCESSING OF PARSED ARGUMENTS:
        
        # Tile format
        # TODO : Passer ces deux variables en options dans optparse_init
        if not self.options.format:
            self.tiledriver = 'PNG'
            self.tileext = 'png'
        else:
            self.tiledriver = self.options.format
            self.tileext = self.tiledriver.lower()
        
        # Number of processes
        nb_process = self.options.nb_process
        num_cpus = multiprocessing.cpu_count()
        if nb_process<=0:
            nb_process=max(1, num_cpus+nb_process)
        self.nb_process = nb_process
        
        # Workaround for old versions of GDAL
        try:
            if (self.options.verbose and self.options.resampling == 'near') or gdal.TermProgress_nocb:
                pass
        except:
            error("This version of GDAL is not supported. Please upgrade to 1.6+.")
        
        # Test output directory, if it doesn't exist
        if os.path.isdir(args[-1]) or ( len(args) > 1 and not os.path.exists(args[-1])):
            self.output = args[-1]
            args = args[:-1]
        
        # More files on the input not directly supported yet
        
        if (len(args) > 1):
            error("Processing of several input files is not supported.",
            """Please first use a tool like gdal_vrtmerge.py or gdal_merge.py on the files:
gdal_vrtmerge.py -o merged.vrt %s""" % " ".join(args))
            # TODO: Call functions from gdal_vrtmerge.py directly
            
        self.input = args[0]
        
        # Default values for not given options
        
        if not self.output:
            # Directory with input filename without extension in actual directory
            self.output = os.path.splitext(os.path.basename( self.input ))[0]
                
        if not self.options.title:
            self.options.title = os.path.basename( self.input )

        if self.options.url and not self.options.url.endswith('/'):
            self.options.url += '/'
        if self.options.url:
            self.options.url += os.path.basename( self.output ) + '/'
        
        self.resampling = None
        
        self.init_resampling()
        
        # KML generation
        self.kml = self.options.kml

        # Output the results

        if self.options.verbose:
            print("Options:", self.options)
            print("Input:", self.input)
            print("Output:", self.output)
            print("Cache: %s MB" % (gdal.GetCacheMax() / 1024 / 1024))
            print('')
            
    def optparse_init(self):
        """Prepare the option parser for input (argv)"""
        # TODO : ajouter les options
        
        from optparse import OptionParser, OptionGroup
        usage = "Usage: %prog [options] input_file(s) [output]"
        p = OptionParser(usage, version="%prog "+ __version__)
        p.add_option("-p", "--profile", dest='profile', type='choice', choices=profile_list,
                          help="Tile cutting profile (%s) - default 'mercator' (Google Maps compatible)" % ",".join(profile_list))
        p.add_option("-r", "--resampling", dest="resampling", type='choice', choices=resampling_list,
                        help="Resampling method (%s) - default 'average'" % ",".join(resampling_list))
        p.add_option('-s', '--s_srs', dest="s_srs", metavar="SRS",
                          help="The spatial reference system used for the source input data")
        p.add_option('-z', '--zoom', dest="zoom",
                          help="Zoom levels to render (format:'2-5' or '10').")
        p.add_option('-e', '--resume', dest="resume", action="store_true",
                          help="Resume mode. Generate only missing files.")
        p.add_option('-a', '--srcnodata', dest="srcnodata", metavar="NODATA",
                          help="NODATA transparency value to assign to the input data")
        p.add_option("-v", "--verbose",
                          action="store_true", dest="verbose",
                          help="Print status messages to stdout")
        p.add_option('-f','--format', dest='format', type='string', help="Format of the tiles - default 'png', see available formats here : http://www.gdal.org/formats_list.html")
        p.add_option('-b','--nb_process', dest ='nb_process',type='int', default=0, help="The number of processes that will be running - default : the number of cores in the computer")
        
        # KML options 
        g = OptionGroup(p, "KML (Google Earth) options", "Options for generated Google Earth SuperOverlay metadata")
        g.add_option("-k", "--force-kml", dest='kml', action="store_true",
                          help="Generate KML for Google Earth - default for 'geodetic' profile and 'raster' in EPSG:4326. For a dataset with different projection use with caution!")
        g.add_option("-n", "--no-kml", dest='kml', action="store_false",
                          help="Avoid automatic generation of KML files for EPSG:4326")
        g.add_option("-u", "--url", dest='url',
                          help="URL address where the generated tiles are going to be published")
        p.add_option_group(g)

        # HTML options
        g = OptionGroup(p, "Web viewer options", "Options for generated HTML viewers a la Google Maps")
        g.add_option("-w", "--webviewer", dest='webviewer', type='choice', choices=webviewer_list,
                          help="Web viewer to generate (%s) - default 'all'" % ",".join(webviewer_list))
        g.add_option("-t", "--title", dest='title',
                          help="Title of the map")
        g.add_option("-c", "--copyright", dest='copyright',
                          help="Copyright for the map")
        g.add_option("-g", "--googlekey", dest='googlekey',
                          help="Google Maps API key from http://code.google.com/apis/maps/signup.html")
        g.add_option("-y", "--yahookey", dest='yahookey',
                          help="Yahoo Application ID from http://developer.yahoo.com/wsregapp/")
        p.add_option_group(g)
        
        # TODO: MapFile + TileIndexes per zoom level for efficient MapServer WMS
        #g = OptionGroup(p, "WMS MapServer metadata", "Options for generated mapfile and tileindexes for MapServer")
        #g.add_option("-i", "--tileindex", dest='wms', action="store_true"
        #                 help="Generate tileindex and mapfile for MapServer (WMS)")
        # p.add_option_group(g)

        p.set_defaults(verbose=False, profile="mercator", kml=False, url='',
        webviewer='all', copyright='', resampling='average', resume=False,
        googlekey='INSERT_YOUR_KEY_HERE', yahookey='INSERT_YOUR_YAHOO_APP_ID_HERE')

        return p


    def init_resampling(self):
            if self.options.resampling == 'average':
                try:
                    if gdal.RegenerateOverview:
                        pass
                except:
                    error("'average' resampling algorithm is not available.", "Please use -r 'near' argument or upgrade to newer version of GDAL.")
            elif self.options.resampling == 'antialias':
                try:
                    if numpy:
                        pass
                except:
                    error("'antialias' resampling algorithm is not available.", "Install PIL (Python Imaging Library) and numpy.")
            elif self.options.resampling == 'near':
                self.resampling = gdal.GRA_NearestNeighbour
                self.querysize_c = TILESIZE
            elif self.options.resampling == 'bilinear':
                self.resampling = gdal.GRA_Bilinear
                self.querysize_c = TILESIZE * 2
            elif self.options.resampling == 'cubic':
                self.resampling = gdal.GRA_Cubic
            elif self.options.resampling == 'cubicspline':
                self.resampling = gdal.GRA_CubicSpline
            elif self.options.resampling == 'lanczos':
                self.resampling = gdal.GRA_Lanczos

    def open_input(self,tile):
        """Initialization of the input raster, reprojection if necessary"""
        
        gdal.AllRegister()

        # Initialize necessary GDAL drivers
        
        self.out_drv,self.mem_drv=init_drv(self.tiledriver)
        
        if not self.out_drv:
            raise Exception("The '%s' driver was not found, is it available in this GDAL build?", self.tiledriver)
        if not self.mem_drv:
            raise Exception("The 'MEM' driver was not found, is it available in this GDAL build?")
        
        # Open the input file
        
        if self.input:
            self.in_ds = gdal.Open(self.input, gdal.GA_ReadOnly)
        else:
            raise Exception("No input file was specified")

        if self.options.verbose:
            print("Input file:", "( %sP x %sL - %s bands)" % (self.in_ds.RasterXSize, self.in_ds.RasterYSize, self.in_ds.RasterCount))

        if not self.in_ds:
            # Note: GDAL prints the ERROR message too
            error("It is not possible to open the input file '%s'." % self.input )
            
        # Read metadata from the input file
        if self.in_ds.RasterCount == 0:
            error( "Input file '%s' has no raster band" % self.input )
            
        if self.in_ds.GetRasterBand(1).GetRasterColorTable():
            # TODO: Process directly paletted dataset by generating VRT in memory
            error( "Please convert this file to RGB/RGBA and run gdal2tiles on the result.",
            """From paletted file you can create RGBA file (temp.vrt) by:
gdal_translate -of vrt -expand rgba %s temp.vrt
then run:
gdal2tiles temp.vrt""" % self.input )

        # Get NODATA value
        self.in_nodata = []
        for i in range(1, self.in_ds.RasterCount+1):
            if self.in_ds.GetRasterBand(i).GetNoDataValue() != None:
                self.in_nodata.append( self.in_ds.GetRasterBand(i).GetNoDataValue() )
        if self.options.srcnodata:
            nds = list(map( float, self.options.srcnodata.split(',')))
            if len(nds) < self.in_ds.RasterCount:
                self.in_nodata = (nds * self.in_ds.RasterCount)[:self.in_ds.RasterCount]
            else:
                self.in_nodata = nds

        if self.options.verbose:
            print("NODATA: %s" % self.in_nodata)

        #
        # Here we should have RGBA input dataset opened in self.in_ds
        #

        if self.options.verbose:
            print("Preprocessed file:", "( %sP x %sL - %s bands)" % (self.in_ds.RasterXSize, self.in_ds.RasterYSize, self.in_ds.RasterCount))

        # Spatial Reference System of the input raster


        self.in_srs,self.in_srs_wkt = def_in_srs(self.options.s_srs,self.input)
        
            #elif self.options.profile != 'raster':
            #   error("There is no spatial reference system info included in the input file.","You should run gdal2tiles with --s_srs EPSG:XXXX or similar.")

        # Spatial Reference System of tiles
        
        self.out_srs=def_out_srs(self.options.profile,self.in_srs)
        
        # Are the reference systems the same? Reproject if necessary.
        
        out_data = init_out_data(self.options.s_srs,self.input,self.options.profile,self.options.verbose,self.in_nodata)
            
        #originX, originY = out_data.out_gt[0], out_data.out_gt[3]
        #pixelSize = out_data.out_gt[1] # = out_data.out_gt[5]
        
        # Test the size of the pixel
        
        # MAPTILER - COMMENTED
        #if out_data.out_gt[1] != (-1 * out_data.out_gt[5]) and self.options.profile != 'raster':
            # TODO: Process corectly coordinates with are have swichted Y axis (display in OpenLayers too)
            #error("Size of the pixel in the output differ for X and Y axes.")
            
        # Report error in case rotation/skew is in geotransform (possible only in 'raster' profile)
        if (out_data.out_gt[2], out_data.out_gt[4]) != (0,0):
            error("Georeference of the raster contains rotation or skew. Such raster is not supported. Please use gdalwarp first.")
            # TODO: Do the warping in this case automaticaly

        # KML test
        isepsg4326 = False
        srs4326 = osr.SpatialReference()
        srs4326.ImportFromEPSG(4326)
        if self.out_srs and srs4326.ExportToProj4() == self.out_srs.ExportToProj4():
            self.kml = True
            isepsg4326 = True
            if self.options.verbose:
                print("KML autotest OK!")

        #
        # Here we expect: pixel is square, no rotation on the raster
        #

        # Output Bounds - coordinates in the output SRS
        tile.ominx = out_data.out_gt[0]
        tile.omaxx = out_data.out_gt[0]+out_data.out_ds.RasterXSize*out_data.out_gt[1]
        tile.omaxy = out_data.out_gt[3]
        tile.ominy = out_data.out_gt[3]-out_data.out_ds.RasterYSize*out_data.out_gt[1]
        # Note: maybe round(x, 14) to avoid the gdal_translate behaviour, when 0 becomes -1e-15

        if self.options.verbose:
            print("Bounds (output srs):", round(tile.ominx, 13), tile.ominy, tile.omaxx, tile.omaxy)

        #
        # Calculating ranges for tiles in different zoom levels
        #
        profile=self.tile_range(out_data,tile,srs4326,isepsg4326)
        
        return out_data, profile

    def tile_range(self,out_data,tile,srs4326,isepsg4326):
        profile = Profile()
        if self.options.profile == 'mercator':
            self.tile_range_mercator(profile,out_data,tile)
        if self.options.profile == 'geodetic':
            self.tile_range_geodetic(profile,out_data,tile)
        if self.options.profile == 'raster':
            self.tile_range_raster(out_data,tile)
            # Function which generates SWNE in LatLong for given tile
            if self.kml and self.in_srs_wkt:
                ct = osr.CoordinateTransformation(self.in_srs, srs4326)
                def rastertileswne(x,y,z):
                    pixelsizex = (2**(self.tmaxz-z) * out_data.out_gt[1]) # X-pixel size in level
                    pixelsizey = (2**(self.tmaxz-z) * out_data.out_gt[1]) # Y-pixel size in level (usually -1*pixelsizex)
                    west = out_data.out_gt[0] + x*TILESIZE*pixelsizex
                    east = west + TILESIZE*pixelsizex
                    south = tile.ominy + y*TILESIZE*pixelsizex
                    north = south + TILESIZE*pixelsizex
                    if not isepsg4326:
                        # Transformation to EPSG:4326 (WGS84 datum)
                        west, south = ct.TransformPoint(west, south)[:2]
                        east, north = ct.TransformPoint(east, north)[:2]
                    return south, west, north, east

                profile.tileswne = rastertileswne
            else:
                profile.tileswne = lambda x, y, z: (0,0,0,0)
        return profile
        
    
    def tile_range_raster(self,out_data,tile):
        log2 = lambda x:math.log10(x) / math.log10(2) # log2 (base 2 logarithm)
        tile.nativezoom = int(max(math.ceil(log2(out_data.out_ds.RasterXSize / float(TILESIZE))), math.ceil(log2(out_data.out_ds.RasterYSize / float(TILESIZE)))))
        if self.options.verbose:
            print "Native zoom of the raster:", tile.nativezoom
    # Get the minimal zoom level (whole raster in one tile)
        if tile.tminz == None:
            tile.tminz = 0
    # Get the maximal zoom level (native resolution of the raster)
        if tile.tmaxz == None:
            tile.tmaxz = tile.nativezoom
    # Generate table with min max tile coordinates for all zoomlevels
        tile.tminmax = list(range(0, tile.tmaxz + 1))
        tile.tsize = list(range(0, tile.tmaxz + 1))
        for tz in range(0, tile.tmaxz + 1):
            tsize = 2.0 ** (tile.nativezoom - tz) * TILESIZE
            tminx, tminy = 0, 0
            tmaxx = int(math.ceil(out_data.out_ds.RasterXSize / tsize)) - 1
            tmaxy = int(math.ceil(out_data.out_ds.RasterYSize / tsize)) - 1
            tile.tsize[tz] = math.ceil(tsize)
            tile.tminmax[tz] = tminx, tminy, tmaxx, tmaxy

    def tile_range_geodetic(self,profile,out_data,tile):
    # from globalmaptiles.py
    # Function which generates SWNE in LatLong for given tile
        profile.tileswne = profile.geodetic.TileLatLonBounds
    # Generate table with min max tile coordinates for all zoomlevels
        tile.tminmax = list(range(0, 32))
        for tz in range(0, 32):
            tminx, tminy = profile.geodetic.LatLonToTile(tile.ominx, tile.ominy, tz)
            tmaxx, tmaxy = profile.geodetic.LatLonToTile(tile.omaxx, tile.omaxy, tz) # crop tiles extending world limits (+-180,+-90)
            tminx, tminy = max(0, tminx), max(0, tminy)
            tmaxx, tmaxy = min(2 ** (tz + 1) - 1, tmaxx), min(2 ** tz - 1, tmaxy)
            tile.tminmax[tz] = tminx, tminy, tmaxx, tmaxy
        
    # TODO: Maps crossing 180E (Alaska?)
    # Get the maximal zoom level (closest possible zoom level up on the resolution of raster)
        if tile.tminz == None:
            tile.tminz = profile.geodetic.ZoomForPixelSize(out_data.out_gt[1] * max(out_data.out_ds.RasterXSize, out_data.out_ds.RasterYSize) / float(TILESIZE))
    # Get the maximal zoom level (closest possible zoom level up on the resolution of raster)
        if tile.tmaxz == None:
            tile.tmaxz = profile.geodetic.ZoomForPixelSize(out_data.out_gt[1])
        if self.options.verbose:
            print "Bounds (latlong):", tile.ominx, tile.ominy, tile.omaxx, tile.omaxy

    def tile_range_mercator(self,profile,out_data,tile):
    # from globalmaptiles.py
    # Function which generates SWNE in LatLong for given tile
        profile.tileswne = profile.mercator.TileLatLonBounds
    # Generate table with min max tile coordinates for all zoomlevels
        tile.tminmax = list(range(0, 32))
        for tz in range(0, 32):
            tminx, tminy = profile.mercator.MetersToTile(tile.ominx, tile.ominy, tz)
            tmaxx, tmaxy = profile.mercator.MetersToTile(tile.omaxx, tile.omaxy, tz) # crop tiles extending world limits (+-180,+-90)
            tminx, tminy = max(0, tminx), max(0, tminy)
            tmaxx, tmaxy = min(2 ** tz - 1, tmaxx), min(2 ** tz - 1, tmaxy)
            tile.tminmax[tz] = tminx, tminy, tmaxx, tmaxy
    # TODO: Maps crossing 180E (Alaska?)
    # Get the minimal zoom level (map covers area equivalent to one tile)
        if tile.tminz == None:
            tile.tminz = profile.mercator.ZoomForPixelSize(out_data.out_gt[1] * max(out_data.out_ds.RasterXSize, out_data.out_ds.RasterYSize) / float(TILESIZE))
    # Get the maximal zoom level (closest possible zoom level up on the resolution of raster)
        if tile.tmaxz == None:
            tile.tmaxz = profile.mercator.ZoomForPixelSize(out_data.out_gt[1])
        if self.options.verbose:
            print "Bounds (latlong):", profile.mercator.MetersToLatLon(tile.ominx, tile.ominy), profile.mercator.MetersToLatLon(tile.omaxx, tile.omaxy)
            print 'MinZoomLevel:', tile.tminz
            print "MaxZoomLevel:", tile.tmaxz, "(", profile.mercator.Resolution(tile.tmaxz), ")"

# =============================================================================

def error(msg, details = "" ):
    """Print an error message and stop the processing"""

    if details:
        sys.stderr.write(msg+ "\n\n" + details)
    else:   
        sys.stderr.write(msg)
    sys.exit(0)
    

def progressbar(complete = 0.0):
    """Print progressbar for float value 0..1"""
    gdal.TermProgress_nocb(complete)


# FONCTIONS METADATAS

def generate_tilemapresource(options,swne,tileext,out_srs,tminz,tmaxz,nativezoom,out_gt):
    """
    Template for tilemapresource.xml. Returns filled string. Expected variables:
      title, north, south, east, west, isepsg4326, projection, publishurl,
      zoompixels, tilesize, tileformat, profile
    """

    args = {}
    args['title'] = options.title
    args['south'], args['west'], args['north'], args['east'] = swne
    args['tilesize'] = TILESIZE
    args['tileformat'] = tileext
    args['publishurl'] = options.url
    args['profile'] = options.profile
    
    if options.profile == 'mercator':
        args['srs'] = "EPSG:900913"
    elif options.profile == 'geodetic':
        args['srs'] = "EPSG:4326"
    elif options.s_srs:
        args['srs'] = options.s_srs
    elif out_srs:
        args['srs'] = out_srs.ExportToWkt()
    else:
        args['srs'] = ""

    s = """<?xml version="1.0" encoding="utf-8"?>
<TileMap version="1.0.0" tilemapservice="http://tms.osgeo.org/1.0.0">
  <Title>%(title)s</Title>
  <Abstract></Abstract>
  <SRS>%(srs)s</SRS>
  <BoundingBox minx="%(south).14f" miny="%(west).14f" maxx="%(north).14f" maxy="%(east).14f"/>
  <Origin x="%(south).14f" y="%(west).14f"/>
  <TileFormat width="%(tilesize)d" height="%(tilesize)d" mime-type="image/%(tileformat)s" extension="%(tileformat)s"/>
  <TileSets profile="%(profile)s">
""" % args
    for z in range(tminz, tmaxz+1):
        if options.profile == 'raster':
            s += """        <TileSet href="%s%d" units-per-pixel="%.14f" order="%d"/>\n""" % (args['publishurl'], z, (2**(nativezoom-z) * out_gt[1]), z)
        elif options.profile == 'mercator':
            s += """        <TileSet href="%s%d" units-per-pixel="%.14f" order="%d"/>\n""" % (args['publishurl'], z, 156543.0339/2**z, z)
        elif options.profile == 'geodetic':
            s += """        <TileSet href="%s%d" units-per-pixel="%.14f" order="%d"/>\n""" % (args['publishurl'], z, 0.703125/2**z, z)
    s += """      </TileSets>
</TileMap>
"""
    return s

def generate_kml(tileext,options,tileswne,tx, ty, tz, children = [], **args):
    """
    Template for the KML. Returns filled string.
    """
    args['tx'], args['ty'], args['tz'] = tx, ty, tz
    args['tileformat'] = tileext
    if 'tilesize' not in args:
        args['tilesize'] = TILESIZE

    if 'minlodpixels' not in args:
        args['minlodpixels'] = int( args['tilesize'] / 2 ) # / 2.56) # default 128
    if 'maxlodpixels' not in args:
        args['maxlodpixels'] = int( args['tilesize'] * 8 ) # 1.7) # default 2048 (used to be -1)
    if children == []:
        args['maxlodpixels'] = -1

    if tx==None:
        tilekml = False
        args['title'] = options.title
    else:
        tilekml = True
        args['title'] = "%d/%d/%d.kml" % (tz, tx, ty)
        args['south'], args['west'], args['north'], args['east'] = tileswne(tx, ty, tz)

    if tx == 0: 
        args['drawOrder'] = 2 * tz + 1 
    elif tx != None: 
        args['drawOrder'] = 2 * tz
    else:
        args['drawOrder'] = 0
        
    url = options.url
    if not url:
        if tilekml:
            url = "../../"
        else:
            url = ""
            
    s = """<?xml version="1.0" encoding="utf-8"?>
<kml xmlns="http://earth.google.com/kml/2.1">
  <Document>
    <Name>%(title)s</Name>
    <Description></Description>
    <Style>
      <ListStyle id="hideChildren">
        <listItemType>checkHideChildren</listItemType>
      </ListStyle>
    </Style>""" % args
    if tilekml:
        s += """
    <Region>
      <Lod>
        <minLodPixels>%(minlodpixels)d</minLodPixels>
        <maxLodPixels>%(maxlodpixels)d</maxLodPixels>
      </Lod>
      <LatLonAltBox>
        <north>%(north).14f</north>
        <south>%(south).14f</south>
        <east>%(east).14f</east>
        <west>%(west).14f</west>
      </LatLonAltBox>
    </Region>
    <GroundOverlay>
      <drawOrder>%(drawOrder)d</drawOrder>
      <Icon>
        <href>%(ty)d.%(tileformat)s</href>
      </Icon>
      <LatLonBox>
        <north>%(north).14f</north>
        <south>%(south).14f</south>
        <east>%(east).14f</east>
        <west>%(west).14f</west>
      </LatLonBox>
    </GroundOverlay>
""" % args  

    for cx, cy, cz in children:
        csouth, cwest, cnorth, ceast = tileswne(cx, cy, cz)
        s += """
    <NetworkLink>
      <name>%d/%d/%d.%s</name>
      <Region>
        <Lod>
          <minLodPixels>%d</minLodPixels>
          <maxLodPixels>-1</maxLodPixels>
        </Lod>
        <LatLonAltBox>
          <north>%.14f</north>
          <south>%.14f</south>
          <east>%.14f</east>
          <west>%.14f</west>
        </LatLonAltBox>
      </Region>
      <Link>
        <href>%s%d/%d/%d.kml</href>
        <viewRefreshMode>onRegion</viewRefreshMode>
        <viewFormat/>
      </Link>
    </NetworkLink>
""" % (cz, cx, cy, args['tileformat'], args['minlodpixels'], cnorth, csouth, ceast, cwest, url, cz, cx, cy)

    s += """      </Document>
</kml>
"""
    return s

def generate_googlemaps(options,swne,tminz,tmaxz,tileext,kml):
    """
    Template for googlemaps.html implementing Overlay of tiles for 'mercator' profile.
    It returns filled string. Expected variables:
    title, googlemapskey, north, south, east, west, minzoom, maxzoom, tilesize, tileformat, publishurl
    """
    args = {}
    args['title'] = options.title
    args['googlemapskey'] = options.googlekey
    args['south'], args['west'], args['north'], args['east'] = swne
    args['minzoom'] = tminz
    args['maxzoom'] = tmaxz
    args['tilesize'] = TILESIZE
    args['tileformat'] = tileext
    args['publishurl'] = options.url
    args['copyright'] = options.copyright

    s = """<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
        <html xmlns="http://www.w3.org/1999/xhtml" xmlns:v="urn:schemas-microsoft-com:vml"> 
          <head>
            <title>%(title)s</title>
            <meta http-equiv="content-type" content="text/html; charset=utf-8"/>
            <meta http-equiv='imagetoolbar' content='no'/>
            <style type="text/css"> v\:* {behavior:url(#default#VML);}
                html, body { overflow: hidden; padding: 0; height: 100%%; width: 100%%; font-family: 'Lucida Grande',Geneva,Arial,Verdana,sans-serif; }
                body { margin: 10px; background: #fff; }
                h1 { margin: 0; padding: 6px; border:0; font-size: 20pt; }
                #header { height: 43px; padding: 0; background-color: #eee; border: 1px solid #888; }
                #subheader { height: 12px; text-align: right; font-size: 10px; color: #555;}
                #map { height: 95%%; border: 1px solid #888; }
            </style>
            <script src='http://maps.google.com/maps?file=api&amp;v=2&amp;key=%(googlemapskey)s' type='text/javascript'></script>
            <script type="text/javascript">
            //<![CDATA[

            /*
             * Constants for given map
             * TODO: read it from tilemapresource.xml
             */

            var mapBounds = new GLatLngBounds(new GLatLng(%(south)s, %(west)s), new GLatLng(%(north)s, %(east)s));
            var mapMinZoom = %(minzoom)s;
            var mapMaxZoom = %(maxzoom)s;

            var opacity = 0.75;
            var map;
            var hybridOverlay;

            /*
             * Create a Custom Opacity GControl
             * http://www.maptiler.org/google-maps-overlay-opacity-control/
             */

            var CTransparencyLENGTH = 58; 
            // maximum width that the knob can move (slide width minus knob width)

            function CTransparencyControl( overlay ) {
                this.overlay = overlay;
                this.opacity = overlay.getTileLayer().getOpacity();
            }
            CTransparencyControl.prototype = new GControl();

            // This function positions the slider to match the specified opacity
            CTransparencyControl.prototype.setSlider = function(pos) {
                var left = Math.round((CTransparencyLENGTH*pos));
                this.slide.left = left;
                this.knob.style.left = left+"px";
                this.knob.style.top = "0px";
            }

            // This function reads the slider and sets the overlay opacity level
            CTransparencyControl.prototype.setOpacity = function() {
                // set the global variable
                opacity = this.slide.left/CTransparencyLENGTH;
                this.map.clearOverlays();
                this.map.addOverlay(this.overlay, { zPriority: 0 });
                if (this.map.getCurrentMapType() == G_HYBRID_MAP) {
                    this.map.addOverlay(hybridOverlay);
                }
            }

            // This gets called by the API when addControl(new CTransparencyControl())
            CTransparencyControl.prototype.initialize = function(map) {
                var that=this;
                this.map = map;

                // Is this MSIE, if so we need to use AlphaImageLoader
                var agent = navigator.userAgent.toLowerCase();
                if ((agent.indexOf("msie") > -1) && (agent.indexOf("opera") < 1)){this.ie = true} else {this.ie = false}

                // create the background graphic as a <div> containing an image
                var container = document.createElement("div");
                container.style.width="70px";
                container.style.height="21px";

                // Handle transparent PNG files in MSIE
                if (this.ie) {
                  var loader = "filter:progid:DXImageTransform.Microsoft.AlphaImageLoader(src='http://www.maptiler.org/img/opacity-slider.png', sizingMethod='crop');";
                  container.innerHTML = '<div style="height:21px; width:70px; ' +loader+ '" ></div>';
                } else {
                  container.innerHTML = '<div style="height:21px; width:70px; background-image: url(http://www.maptiler.org/img/opacity-slider.png)" ></div>';
                }

                // create the knob as a GDraggableObject
                // Handle transparent PNG files in MSIE
                if (this.ie) {
                  var loader = "progid:DXImageTransform.Microsoft.AlphaImageLoader(src='http://www.maptiler.org/img/opacity-slider.png', sizingMethod='crop');";
                  this.knob = document.createElement("div"); 
                  this.knob.style.height="21px";
                  this.knob.style.width="13px";
              this.knob.style.overflow="hidden";
                  this.knob_img = document.createElement("div"); 
                  this.knob_img.style.height="21px";
                  this.knob_img.style.width="83px";
                  this.knob_img.style.filter=loader;
              this.knob_img.style.position="relative";
              this.knob_img.style.left="-70px";
                  this.knob.appendChild(this.knob_img);
                } else {
                  this.knob = document.createElement("div"); 
                  this.knob.style.height="21px";
                  this.knob.style.width="13px";
                  this.knob.style.backgroundImage="url(http://www.maptiler.org/img/opacity-slider.png)";
                  this.knob.style.backgroundPosition="-70px 0px";
                }
                container.appendChild(this.knob);
                this.slide=new GDraggableObject(this.knob, {container:container});
                this.slide.setDraggableCursor('pointer');
                this.slide.setDraggingCursor('pointer');
                this.container = container;

                // attach the control to the map
                map.getContainer().appendChild(container);

                // init slider
                this.setSlider(this.opacity);

                // Listen for the slider being moved and set the opacity
                GEvent.addListener(this.slide, "dragend", function() {that.setOpacity()});
                //GEvent.addListener(this.container, "click", function( x, y ) { alert(x, y) });

                return container;
              }

              // Set the default position for the control
              CTransparencyControl.prototype.getDefaultPosition = function() {
                return new GControlPosition(G_ANCHOR_TOP_RIGHT, new GSize(7, 47));
              }

            /*
             * Full-screen Window Resize
             */

            function getWindowHeight() {
                if (self.innerHeight) return self.innerHeight;
                if (document.documentElement && document.documentElement.clientHeight)
                    return document.documentElement.clientHeight;
                if (document.body) return document.body.clientHeight;
                return 0;
            }

            function getWindowWidth() {
                if (self.innerWidth) return self.innerWidth;
                if (document.documentElement && document.documentElement.clientWidth)
                    return document.documentElement.clientWidth;
                if (document.body) return document.body.clientWidth;
                return 0;
            }

            function resize() {  
                var map = document.getElementById("map");  
                var header = document.getElementById("header");  
                var subheader = document.getElementById("subheader");  
                map.style.height = (getWindowHeight()-80) + "px";
                map.style.width = (getWindowWidth()-20) + "px";
                header.style.width = (getWindowWidth()-20) + "px";
                subheader.style.width = (getWindowWidth()-20) + "px";
                // map.checkResize();
            } 


            /*
             * Main load function:
             */

            function load() {

               if (GBrowserIsCompatible()) {

                  // Bug in the Google Maps: Copyright for Overlay is not correctly displayed
                  var gcr = GMapType.prototype.getCopyrights;
                  GMapType.prototype.getCopyrights = function(bounds,zoom) {
                      return ["%(copyright)s"].concat(gcr.call(this,bounds,zoom));
                  }

                  map = new GMap2( document.getElementById("map"), { backgroundColor: '#fff' } );

                  map.addMapType(G_PHYSICAL_MAP);
                  map.setMapType(G_PHYSICAL_MAP);

                  map.setCenter( mapBounds.getCenter(), map.getBoundsZoomLevel( mapBounds ));

                  hybridOverlay = new GTileLayerOverlay( G_HYBRID_MAP.getTileLayers()[1] );
                  GEvent.addListener(map, "maptypechanged", function() {
                    if (map.getCurrentMapType() == G_HYBRID_MAP) {
                        map.addOverlay(hybridOverlay);
                    } else {
                       map.removeOverlay(hybridOverlay);
                    }
                  } );

                  var tilelayer = new GTileLayer(GCopyrightCollection(''), mapMinZoom, mapMaxZoom);
                  var mercator = new GMercatorProjection(mapMaxZoom+1);
                  tilelayer.getTileUrl = function(tile,zoom) {
                      if ((zoom < mapMinZoom) || (zoom > mapMaxZoom)) {
                          return "http://www.maptiler.org/img/none.png";
                      } 
                      var ymax = 1 << zoom;
                      var y = ymax - tile.y -1;
                      var tileBounds = new GLatLngBounds(
                          mercator.fromPixelToLatLng( new GPoint( (tile.x)*256, (tile.y+1)*256 ) , zoom ),
                          mercator.fromPixelToLatLng( new GPoint( (tile.x+1)*256, (tile.y)*256 ) , zoom )
                      );
                      if (mapBounds.intersects(tileBounds)) {
                          return zoom+"/"+tile.x+"/"+y+"."""+tileext+"""";
                      } else {
                          return "http://www.maptiler.org/img/none.png";
                      }
                  }
                  // IE 7-: support for PNG alpha channel
                  // Unfortunately, the opacity for whole overlay is then not changeable, either or...
                  tilelayer.isPng = function() { return true;};
                  tilelayer.getOpacity = function() { return opacity; }

                  overlay = new GTileLayerOverlay( tilelayer );
                  map.addOverlay(overlay);

                  map.addControl(new GLargeMapControl());
                  map.addControl(new GHierarchicalMapTypeControl());
                  map.addControl(new CTransparencyControl( overlay ));
    """ % args
    if kml:
        s += """
                  map.addMapType(G_SATELLITE_3D_MAP);
                  map.getEarthInstance(getEarthInstanceCB);
    """
    s += """

                  map.enableContinuousZoom();
                  map.enableScrollWheelZoom();

                  map.setMapType(G_HYBRID_MAP);
               }
               resize();
            }
    """
    if kml:
        s += """
            function getEarthInstanceCB(object) {
               var ge = object;

               if (ge) {
                   var url = document.location.toString();
                   url = url.substr(0,url.lastIndexOf('/'))+'/doc.kml';
                   var link = ge.createLink("");
                   if ("%(publishurl)s") { link.setHref("%(publishurl)s/doc.kml") }
                   else { link.setHref(url) };
                   var networkLink = ge.createNetworkLink("");
                   networkLink.setName("TMS Map Overlay");
                   networkLink.setFlyToView(true);  
                   networkLink.setLink(link);
                   ge.getFeatures().appendChild(networkLink);
               } else {
                   // alert("You should open a KML in Google Earth");
                   // add div with the link to generated KML... - maybe JavaScript redirect to the URL of KML?
               }
            }
    """ % args
    s += """
            onresize=function(){ resize(); };

            //]]>
            </script>
          </head>
          <body onload="load()">
              <div id="header"><h1>%(title)s</h1></div>
              <div id="subheader">Generated by <a href="http://www.maptiler.org/">MapTiler</a>/<a href="http://www.klokan.cz/projects/gdal2tiles/">GDAL2Tiles</a>, Copyright &copy; 2008 <a href="http://www.klokan.cz/">Klokan Petr Pridal</a>,  <a href="http://www.gdal.org/">GDAL</a> &amp; <a href="http://www.osgeo.org/">OSGeo</a> <a href="http://code.google.com/soc/">GSoC</a>
        <!-- PLEASE, LET THIS NOTE ABOUT AUTHOR AND PROJECT SOMEWHERE ON YOUR WEBSITE, OR AT LEAST IN THE COMMENT IN HTML. THANK YOU -->
              </div>
               <div id="map"></div>
          </body>
        </html>
    """ % args

    return s

def generate_openlayers(options,swne,tminz,tmaxz,tileext,nativezoom,out_gt):
    """
    Template for openlayers.html implementing overlay of available Spherical Mercator layers.

    It returns filled string. Expected variables:
    title, googlemapskey, yahooappid, north, south, east, west, minzoom, maxzoom, tilesize, tileformat, publishurl
    """

    args = {}
    args['title'] = options.title
    args['googlemapskey'] = options.googlekey
    args['yahooappid'] = options.yahookey
    args['south'], args['west'], args['north'], args['east'] = swne
    args['minzoom'] = tminz
    args['maxzoom'] = tmaxz
    args['tilesize'] = TILESIZE
    args['tileformat'] = tileext
    args['publishurl'] = options.url
    args['copyright'] = options.copyright
    if options.profile == 'raster':
        args['rasterzoomlevels'] = tmaxz+1
        args['rastermaxresolution'] = 2**(nativezoom) * out_gt[1]

    s = """<!DOCTYPE html PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" "http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd">
    <html xmlns="http://www.w3.org/1999/xhtml"
      <head>
        <title>%(title)s</title>
        <meta http-equiv='imagetoolbar' content='no'/>
        <style type="text/css"> v\:* {behavior:url(#default#VML);}
            html, body { overflow: hidden; padding: 0; height: 100%%; width: 100%%; font-family: 'Lucida Grande',Geneva,Arial,Verdana,sans-serif; }
            body { margin: 10px; background: #fff; }
            h1 { margin: 0; padding: 6px; border:0; font-size: 20pt; }
            #header { height: 43px; padding: 0; background-color: #eee; border: 1px solid #888; }
            #subheader { height: 12px; text-align: right; font-size: 10px; color: #555;}
            #map { height: 95%%; border: 1px solid #888; }
        </style>""" % args
    
    if options.profile == 'mercator':
        s += """
        <script src='http://dev.virtualearth.net/mapcontrol/mapcontrol.ashx?v=6.1'></script>
        <script src='http://maps.google.com/maps?file=api&amp;v=2&amp;key=%(googlemapskey)s' type='text/javascript'></script>
        <script src="http://api.maps.yahoo.com/ajaxymap?v=3.0&amp;appid=%(yahooappid)s"></script>""" % args

    s += """
        <script src="http://www.openlayers.org/api/2.7/OpenLayers.js" type="text/javascript"></script>
        <script type="text/javascript">
            var map;
            var mapBounds = new OpenLayers.Bounds( %(west)s, %(south)s, %(east)s, %(north)s);
            var mapMinZoom = %(minzoom)s;
            var mapMaxZoom = %(maxzoom)s;

            // avoid pink tiles
            OpenLayers.IMAGE_RELOAD_ATTEMPTS = 3;
            OpenLayers.Util.onImageLoadErrorColor = "transparent";

            function init(){""" % args

    if options.profile == 'mercator':
        s += """
            var options = {
                controls: [],
                projection: new OpenLayers.Projection("EPSG:900913"),
                displayProjection: new OpenLayers.Projection("EPSG:4326"),
                units: "m",
                maxResolution: 156543.0339,
                maxExtent: new OpenLayers.Bounds(-20037508, -20037508, 20037508, 20037508.34)
                };
            map = new OpenLayers.Map('map', options);

            // create Google Mercator layers
            var gmap = new OpenLayers.Layer.Google("Google Streets",
                { sphericalMercator: true, numZoomLevels: 20} );
            var gsat = new OpenLayers.Layer.Google("Google Satellite",
                {type: G_SATELLITE_MAP, sphericalMercator: true, numZoomLevels: 20} );
            var ghyb = new OpenLayers.Layer.Google("Google Hybrid",
                {type: G_HYBRID_MAP, sphericalMercator: true, numZoomLevels: 20});
            var gter = new OpenLayers.Layer.Google("Google Terrain",
                {type: G_PHYSICAL_MAP, sphericalMercator: true, numZoomLevels: 20 });

            // create Virtual Earth layers
            OpenLayers.Layer.VirtualEarth.prototype.MAX_ZOOM_LEVEL=19;
            OpenLayers.Layer.VirtualEarth.prototype.RESOLUTIONS=OpenLayers.Layer.Google.prototype.RESOLUTIONS
            var veroad = new OpenLayers.Layer.VirtualEarth("Virtual Earth Roads",
                {'type': VEMapStyle.Road, 'sphericalMercator': true, numZoomLevels: 20});
            var veaer = new OpenLayers.Layer.VirtualEarth("Virtual Earth Aerial",
                {'type': VEMapStyle.Aerial, 'sphericalMercator': true, numZoomLevels: 20 });
            var vehyb = new OpenLayers.Layer.VirtualEarth("Virtual Earth Hybrid",
                {'type': VEMapStyle.Hybrid, 'sphericalMercator': true});

            // create Yahoo layer
            var yahoo = new OpenLayers.Layer.Yahoo("Yahoo Street",
                {'sphericalMercator': true});
            var yahoosat = new OpenLayers.Layer.Yahoo("Yahoo Satellite",
                {'type': YAHOO_MAP_SAT, 'sphericalMercator': true});
            var yahoohyb = new OpenLayers.Layer.Yahoo("Yahoo Hybrid",
                {'type': YAHOO_MAP_HYB, 'sphericalMercator': true});

            // create OSM/OAM layer
            var osm = new OpenLayers.Layer.TMS( "OpenStreetMap",
                "http://tile.openstreetmap.org/",
                { type: '"""+tileext+"""', getURL: osm_getTileURL, displayOutsideMaxExtent: true, attribution: '<a href="http://www.openstreetmap.org/">OpenStreetMap</a>'} );
            var oam = new OpenLayers.Layer.TMS( "OpenAerialMap",
                "http://tile.openaerialmap.org/tiles/1.0.0/openaerialmap-900913/",
                { type: '"""+tileext+"""', getURL: osm_getTileURL } );

            // create TMS Overlay layer
            var tmsoverlay = new OpenLayers.Layer.TMS( "TMS Overlay", "",
                {   // url: '', serviceVersion: '.', layername: '.',
                    type: '"""+tileext+"""', getURL: overlay_getTileURL, alpha: true, 
                    isBaseLayer: false
                });
            if (OpenLayers.Util.alphaHack() == false) { tmsoverlay.setOpacity(0.7); }

            map.addLayers([gmap, gsat, ghyb, gter, veroad, veaer, vehyb,
                           yahoo, yahoosat, yahoohyb, osm, oam,
                           tmsoverlay]);

            var switcherControl = new OpenLayers.Control.LayerSwitcher();
            map.addControl(switcherControl);
            switcherControl.maximizeControl();

            map.zoomToExtent( mapBounds.transform(map.displayProjection, map.projection ) );
        """ % args

    elif options.profile == 'geodetic':
        s += """
            var options = {
                controls: [],
                projection: new OpenLayers.Projection("EPSG:4326"),
                maxResolution: 0.703125,
                maxExtent: new OpenLayers.Bounds(-180, -90, 180, 90)
                };
            map = new OpenLayers.Map('map', options);

            layer = new OpenLayers.Layer.WMS( "Blue Marble",
                    "http://labs.metacarta.com/wms-c/Basic.py?", {layers: 'satellite' } );
            map.addLayer(layer);
            wms = new OpenLayers.Layer.WMS( "VMap0",
                    "http://labs.metacarta.com/wms-c/Basic.py?", {layers: 'basic', format: 'image/png' } );
            map.addLayer(wms);
            
            var tmsoverlay = new OpenLayers.Layer.TMS( "TMS Overlay", "",
                {
                    serviceVersion: '.', layername: '.', alpha: true,
                    type: '"""+tileext+"""', getURL: overlay_getTileURL,
                    isBaseLayer: false
                });
            map.addLayer(tmsoverlay);
            if (OpenLayers.Util.alphaHack() == false) { tmsoverlay.setOpacity(0.7); }

            var switcherControl = new OpenLayers.Control.LayerSwitcher();
            map.addControl(switcherControl);
            switcherControl.maximizeControl();

            map.zoomToExtent( mapBounds );
        """
        
    elif options.profile == 'raster':
        s += """
            var options = {
                controls: [],
                maxExtent: new OpenLayers.Bounds(  %(west)s, %(south)s, %(east)s, %(north)s ),
                maxResolution: %(rastermaxresolution)f,
                numZoomLevels: %(rasterzoomlevels)d
                };
            map = new OpenLayers.Map('map', options);

            var layer = new OpenLayers.Layer.TMS( "TMS Layer","",
                {  url: '', serviceVersion: '.', layername: '.', alpha: true,
                    type: '"""+tileext+"""', getURL: overlay_getTileURL 
                });
            map.addLayer(layer);
            map.zoomToExtent( mapBounds );  
    """ % args


    s += """
            map.addControl(new OpenLayers.Control.PanZoomBar());
            map.addControl(new OpenLayers.Control.MousePosition());
            map.addControl(new OpenLayers.Control.MouseDefaults());
            map.addControl(new OpenLayers.Control.KeyboardDefaults());
        }
        """ % args
    
    if options.profile == 'mercator':
        s += """
        function osm_getTileURL(bounds) {
            var res = this.map.getResolution();
            var x = Math.round((bounds.left - this.maxExtent.left) / (res * this.tileSize.w));
            var y = Math.round((this.maxExtent.top - bounds.top) / (res * this.tileSize.h));
            var z = this.map.getZoom();
            var limit = Math.pow(2, z);

            if (y < 0 || y >= limit) {
                return "http://www.maptiler.org/img/none.png";
            } else {
                x = ((x %% limit) + limit) %% limit;
                return this.url + z + "/" + x + "/" + y + "." + this.type;
            }
        }

        function overlay_getTileURL(bounds) {
            var res = this.map.getResolution();
            var x = Math.round((bounds.left - this.maxExtent.left) / (res * this.tileSize.w));
            var y = Math.round((bounds.bottom - this.tileOrigin.lat) / (res * this.tileSize.h));
            var z = this.map.getZoom();
            if (this.map.baseLayer.name == 'Virtual Earth Roads' || this.map.baseLayer.name == 'Virtual Earth Aerial' || this.map.baseLayer.name == 'Virtual Earth Hybrid') {
               z = z + 1;
            }
            if (mapBounds.intersectsBounds( bounds ) && z >= mapMinZoom && z <= mapMaxZoom ) {
               //console.log( this.url + z + "/" + x + "/" + y + "." + this.type);
               return this.url + z + "/" + x + "/" + y + "." + this.type;
            } else {
               return "http://www.maptiler.org/img/none.png";
            }
        }       
        """ % args
        
    elif options.profile == 'geodetic':
        s += """
        function overlay_getTileURL(bounds) {
            bounds = this.adjustBounds(bounds);
            var res = this.map.getResolution();
            var x = Math.round((bounds.left - this.tileOrigin.lon) / (res * this.tileSize.w));
            var y = Math.round((bounds.bottom - this.tileOrigin.lat) / (res * this.tileSize.h));
            var z = this.map.getZoom();
            var path = this.serviceVersion + "/" + this.layername + "/" + z + "/" + x + "/" + y + "." + this.type;
            var url = this.url;
            if (mapBounds.intersectsBounds( bounds ) && z >= mapMinZoom && z <= mapMaxZoom) {
               // console.log( this.url + z + "/" + x + "/" + y + "." + this.type);
               return this.url + z + "/" + x + "/" + y + "." + this.type;
            } else {
               return "http://www.maptiler.org/img/none.png";
            }
        }
        """ % args
        
    elif options.profile == 'raster':
        s += """
        function overlay_getTileURL(bounds) {
            var res = this.map.getResolution();
            var x = Math.round((bounds.left - this.maxExtent.left) / (res * this.tileSize.w));
            var y = Math.round((bounds.bottom - this.maxExtent.bottom) / (res * this.tileSize.h));
            var z = this.map.getZoom();
            if (x >= 0 && y >= 0) {
                return this.url + z + "/" + x + "/" + y + "." + this.type;              
            } else {
                return "http://www.maptiler.org/img/none.png";
            }
        }
        """ % args
    
    s += """
       function getWindowHeight() {
            if (self.innerHeight) return self.innerHeight;
            if (document.documentElement && document.documentElement.clientHeight)
                return document.documentElement.clientHeight;
            if (document.body) return document.body.clientHeight;
                return 0;
        }

        function getWindowWidth() {
            if (self.innerWidth) return self.innerWidth;
            if (document.documentElement && document.documentElement.clientWidth)
                return document.documentElement.clientWidth;
            if (document.body) return document.body.clientWidth;
                return 0;
        }

        function resize() {  
            var map = document.getElementById("map");  
            var header = document.getElementById("header");  
            var subheader = document.getElementById("subheader");  
            map.style.height = (getWindowHeight()-80) + "px";
            map.style.width = (getWindowWidth()-20) + "px";
            header.style.width = (getWindowWidth()-20) + "px";
            subheader.style.width = (getWindowWidth()-20) + "px";
            if (map.updateSize) { map.updateSize(); };
        } 

        onresize=function(){ resize(); };

        </script>
      </head>
      <body onload="init()">
        <div id="header"><h1>%(title)s</h1></div>
        <div id="subheader">Generated by <a href="http://www.maptiler.org/">MapTiler</a>/<a href="http://www.klokan.cz/projects/gdal2tiles/">GDAL2Tiles</a>, Copyright &copy; 2008 <a href="http://www.klokan.cz/">Klokan Petr Pridal</a>,  <a href="http://www.gdal.org/">GDAL</a> &amp; <a href="http://www.osgeo.org/">OSGeo</a> <a href="http://code.google.com/soc/">GSoC</a>
        <!-- PLEASE, LET THIS NOTE ABOUT AUTHOR AND PROJECT SOMEWHERE ON YOUR WEBSITE, OR AT LEAST IN THE COMMENT IN HTML. THANK YOU -->
        </div>
        <div id="map"></div>
        <script type="text/javascript" >resize()</script>
      </body>
    </html>""" % args

    return s


def generate_metadata(config,profile,tile,out_data):
    """Generation of main metadata files and HTML viewers (metadata related to particular tiles are generated during the tile processing)."""
    
    if not os.path.exists(config.output):
        os.makedirs(config.output)

    if config.options.profile == 'mercator':
        
        south, west = profile.mercator.MetersToLatLon(tile.ominx, tile.ominy)
        north, east = profile.mercator.MetersToLatLon(tile.omaxx, tile.omaxy)
        south, west = max(-85.05112878, south), max(-180.0, west)
        north, east = min(85.05112878, north), min(180.0, east)
        swne = (south, west, north, east)

        # Generate googlemaps.html
        if config.options.webviewer in ('all','google') and config.options.profile == 'mercator':
            if not config.options.resume or not os.path.exists(os.path.join(config.output, 'googlemaps.html')):
                f = open(os.path.join(config.output, 'googlemaps.html'), 'w')
                f.write( generate_googlemaps(config.options,swne,tile.tminz,tile.tmaxz,config.tileext,config.kml) )
                f.close()

        # Generate openlayers.html
        if config.options.webviewer in ('all','openlayers'):
            if not config.options.resume or not os.path.exists(os.path.join(config.output, 'openlayers.html')):
                f = open(os.path.join(config.output, 'openlayers.html'), 'w')
                f.write( generate_openlayers(config.options,swne,tile.tminz,tile.tmaxz,config.tileext,tile.nativezoom,out_data.out_gt) )
                f.close()

    elif config.options.profile == 'geodetic':
        
        west, south = tile.ominx, tile.ominy
        east, north = tile.omaxx, tile.omaxy
        south, west = max(-90.0, south), max(-180.0, west)
        north, east = min(90.0, north), min(180.0, east)
        swne = (south, west, north, east)
        
        # Generate openlayers.html
        if config.options.webviewer in ('all','openlayers'):
            if not config.options.resume or not os.path.exists(os.path.join(config.output, 'openlayers.html')):
                f = open(os.path.join(config.output, 'openlayers.html'), 'w')
                f.write( generate_openlayers(config.options,swne,tile.tminz,tile.tmaxz,config.tileext,tile.nativezoom,out_data.out_gt) )
                f.close()           

    elif config.options.profile == 'raster':
        
        west, south = tile.ominx, tile.ominy
        east, north = tile.omaxx, tile.omaxy

        swne = (south, west, north, east)
        
        # Generate openlayers.html
        if config.options.webviewer in ('all','openlayers'):
            if not config.options.resume or not os.path.exists(os.path.join(config.output, 'openlayers.html')):
                f = open(os.path.join(config.output, 'openlayers.html'), 'w')
                f.write( generate_openlayers(config.options,swne,tile.tminz,tile.tmaxz,config.tileext,tile.nativezoom,out_data.out_gt) )
                f.close()           


    # Generate tilemapresource.xml.
    if not config.options.resume or not os.path.exists(os.path.join(config.output, 'tilemapresource.xml')):
        f = open(os.path.join(config.output, 'tilemapresource.xml'), 'w')
        f.write( generate_tilemapresource(config.options,swne,config.tileext,config.out_srs,tile.tminz,tile.tmaxz,tile.nativezoom,out_data.out_gt))
        f.close()

    if config.kml:
        # TODO: Maybe problem for not automatically generated tminz
        # The root KML should contain links to all tiles in the tminz level
        children = []
        xmin, ymin, xmax, ymax = tile.tminmax[tile.tminz]
        for x in range(xmin, xmax+1):
            for y in range(ymin, ymax+1):
                children.append( [ x, y, tile.tminz ] ) 
        # Generate Root KML
        if config.kml:
            if not config.options.resume or not os.path.exists(os.path.join(config.output, 'doc.kml')):
                f = open(os.path.join(config.output, 'doc.kml'), 'w')
                f.write( generate_kml(config.tileext,config.options,profile.tileswne, None, None, None, children) )
                f.close()


#FONCTIONS TILES
def geo_query(ds, ulx, uly, lrx, lry, querysize = 0):
    """For given dataset and query in cartographic coordinates
    returns parameters for ReadRaster() in raster coordinates and
    x/y shifts (for border tiles). If the querysize is not given, the
    extent is returned in the native resolution of dataset ds."""

    geotran = ds.GetGeoTransform()
    rx= int((ulx - geotran[0]) / geotran[1] + 0.001)
    ry= int((uly - geotran[3]) / geotran[5] + 0.001)
    rxsize= int((lrx - ulx) / geotran[1] + 0.5)
    rysize= int((lry - uly) / geotran[5] + 0.5)

    if not querysize:
        wxsize, wysize = rxsize, rysize
    else:
        wxsize, wysize = querysize, querysize

    # Coordinates should not go out of the bounds of the raster
    wx = 0
    if rx < 0:
        rxshift = abs(rx)
        wx = int( wxsize * (float(rxshift) / rxsize) )
        wxsize = wxsize - wx
        rxsize = rxsize - int( rxsize * (float(rxshift) / rxsize) )
        rx = 0
    if rx+rxsize > ds.RasterXSize:
        wxsize = int( wxsize * (float(ds.RasterXSize - rx) / rxsize) )
        rxsize = ds.RasterXSize - rx

    wy = 0
    if ry < 0:
        ryshift = abs(ry)
        wy = int( wysize * (float(ryshift) / rysize) )
        wysize = wysize - wy
        rysize = rysize - int( rysize * (float(ryshift) / rysize) )
        ry = 0
    if ry+rysize > ds.RasterYSize:
        wysize = int( wysize * (float(ds.RasterYSize - ry) / rysize) )
        rysize = ds.RasterYSize - ry

    return (rx, ry, rxsize, rysize), (wx, wy, wxsize, wysize)

def scale_query_to_tile(options,tiledriver,resampling,dsquery, dstile, tilefilename=''):
    """Scales down query dataset to the tile dataset"""

    querysize = dsquery.RasterXSize
    tilesize = dstile.RasterXSize
    tilebands = dstile.RasterCount

    if options.resampling == 'average':

        # Function: gdal.RegenerateOverview()
        for i in range(1,tilebands+1):
            # Black border around NODATA
            #if i != 4:
            #   dsquery.GetRasterBand(i).SetNoDataValue(0)
            res = gdal.RegenerateOverview( dsquery.GetRasterBand(i),
                dstile.GetRasterBand(i), 'average' )
            if res != 0:
                error("RegenerateOverview() failed on %s, error %d" % (tilefilename, res))

    elif options.resampling == 'antialias':

        # Scaling by PIL (Python Imaging Library) - improved Lanczos
        array = numpy.zeros((querysize, querysize, tilebands), numpy.uint8)
        for i in range(tilebands):
            array[:,:,i] = gdalarray.BandReadAsArray(dsquery.GetRasterBand(i+1), 0, 0, querysize, querysize)
        im = Image.fromarray(array, 'RGBA') # Always four bands
        im1 = im.resize((tilesize,tilesize), Image.ANTIALIAS)
        if os.path.exists(tilefilename):
            im0 = Image.open(tilefilename)
            im1 = Image.composite(im1, im0, im1) 
        im1.save(tilefilename,tiledriver)
        
    else:

        # Other algorithms are implemented by gdal.ReprojectImage().
        dsquery.SetGeoTransform( (0.0, tilesize / float(querysize), 0.0, 0.0, 0.0, tilesize / float(querysize)) )
        dstile.SetGeoTransform( (0.0, 1.0, 0.0, 0.0, 0.0, 1.0) )

        res = gdal.ReprojectImage(dsquery, dstile, None, None, resampling)
        if res != 0:
            error("ReprojectImage() failed on %s, error %d" % (tilefilename, res))

def tile_bounds(tmaxx, tmaxy, ds, querysize, tz, ty, tx,options,mercator,geodetic,tsize,out_ds,nativezoom):
    "Considering the profile option, initialize the bounds for reading and writing the tile"
    if options.profile == 'mercator':
        # Tile bounds in EPSG:900913
        b = mercator.TileBounds(tx, ty, tz)
    elif options.profile == 'geodetic':
        b = geodetic.TileBounds(tx, ty, tz)

    if options.profile in ('mercator', 'geodetic'):
        rb, wb = geo_query(ds, b[0], b[3], b[2], b[1])
        nativesize = wb[0] + wb[2]
        # Pixel size in the raster covering query geo extent
        if options.verbose:
            print "\tNative Extent (querysize", nativesize, "): ", rb, wb
        # Tile bounds in raster coordinates for ReadRaster query
        rb, wb = geo_query(ds, b[0], b[3], b[2], b[1], querysize=querysize)
    else:
        tsize = int(tsize[tz]) # tilesize in raster coordinates for actual zoom
        xsize = out_ds.RasterXSize # size of the raster in pixels
        ysize = out_ds.RasterYSize
        if tz >= nativezoom:
            querysize = TILESIZE # int(2**(nativezoom-tz) * tilesize)
        rx = (tx) * tsize
        rxsize = 0
        if tx == tmaxx:
            rxsize = xsize % tsize
        if rxsize == 0:
            rxsize = tsize
        rysize = 0
        if ty == tmaxy:
            rysize = ysize % tsize
        if rysize == 0:
            rysize = tsize
        ry = ysize - (ty * tsize) - rysize
        wx, wy = 0, 0
        wxsize, wysize = int(rxsize / float(tsize) * TILESIZE), int(rysize / float(tsize) * TILESIZE)
        if wysize != TILESIZE:
            wy = TILESIZE - wysize
        rb = (rx,ry,rxsize,rysize)
        wb = (wx,wy,wxsize,wysize)
    return rb, wb, querysize

def generate_base_tile(ds, tilebands, querysize, tz, ty, tx, tilefilename, rb, wb, options, tiledriver, resampling, mem_drv, out_drv, out_data, tile):
    "Create 1 base tile"
    if not (options.resume and os.path.exists(tilefilename)):
        # Query is in 'nearest neighbour' but can be bigger in then the tilesize
        # We scale down the query to the tilesize by supplied algorithm.
        
        # Tile dataset in memory
        rx, ry, rxsize, rysize = rb
        wx, wy, wxsize, wysize = wb
        dstile = mem_drv.Create('', TILESIZE, TILESIZE, tilebands)
        data = ds.ReadRaster(rx, ry, rxsize, rysize, wxsize, wysize, band_list=list(range(1, out_data.dataBandsCount + 1)))
        alpha = out_data.alphaband.ReadRaster(rx, ry, rxsize, rysize, wxsize, wysize)
        if TILESIZE == querysize:
            # Use the ReadRaster result directly in tiles ('nearest neighbour' query)
            dstile.WriteRaster(wx, wy, wxsize, wysize, data, band_list=list(range(1, out_data.dataBandsCount + 1)))
            dstile.WriteRaster(wx, wy, wxsize, wysize, alpha, band_list=[tilebands])
        else:
            dsquery = mem_drv.Create('', querysize, querysize, tilebands)
            # TODO: fill the null value in case a tile without alpha is produced (now only png tiles are supported)
            #for i in range(1, tilebands+1):
            #   dsquery.GetRasterBand(1).Fill(tilenodata)
            dsquery.WriteRaster(wx, wy, wxsize, wysize, data, band_list=list(range(1, out_data.dataBandsCount + 1)))
            dsquery.WriteRaster(wx, wy, wxsize, wysize, alpha, band_list=[tilebands])
            scale_query_to_tile(options,tiledriver,resampling,dsquery, dstile, tilefilename)
            del dsquery
        # Big ReadRaster query in memory scaled to the tilesize - all but 'near' algo
        # Note: For source drivers based on WaveLet compression (JPEG2000, ECW, MrSID)
        # the ReadRaster function returns high-quality raster (not ugly nearest neighbour)
        # TODO: Use directly 'near' for WaveLet files
        del data
        if options.resampling != 'antialias':
            # Write a copy of tile to png/jpg
            out_drv.CreateCopy(tilefilename, dstile, strict=0)
        del dstile
    
    
def process_base_tile_job(out_data,tiledriver, options, resampling, tilebands, querysize, tz, ty, tx, tile,output,tileext,tmaxx,tmaxy,mercator,geodetic):
    "Create the drivers and the tile bounds, then generate the base tile"
    out_drv,mem_drv = init_drv(tiledriver)
    ds=out_data.out_ds
    tilefilename = os.path.join(output, str(tz), str(tx), "%s.%s" % (ty, tileext))
    
    rb,wb,querysize= tile_bounds(tmaxx, tmaxy, ds, querysize, tz, ty, tx,options,mercator,geodetic,tile.tsize,out_data.out_ds,tile.nativezoom)

    generate_base_tile(ds, tilebands, querysize, tz, ty, tx, tilefilename, rb, wb, options, tiledriver, resampling, mem_drv, out_drv, out_data, tile)


def generate_base_tiles(config,profile,tile,out_data):
    """Generation of the base tiles (the lowest in the pyramid) directly from the input raster"""
    
    print("Generating Base Tiles:")
    
    if config.options.verbose:
        print('')
        print("Tiles generated from the max zoom level:")
        print("----------------------------------------")
        print('')


    # Set the bounds
    tminx, tminy, tmaxx, tmaxy = tile.tminmax[tile.tmaxz]

    # Just the center tile
    #tminx = tminx+ (tmaxx - tminx)/2
    #tminy = tminy+ (tmaxy - tminy)/2
    #tmaxx = tminx
    #tmaxy = tminy
    
    tilebands = out_data.dataBandsCount + 1
    querysize = config.querysize_c
    
    if config.options.verbose:
        print("dataBandsCount: ", out_data.dataBandsCount)
        print("tilebands: ", tilebands)
    
    tcount = (1+abs(tmaxx-tminx)) * (1+abs(tmaxy-tminy))
    
    tz = tile.tmaxz
    
    # Create directories for the tiles
    for tx in range(tminx,tmaxx+1):
        tilefolder = os.path.join(config.output, str(tz), str(tx))
        if not os.path.exists(tilefolder):
            os.makedirs(tilefolder)
    
    # To know were we stand in the creation of the tiles, we use a shared counter ti
    # that all processes can modify.
    ti=multiprocessing.Value('f',0.0)
    
    multiprocess = MultiProcess(config.nb_process)
    multiprocess.processBase(config.options.s_srs, config.input, config.options.profile, config.options.verbose, config.in_nodata,ti)
    for ty in range(tmaxy, tminy-1, -1): #range(tminy, tmaxy+1):
               
        lx=range(tminx,tmaxx+1)
        multiprocess.sendJob((config.tiledriver, config.options, config.resampling, tilebands, querysize, tz, ty, lx, tile,config.output,config.tileext,tmaxx,tmaxy,profile.mercator,profile.geodetic))                    
    multiprocess.finish(ti,tcount)
    
    if config.kml:
        for ty in range(tmaxy,tminy-1,-1):
            for tx in range(tminx, tmaxx+1):
            # Create a KML file for this tile.
                kmlfilename = os.path.join(config.output, str(tz), str(tx), '%d.kml' % ty)
                if not config.options.resume or not os.path.exists(kmlfilename):
                    f = open(kmlfilename, 'w')
                    f.write(generate_kml(config.tileext,config.options,profile.tileswne,tx, ty, tz))
                    f.close()


def read_tile(tz, output, tileext, y, x):
    "Read a tile"
    dsquerytile = gdal.Open(os.path.join(output, str(tz + 1), str(x), "%s.%s" % (y, tileext)), gdal.GA_ReadOnly)
    tile_r = dsquerytile.ReadRaster(0, 0, TILESIZE, TILESIZE)
    return tile_r

def init_children(tile,tz,ty,tx):
    "Initialize a table for kml generation"
    children = []
    for y in range(2 * ty, 2 * ty + 2):
        for x in range(2 * tx, 2 * tx + 2):
            minx, miny, maxx, maxy = tile.tminmax[tz + 1]
            if x >= minx and x <= maxx and y >= miny and y <= maxy:
                children.append([x, y, tz + 1])
    return children


def generate_overview_tile(tilebands, tz, ty, tx, tilefilename,mem_drv,out_drv,output,options,tiledriver,resampling,tileext,tile):
    "Generate 1 tile from 4 underlying tiles"
    if not (options.resume and os.path.exists(tilefilename)):
        dsquery = mem_drv.Create('', 2 * TILESIZE, 2 * TILESIZE, tilebands)
        # TODO: fill the null value
        #for i in range(1, tilebands+1):
        #   dsquery.GetRasterBand(1).Fill(tilenodata)
        dstile = mem_drv.Create('', TILESIZE, TILESIZE, tilebands)
        # TODO: Implement more clever walking on the tiles with cache functionality
        # probably walk should start with reading of four tiles from top left corner
        # Hilbert curve...
        
    # Read the tiles and write them to query window
        for y in range(2 * ty, 2 * ty + 2):
            for x in range(2 * tx, 2 * tx + 2):
                minx, miny, maxx, maxy = tile.tminmax[tz + 1]
                if x >= minx and x <= maxx and y >= miny and y <= maxy:
                    tile_r = read_tile(tz, output, tileext, y, x)
                    if (ty == 0 and y == 1) or (ty != 0 and (y % (2 * ty)) != 0):
                        tileposy = 0
                    else:
                        tileposy = TILESIZE
                    if tx:
                        tileposx = x % (2 * tx) * TILESIZE
                    elif tx == 0 and x == 1:
                        tileposx = TILESIZE
                    else:
                        tileposx = 0
                    dsquery.WriteRaster(tileposx, tileposy, TILESIZE, TILESIZE,tile_r, band_list=list(range(1, tilebands + 1)))
        
        scale_query_to_tile(options,tiledriver,resampling,dsquery, dstile, tilefilename)
        # Write a copy of tile to png/jpg
        if options.resampling != 'antialias':
            # Write a copy of tile to png/jpg
            out_drv.CreateCopy(tilefilename, dstile, strict=0)
        if options.verbose:
            print "\tbuild from zoom", tz + 1, " tiles:", 2 * tx, 2 * ty, 2 * tx + 1, 2 * ty, 2 * tx, 2 * ty + 1, 2 * tx + 1, 2 * ty + 1 # Create a KML file for this tile.


def process_overview_tile_job(output,options,tiledriver,resampling,tileext, tile, tilebands, tz, ty, tx):
    "Create the drivers and generate the overview tile"
    out_drv,mem_drv = init_drv(tiledriver)
    tilefilename = os.path.join(output, str(tz), str(tx), "%s.%s" % (ty, tileext))
    generate_overview_tile(tilebands, tz, ty, tx, tilefilename,mem_drv,out_drv,output,options,tiledriver,resampling,tileext,tile)

def generate_overview_tiles(config,profile,tile,out_data):
    """Generation of the overview tiles (higher in the pyramid) based on existing tiles"""
    
    print("Generating Overview Tiles:")
    
    tilebands = out_data.dataBandsCount + 1
    
    # Create directories for the tiles
    for tz in range(tile.tmaxz-1,tile.tminz-1,-1):
        tminx, tminy, tmaxx, tmaxy = tile.tminmax[tz]
        for tx in range(tminx,tmaxx+1):
                tilefolder = os.path.join(config.output, str(tz), str(tx))
                if not os.path.exists(tilefolder):
                    os.makedirs(tilefolder)    
    
    # Usage of existing tiles: from 4 underlying tiles generate one as overview.
    
    tcount = 0
    for tz in range(tile.tmaxz-1, tile.tminz-1, -1):
        tminx, tminy, tmaxx, tmaxy = tile.tminmax[tz]
        tcount += (1+abs(tmaxx-tminx)) * (1+abs(tmaxy-tminy))

    # To know were we stand in the creation of the tiles, we use a shared counter ti
    # that all processes can modify.
    ti = multiprocessing.Value('f',0.0)
    
    # querysize = TILESIZE * 2
    for tz in range(tile.tmaxz-1, tile.tminz-1, -1):
        multiprocess=MultiProcess(config.nb_process)
        multiprocess.processOverview(ti)
        tminx, tminy, tmaxx, tmaxy = tile.tminmax[tz]
        for ty in range(tmaxy, tminy-1, -1): #range(tminy, tmaxy+1):
            
            # If there is less jobs to calculate than processes, we send one tile per job.
            # Else, we send a list of tiles in one job.
            if tmaxy-tminy<config.nb_process:            
                for tx in range(tminx, tmaxx+1):                    
                    multiprocess.sendJob((config.output, config.options, config.tiledriver, config.resampling, config.tileext, tile, tilebands, tz, ty, tx))
            else:
                lx=range(tminx,tmaxx+1)
                multiprocess.sendJob((config.output, config.options, config.tiledriver, config.resampling, config.tileext, tile, tilebands, tz, ty, lx))
        multiprocess.finish(ti,tcount)
        if config.kml:
            for ty in range(tmaxy,tminy-1,-1):
                for tx in range(tminx, tmaxx+1):
                    children=init_children(tile,tz, ty, tx)
                    f = open(os.path.join(config.output, '%d\\%d\\%d.kml' % (tz, tx, ty)), 'w')
                    f.write(generate_kml(config.tileext,config.options,profile.tileswne,tx, ty, tz, children))
                    f.close()



def process(config,tile):
    """The main processing function, runs all the main steps of processing"""
    
    # Opening and preprocessing of the input file
    out_data,profile=config.open_input(tile)
    
    # Generation of main metadata files and HTML viewers
    generate_metadata(config,profile,tile,out_data)

    # Generation of the lowest tiles
    generate_base_tiles(config,profile,tile,out_data)
    
    # Generation of the overview tiles (higher in the pyramid)
    generate_overview_tiles(config,profile,tile,out_data)


# =============================================================================
# =============================================================================

if __name__=='__main__':
    argv = gdal.GeneralCmdLineProcessor( sys.argv )
    if argv:
        c1 = Configuration(argv[1:])
        tile=c1.create_tile()
        process(c1,tile)