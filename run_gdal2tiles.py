import os
from time import clock

profile_list = ('mercator','geodetic','raster')
resampling_list = ('average','near','antialias')

chemin = 'D:\\Stage\\gdal2tiles\\test_image\\'

img = 'D:\Stage\gdal2tiles\image_test_miniature\image.vrt'


for p in profile_list:
    for r in resampling_list:
        if (p=='raster' and r=='antialias')!=1:
            print p+" "+r
            path=chemin+p+'_'+r
            os.chdir(path)
            cmd="gdal2tiles.py -p %s -r %s %s" % (p, r, img)
            t0=clock()
            os.system(cmd)
            t1=clock()
            #f = os.open(p+'_'+r,'a+')
            #f.write(t1-t0+'\n')
            print t1-t0