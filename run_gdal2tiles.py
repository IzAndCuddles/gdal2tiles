import os
from time import clock

profile_list = ('mercator','geodetic','raster')
resampling_list = ('average','near','antialias')

chemin = 'D:\\Stage\\gdal2tiles\\test_image\\'

img = 'D:\Stage\gdal2tiles\image_test_miniature\image.vrt'

f = open('D:\\Stage\\gdal2tiles\\stats.txt','w')

for p in profile_list:
    for r in resampling_list:
        name=p+'-'+r
        if not (p=='raster' and r=='antialias'):
            f.write(name)
            f = open('D:\\Stage\\gdal2tiles\\stats.txt','a')
            print name
            for i in range (10):
                path=chemin+p+'_'+r
                os.chdir(path)
                cmd="gdal2tiles.py -p %s -r %s %s" % (p, r, img)
                t0=clock()
                os.system(cmd)
                t1=clock()
                t=t1-t0
                f.write(' '+str(t))
                print t
            f.write('\r\n')
        else:
            print name+' is not available.'