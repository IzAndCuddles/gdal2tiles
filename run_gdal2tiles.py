import os
from time import clock

CFG_NUM_ITERATIONS = 1
THIS_DIR = os.path.dirname( __file__ )

profile_list = ('mercator','geodetic','raster')
resampling_list = ('average','near','antialias')

chemin = os.path.join( THIS_DIR, "test_image" )

img = os.path.join( THIS_DIR, "image_test_miniature", "image.vrt" )

stats = os.path.join( THIS_DIR, "stats.csv" )
with open(stats, 'w') as f:

    for p in profile_list:
        for r in resampling_list:
            name=p+'-'+r
            if (p,r) != ('raster', 'antialias'):
                f.write(name)
                print name
                for i in range (CFG_NUM_ITERATIONS):
                    path=os.path.join(chemin, p+'_'+r)
                    if not os.path.isdir(path):
                        os.makedirs(path)
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
