import os
from time import clock

BRANCH = 'multiprocess(process=4)'
CFG_NUM_ITERATIONS = 1
THIS_DIR = os.path.dirname( __file__ )

profile_list = ('mercator','geodetic','raster')
resampling_list = ('average','near','antialias')

chemin = os.path.join( THIS_DIR, "test_image" )

img = os.path.join( THIS_DIR, "image_test_miniature", "image.vrt" )

stats = os.path.join( THIS_DIR, "stats.csv" )


with open(stats, 'a') as f:
    f.write(BRANCH)
    for p in profile_list:
        for r in resampling_list:
            name=p+'-'+r
            t_list=[]
            if (p,r) != ('raster', 'antialias'):
                print name
                for i in range (CFG_NUM_ITERATIONS):
                    path=os.path.join(chemin, p+'_'+r)
                    if not os.path.exists(path):
                        os.makedirs(path)
                    os.chdir(path)
                    cmd=r"%s\gdal2tiles.py -p %s -r %s %s" % (THIS_DIR,p, r, img)
                    t0=clock()
                    os.system(cmd)
                    t1=clock()
                    t=t1-t0
                    t_list[len(t_list):]=[t]
                    print t
                f.write(' '+str(min(t_list)))
            else:
                print name+' is not available.'
    f.write('\n')