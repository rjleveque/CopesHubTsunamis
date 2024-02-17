from pylab import *
from clawpack.geoclaw import topotools

topo = topotools.Topography('/Users/rjl/topo/topofiles/etopo1_-163_-122_38_63.asc',3)

ygs = arange(48.75,50.9,0.15)
topovi = topo.crop([-130,-124,48.5,52])
xgs = zeros(ygs.shape)
for k,yg in enumerate(ygs):
    j = where(topovi.y>yg)[0].min()
    Bj = topovi.Z[j,:]
    i1 = where(Bj>0)[0].min()
    i = where(Bj[:i1]<=-100)[0].max()
    xgs[k] = xg = topovi.x[i]
    Bg = topovi.Z[j,i]
    print('x = %.4f, y = %.4f, B = %.1f' % (xg,yg,Bg))

fig,ax = subplots(figsize=(7,7))
topovi.plot(axes=ax)
plot(xgs,ygs,'ko')

fname = 'gaugesVI.txt'
with open(fname,'w') as f:
    f.write('# 100m depth gauges along Vancouver Island\n')
    for k in range(len(xgs)):
        f.write('%12.4f  %10.4f\n' % (xgs[k],ygs[k]))
    
