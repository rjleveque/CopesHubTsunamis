
from pylab import *
from clawpack.visclaw import animation_tools

plotdir='plots_GH_SeismicVsNosub_46-95'

images = animation_tools.read_images(plotdir=plotdir,
                                     fname_pattern='dtopo_transect_GH_SeismicVsNosub_*.png')

anim = animation_tools.animate_images(images, figsize=(12,7))

file_name = f'{plotdir}/{plotdir}.mp4'
animation_tools.make_mp4(anim, file_name=file_name)
print('Created ',file_name)
