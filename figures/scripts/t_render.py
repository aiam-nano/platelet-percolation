import numpy as np, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from figgen import *; from render import *
cfg=make_cfg('cyl',L=16,seed=3,A=30)
c0,lab,hist,wr=run_trace(cfg)
y0,y1=6.0,10.0; ids=slab_ids(cfg,c0,y0,y1); print(len(ids))
fig,ax=plt.subplots(figsize=(5,5),dpi=150)
allp,r,u,v=draw(ax,cfg,lab,wr,ids,view=(-28,24))
b=frame_box(ax,16,y0,y1,(-28,24)); 
a=np.vstack([allp,b]); ax.set_xlim(a[:,0].min()-.3,a[:,0].max()+.3); ax.set_ylim(a[:,1].min()-.3,a[:,1].max()+.3); ax.set_aspect('equal'); ax.axis('off')
fig.savefig('t_render.png',bbox_inches='tight')
