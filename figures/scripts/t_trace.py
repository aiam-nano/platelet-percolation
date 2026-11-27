import time, numpy as np, math
from figgen import *
t=time.time()
cfg=make_cfg('cyl',L=16,seed=3,A=30)
c0,lab,hist,wr=run_trace(cfg)
print('count',c0,time.time()-t)
nc=c0/16**3; print('phi',nc*math.pi*2/30)
# compare to production
print('prod', run_once(cfg['mode'],cfg['pos'],cfg['rot'],cfg['rad'],cfg['hth'],16.0,cfg['brad']))
print((lab==wr).sum(), hist[-1], len(np.unique(lab)))
