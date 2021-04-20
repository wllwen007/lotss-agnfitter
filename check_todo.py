import numpy as np
import glob

#cd EN1_v1.0/agnfitter_out/

t1 = glob.glob('*/done*')
t1 = np.array([int(ti[:ti.find('/done')]) for ti in t1])
t1.sort()

todo = np.array([ti for ti in  np.arange(t1.max()) if ti not in t1])

t2 = glob.glob('*/na*')
t2 = np.array([int(ti[:ti.find('/na')]) for ti in t2])
t2.sort()

todo = np.array([ti for ti in todo if ti not in t2])


t3 = glob.glob('*/fail*')
t3 = np.array([int(ti[:ti.find('/fail')]) for ti in t3])
t3.sort()

todo = np.array([ti for ti in todo if ti not in t3])

print(len(t1), 'done')
print(len(t2), 'na')
print(len(t3), 'failed')
print(len(todo), 'not done/na/failed')

s = ''
for t in todo: s+=str(t)+','


#s = ''
#for t in todo: 
    #if (int(str(t)[-2:]) > 95) and  (int(str(t)[-2:]) < 100) :
        #s+=str(t)+','

