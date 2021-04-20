import pylab as pl
import numpy as np
import os

filter_res_file = 'ndwfs_brown.FILTER.RES'

if not os.path.exists('filter_plots'):
    os.mkdir('filter_plots')


with open(filter_res_file) as f:
    alllines = f.readlines()
    #while read:
    filternames = []
    filterlc = []
    filterstart = []
    filterstop= []
    for iline, line in enumerate(alllines):
        C = line.strip()
        if 'lambda' in C:
            D = C.split()
            nlines = np.int(D[0])
            print C, nlines, iline
            if '/' in D[1]:
                D[1] = D[1].replace('/','_')
            filternames.append(D[1])
            filterlc.append(D[1])
            filterstart.append(iline+1)
            filterstop.append(iline+1+nlines)
    filterlambda = []
    filtertrans = []
    for fi in range(len(filternames)):
        lam = []
        trans = []
        for iline in range(filterstart[fi], filterstop[fi]):
            C = alllines[iline]
            D = C.strip().split()
            if len(D) == 3:
                lam.append(np.float(D[1]))
                trans.append(np.float(D[2]))
            elif len(D) == 2:
                lam.append(np.float(D[0]))
                trans.append(np.float(D[1]))
                print iline-filterstart[fi]+1, D[0], D[1]
        lam = np.array(lam)
        trans= np.array(trans)
        filterlambda.append(lam)
        filtertrans.append(trans)
        
        
pl.figure()
for fi in range(len(filternames)):
    pl.semilogx(filterlambda[fi], filtertrans[fi])
    lamc = np.nanmean(filterlambda[fi])
    pl.xlim(1e3,1e6)
    pl.ylim(-0.01, 1.01)
    #pl.title(filternames[fi])
    #pl.savefig('filter_plots/'+filternames[fi]+'.png')
    pl.savefig('filter_plots/'+'filters'+'.png')
    
    np.savez('filter_plots/'+filternames[fi]+'.npz', l=filterlambda[fi], r=filtertrans[fi])
    #pl.text(lamc, 0.1, str(fi))

    
        
            

