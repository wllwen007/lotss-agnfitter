import pylab as plt
import numpy as np
import os


filter_res_files = ['filter.bootes_mbrown_2014a.res',
                    'FIR_filters.res',
                    'EN1_filters.res',
                    'Lockman-SWIRE_filters.res']


for filter_res_file in filter_res_files:
    print filter_res_file
    with open(filter_res_file) as f:
        alllines = f.readlines()
        #while read:
        filternames = []
        filterlc = []
        filterstart = []
        filterstop= []
        for iline, line in enumerate(alllines):
            C = line.strip()
            #if 'lambda' in C:
            if ('filter' in C) or ('res' in C)or ('txt' in C):
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
        
    #figdir = filter_res_file.replace('.res','_plots')
    filtdir = filter_res_file.replace('.res','_filters')
    #if not os.path.isdir(figdir):
        #os.mkdir(figdir)
    if not os.path.isdir(filtdir):
        os.mkdir(filtdir)
        
    plt.figure(figsize=(12,6))
    for fi in range(len(filternames)):
        plt.semilogx(filterlambda[fi], filtertrans[fi])
        lamc = np.nanmean(filterlambda[fi])
        plt.xlim(1e3,1e7)
        plt.ylim(-0.01, 1.01)
        #plt.title(filternames[fi])
        #plt.savefig('filter_plots/'+filternames[fi]+'.png')
        plt.savefig(filtdir+'.png')
        
        np.savez(filtdir+'/'+filternames[fi]+'.npz', l=filterlambda[fi], r=filtertrans[fi])
        np.savetxt(filtdir+'/'+filternames[fi]+'.filter', np.array([filterlambda[fi], filtertrans[fi]]).T)
        #plt.text(lamc, 0.1, str(fi))

    
        
            

