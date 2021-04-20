import sys
import os
from astropy.table import Table, Column
import astropy.units as au
import numpy as np
import utils.plot_util as pp

#field = 'Bootes'
#field = 'LH'
#field = 'EN1'
field = sys.argv[1]

if field=='all':
    fields = ['Bootes','EN1','LH']
else:
    fields = [field]
    

# filter files should be in Angstrom
c=    2.997e8
Angstrom = 1e10

filterpathroot = '/data2/wwilliams/projects/lofar_surveys/deep/science_ready_catalogs/filter_information/'
outfilterpath = 'agnfitter/filters/'

for field in fields:
    
    if field == 'Bootes':
        datfile = 'agnfitter_cols_bootes.csv'
        filterpath = 'filter.bootes_mbrown_2014a_filters/'
    elif field == 'EN1':
        datfile = 'agnfitter_cols_en1.csv'
        filterpath = 'EN1_filters_filters/'
    elif field == 'LH':
        datfile = 'agnfitter_cols_lockman.csv'
        filterpath = 'Lockman-SWIRE_filters_filters/'
        
    filterpathFIR = 'FIR_filters_filters/'
    outfilterpath = 'agnfitter/filters/'+field+'/'
        
    filterdat = Table.read('agnfitter/filters/'+datfile)
    filterdat['wavelength'].dtype=float

    
    for f in filterdat:
        ffile = filterpathroot+filterpath+f['filtername']+'.filter'
        
        if not os.path.isfile(ffile):
            print(ffile)
            
            ffile = filterpathroot+filterpathFIR+f['filtername']+'.filter'
            if not os.path.isfile(ffile):
                print(ffile)
                print("I don't know what to do")
                sys.exit()
        
        band_lambda, band_factor =  np.loadtxt(ffile, usecols=(0,1),unpack= True)
        os.system('cp '+ffile+' '+outfilterpath)
        print(('cp '+ffile+' '+outfilterpath))

        central_lamb = np.sum(band_lambda*band_factor)/np.sum(band_factor)
        central_nu = float(np.log10((Angstrom*c)/central_lamb))
        f['wavelength'] = central_lamb
        
    # deal with the duplicate bands - by shifting the bandpass by -0.1, 0, 1 angstrom (incase of 3)
    filterdat.add_column(Column(data=filterdat['wavelength'], name='af_wavelength'))
    filterdat.add_column(Column(data=filterdat['filtername'], name='af_filtername'))
        
    w,nw = np.unique(filterdat['wavelength'], return_counts=True)
    for wi, nwi in zip(w,nw):
        if nwi > 1:
            print((wi, nwi))
            fwi = np.where(filterdat['wavelength'] == wi)[0]
            
            
            for i in range(len(fwi)):
                filterdat['af_filtername'][fwi[i]] = 'AF'+str(i)+'_'+filterdat['filtername'][fwi[i]]
                
                
                ffile = filterpathroot+filterpath+filterdat['filtername'][fwi[i]]+'.filter'
                affile = filterpathroot+filterpath+filterdat['af_filtername'][fwi[i]]+'.filter'
            
      
                band_lambda, band_factor =  np.loadtxt(ffile, usecols=(0,1),unpack= True)
                
                with open(affile, 'w') as f:
                    for bl, bf in zip(band_lambda, band_factor):
                        f.write('%f\t%f\n'%(bl-0.001+0.001*i, bf))
                        
                band_lambda, band_factor =  np.loadtxt(affile, usecols=(0,1),unpack= True)
            
            
                os.system('cp '+affile+' '+outfilterpath)
                print(('cp '+affile+' '+outfilterpath))

                central_lamb = np.sum(band_lambda*band_factor)/np.sum(band_factor)
                central_nu = float(np.log10((Angstrom*c)/central_lamb))
                filterdat['af_wavelength'][fwi[i]] = central_lamb
            
        else:
            print('{:f} {:d} is ok'.format(wi, nwi))
            
    
    w,nw = np.unique(filterdat['af_wavelength'], return_counts=True)
    for wi, nwi in zip(w,nw):
        if nwi > 1:
            print ('Still have duplicates!')
        
    filterdat.write('agnfitter/filters/'+datfile.replace('.csv','_wl.csv'),overwrite=True)
