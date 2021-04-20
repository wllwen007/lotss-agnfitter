from astropy.table import Table, Column
import numpy as np
import utils.plot_util as pp

#bands = ['FUV', 'NUV', 'u', 'Bw', 'R', 'I', 'z_Subaru', 'y', 'J', 'H', 'Ks', 'ch1', 'ch2', 'ch3', 'ch4', '24', 'F250', 'F350', 'F500']
# bootes - leave out for agnfitter 'z' and 'K' - can't do overlapping filters

'''
### Bootes
----------

FLAG_CLEAN == 1
FLAG_DEEP != 0

I_fluxerr > 0
ch2_fluxerr > 0

### ELAIS-N1
------------

FLAG_CLEAN == 1

i_fluxerr > 0
ch2_swire_fluxerr > 0

or

(Note, the FLAG_OVERLAP == 7 being applied for LR matching will supercede these cuts)

### Lockman-Hole
----------------

FLAG_CLEAN == 1

r_fluxerr > 0
ch2_swire_fluxerr > 0
'''




bandfiles = {'Bw': 'Bw.filter',
'H': 'H.filter',
'I': 'I.filter',
'ch1': 'IRAC_ch1_total_response.filter',
'ch2': 'IRAC_ch2_total_response.filter',
'ch3': 'IRAC_ch3_total_response.filter',
'ch4': 'IRAC_ch4_total_response.filter',
'J': 'J.filter',
'Ks': 'Ks.filter',
#'': 'PACS_100mu.txt',
#'': 'PACS_160mu.txt',
'R': 'R.filter',
'F250': 'SPIRE_250mu.txt',
'F350': 'SPIRE_350mu.txt',
'F500': 'SPIRE_500mu.txt',
'u': 'U.filter',
'y': 'Y.filter',
#'FUV': 'galexFUV.filter',
#'NUV': 'galexNUV.filter',
'F24': 'mips24.filter',
'z_Subaru': 'subaru_z.filter' }
bandwaves = {}

f,ax = pp.paper_single_ax()

with open('LOFAR-Deep-BootesFilters-select.lst','w') as f:
    
    for band in bandfiles.keys():
        bandfile = bandfiles[band]
        #tt = f.readlines()
        print band, bandfile, 
        band_lambda, band_factor =  np.loadtxt('LOFAR-Deep-BootesFilters/'+bandfile, usecols=(0,1),unpack= True)
        
        band_factor = band_factor/np.max(band_factor)
        bandwl = np.average(band_lambda, weights=band_factor) 
        print bandwl*1e-10 /1e-6
        
        bandwaves[band] = bandwl
        
        ax.plot(band_lambda*1e-10 /1e-6, band_factor)
        ax.vlines(bandwl, 0, 1, 'k')
        
        f.write('LOFAR-Deep-BootesFilters/'+bandfile + '\n')


#tt = Table.read('data/ML_RUN1_srl_ionly_AllData.fits')
tt = Table.read('data/bootes_cat.fits')
tt = Table.read('/beegfs/lofar/deepfields/science_ready_catalogs/Bootes_opt_spitzer_merged_vac_opt3as_irac4as_all_hpx_public.fits')



tt.add_column(Column(tt['z1_median'], 'z_best'))
spec = np.isfinite(tt['z_spec']) & (tt['z_spec'] >0)
tt['z_best'][spec]  = tt['z_spec'][spec]

keep_list = ['ID','z_best','z1_median','z_spec']
for band in bandfiles.keys():
    keep_list.append(band+'_wl')
    
    bandwl = bandwaves[band]
    tt.add_column( Column(bandwl*np.ones(len(tt)), band+'_wl' ))    

    # add fluxes
    keep_list.append(band+'_f')
    if band not in ['F24', 'F250', 'F350', 'F500']:
        tt.add_column( Column(tt[band+'_flux_4'], band+'_f' )) 
    else:
        tt.add_column( Column(tt[band], band+'_f' ))    
        
    # scale to Jy
    if band not in ['F250', 'F350', 'F500']:   
        tt[band+'_f'] *= 1e-6  ## to Jy
    else:
        tt[band+'_f'] *= 1e-3  ## to Jy
        
    # add flux errors
    keep_list.append(band+'_e')
    if band not in ['F24', 'F250', 'F350', 'F500']:
        tt.add_column( Column(tt[band+'_fluxerr_4'], band+'_e' ))  
    else:
        tt.add_column( Column(tt['e_'+band], band+'_e' ))   
        
    # check for -99 values
    tt[band+'_f'][tt[band+'_e']==-99] = np.nan
    tt[band+'_e'][tt[band+'_e']==-99] = np.nan
    
    # scale to Jy
    if band not in ['F250', 'F350', 'F500']:     
        tt[band+'_e'] *= 1e-6  ## to Jy
    else:
        tt[band+'_e'] *= 1e-3  ## to Jy
        
    # add flux error of 10% of flux in quadrature to flux error
    tt[band+'_e'] = np.sqrt( (tt[band+'_e'])**2. + (0.1*tt[band+'_f'])**2.)
    
    tt[band+'_f'].fill_value = np.nan
    tt[band+'_e'].fill_value = np.nan
    
ttout = tt.copy()
ttout.keep_columns(keep_list)

#ttout = ttout[np.isfinite(ttout['z_best']) & (ttout['z_best'] > 0)]

ttout.write('data/bootes_cat_agnfitter_in.fits', overwrite=True)
