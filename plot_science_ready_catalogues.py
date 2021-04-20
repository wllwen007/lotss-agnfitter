import numpy as np
from astropy.table import Table, Column
import utils.plot_util as pp

path = '/beegfs/lofar/deepfields/'
cat_fnames = {'Bootes': 'science_ready_catalogs/Bootes_opt_spitzer_merged_vac_opt3as_irac4as_all_hpx_public.fits',
              'LH': 'science_ready_catalogs/LH_opt_spitzer_merged_vac_opt3as_irac4as_all_hpx_public.fits',
              'EN1': 'science_ready_catalogs/EN1_opt_spitzer_merged_vac_opt3as_irac4as_all_hpx_public.fits'}

cats = ['Bootes', 'LH', 'EN1']
cats = ['EN1']

def add_zbest(cat):
    zbest = cat['z_spec']
    zph = np.isnan(zbest)
    zbest[zph] = cat['z1_median'][zph]
    cat.add_column(Column(zbest,'z_best'))
    return cat

def add_zerr(cat):
    zerr =  0.5*(cat['z1_max']-cat['z1_min'])/(1+cat['z1_median'])
    cat.add_column(Column(zerr,'z1_err'))
    return 

def make_selection(cat):
    cat = cat[np.isfinite(cat['z_best'])]
    goodsel = np.isfinite(cat['z_spec']) | (cat['z1_err'] <0.1)
    cat = cat[goodsel]
    return cat

def col_num_to_nan(cat, col, num):
    cat[(cat[col] == num)][col] *= np.nan   
    return cat 

def save_cat(cat, name):
    savename  = 'science_ready_catalogs/{catname}_opt_zselect.fits'.format(catname=name)
    cat.write(savename, overwrite=True)
    return
    

def select_good(cat, name):
    if name == 'Bootes':
        select = (cat['FLAG_CLEAN'] == 1) & (cat['FLAG_DEEP'] != 0) & (cat['I_fluxerr'] > 0) & (cat['ch2_fluxerr'] > 0 )
    elif name == 'LH':
        select = (cat['FLAG_CLEAN'] == 1) & (cat['i_fluxerr'] > 0) & (cat['ch2_swire_fluxerr'] > 0)
    elif name == 'EN1':
        select = (cat['FLAG_CLEAN'] == 1) & (cat['r_fluxerr'] > 0) & (cat['ch2_swire_fluxerr'] > 0)
    else:
        print ('not implemented')
        sys.exit(1)
        
    return cat[select]

for catname in cats:
    mcat = Table.read(path+cat_fnames[catname])
    
    zz = mcat['Z_SPEC']
    zz[zz==-99] = np.nan
    mcat.remove_column('Z_SPEC')
    mcat.add_column(Column(zz,'z_spec', dtype='float64'))
    
    #add_zbest(mcat)
    zbest = np.array(mcat['z_spec'])
    zph = np.isnan(zbest)
    zbest[zph] = mcat['z1_median'][zph]
    mcat.add_column(Column(zbest,'z_best', dtype='float64'))
    
    #add_zerr(mcat)
    zerr =  0.5*(mcat['z1_max']-mcat['z1_min'])/(1+mcat['z1_median'])
    mcat.add_column(Column(zerr,'z1_err'))
    
    print('{c}: {l}'.format(c=catname, l=len(mcat)))
    
    #mcat_sel = make_selection(mcat)
    mcat_sel = mcat[np.isfinite(mcat['z_best'])]
    goodsel = np.isfinite(mcat_sel['z_spec']) | (mcat_sel['z1_err'] <0.1)
    mcat_sel = mcat_sel[goodsel]
    
    
    print('{c}: {l} with zspec or good zphot'.format(c=catname, l=len(mcat_sel)))
    print('{c}: {l} with zspec'.format(c=catname, l=len(mcat[np.isfinite(mcat['z_spec'])])))
    
    mcat_sel = select_good(mcat_sel, catname)
    print('{c}: {l} good selection'.format(c=catname, l=len(mcat_sel)))
    
    mcat = select_good(mcat, catname)
    print('{c}: {l} good selection from all'.format(c=catname, l=len(mcat)))
    
    save_cat(mcat_sel,catname)

