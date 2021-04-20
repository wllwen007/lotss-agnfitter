import numpy as np
from astropy.table import Table, Column
import utils.plot_util as pp

path = '/beegfs/lofar/deepfields/'
cat_fnames = {'Bootes': 'science_ready_catalogs/Bootes_opt_spitzer_merged_vac_opt3as_irac4as_all_hpx_public.fits',
              'LH': 'science_ready_catalogs/LH_opt_spitzer_merged_vac_opt3as_irac4as_all_hpx_public.fits',
              'EN1': 'science_ready_catalogs/EN1_opt_spitzer_merged_vac_opt3as_irac4as_all_hpx_public.fits'}

cats = ['Bootes', 'LH', 'EN1']
#cats = ['LH']

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

def select_good(cat, name):
    if name == 'Bootes':
        select = (cat['FLAG_CLEAN'] == 1) & (cat['FLAG_DEEP'] != 0) & (cat['I_fluxerr'] > 0) & (cat['ch2_fluxerr'] > 0 )
    elif name == 'LH':
        select = (cat['FLAG_CLEAN'] == 1) & (cat['r_fluxerr'] > 0) & (cat['ch2_swire_fluxerr'] > 0)
    elif name == 'EN1':
        select = (cat['FLAG_CLEAN'] == 1) & (cat['i_fluxerr'] > 0) & (cat['ch2_swire_fluxerr'] > 0)
    else:
        print ('not implemented')
        sys.exit(1)
        
    return cat[select]

for catname in cats:
    mcat = Table.read(path+cat_fnames[catname])
    
    #mcat.add_column(Column(mcat['Z_SPEC'],'z_spec', dtype='float64'))
    #mcat.rename_column('Z_SPEC','z_spec')
    #mcat[(mcat['z_spec'] == -99)]['z_spec'] *= np.nan   
    #mcat['z_spec'][np.where(mcat['z_spec']==-99)[0]] = np.nan
    
    #mcat = col_num_to_nan(mcat, 'z_spec', -99)
    
    #add_zbest(mcat)
    zspec = np.array(mcat['Z_SPEC'], dtype='float64')
    mcat.remove_column('Z_SPEC')
    
    zspec[zspec == -99] = np.nan
    mcat.add_column(Column(zspec,'z_spec'))
    
    # hack to get rid of those damn -99
    z1med = np.array(mcat['z1_median'], dtype='float64')
    mcat.remove_column('z1_median')
    z1med[z1med == -99] = np.nan
    mcat.add_column(Column(z1med,'z1_median'))
    
    zbest = zspec.copy()
    zphi = np.isnan(zbest)
    zbest[zphi] = z1med[zphi]
    mcat.add_column(Column(zbest,'z_best'))
    
    
    #add_zerr(mcat)
    z1max = np.array(mcat['z1_max'], dtype='float64')
    z1min = np.array(mcat['z1_min'], dtype='float64')
    z1max[z1max == -99] = np.nan
    z1min[z1min == -99] = np.nan
    zerr =  0.5*(z1max-z1min)/(1+z1med)
    mcat.add_column(Column(zerr,'z1_err'))
    
    print('{c}: {l}'.format(c=catname, l=len(mcat)))
    
    #mcat_sel = make_selection(mcat)
    #mcat_sel = mcat[np.isfinite(mcat['z_best'])]
    goodzsel = np.isfinite(zspec) | (zerr <0.1)
    mcat_sel = mcat[goodzsel]
    
    
    print('{c}: {l} with zspec or good zphot'.format(c=catname, l=len(mcat_sel)))
    print('{c}: {l} with zspec'.format(c=catname, l=np.sum(np.isfinite(mcat['z_spec']))))
    
    mcat_sel = select_good(mcat_sel, catname)
    print('{c}: {l} good selection from zsel'.format(c=catname, l=len(mcat_sel)))
    
    mcat = select_good(mcat, catname)
    print('{c}: {l} good selection from all'.format(c=catname, l=len(mcat)))

    mcat_sel.write('science_ready_catalogs/'+catname+'_select.fits', overwrite=True)

