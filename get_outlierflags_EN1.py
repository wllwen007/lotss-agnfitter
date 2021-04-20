from astropy.table import Table
import numpy as np
import os

ft = '/beegfs/lofar/wwilliams/lofar_surveys/deep/agnfitter/EN1_final_cross_match_catalogue-v0.5_mergenew.fits'
ftoutlier = '/beegfs/lofar/deepfields/science_ready_catalogs/EN1_opt_spitzer_merged_vac_opt3as_irac4as_all_hpx_public.fits'
ftmatch = 'temp_match.fits'



os.system('stilts tmatch2 in1={in1} in2={in2} out={out} values1=ID values2=id join=all1 find=best suffix2=outlier matcher=exact'.format(in1=ft, in2=ftoutlier, out=ftmatch))


t = Table.read(ft)
toutlier = Table.read(ftmatch)


for col in t.colnames:
    if 'flux_corr' in col:
        col1 = col.replace('_corr','')
        Table([t[col], toutlier[col1]])
        print col1, col, np.sum(toutlier[col1] == -90),np.sum(np.isnan(t[col])),
        t[col][toutlier[col1] == -90] = np.nan
        print np.sum(np.isnan(t[col]))

t.write('/beegfs/lofar/wwilliams/lofar_surveys/deep/agnfitter/EN1_final_cross_match_catalogue-v0.5_mergenew_outlier.fits',overwrite=True)


