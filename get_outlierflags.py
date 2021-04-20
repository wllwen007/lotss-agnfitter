from astropy.table import Table
import numpy as np
import os
import sys
field = sys.argv[1]

if field=='all':
    fields = ['Bootes','EN1','LH']
else:
    fields = [field]

for field in fields:    

    if field == 'EN1':
        # mergenew comes from  combine_new_fir
        ft = '/beegfs/lofar/wwilliams/lofar_surveys/deep/agnfitter/en1_final_cross_match_catalogue-v0.7.fits'
        ftoutlier = '/beegfs/lofar/deepfields/science_ready_catalogs/EN1_opt_spitzer_merged_vac_opt3as_irac4as_all_hpx_public.fits'
        
        ftout = '/beegfs/lofar/wwilliams/lofar_surveys/deep/agnfitter/EN1_final_cross_match_catalogue-v0.7_outlier.fits'
        
    elif field == 'LH':
        ft = '/beegfs/lofar/wwilliams/lofar_surveys/deep/agnfitter/lockman_final_cross_match_catalogue-v0.7.fits'
        ftoutlier = '/beegfs/lofar/deepfields/science_ready_catalogs/LH_opt_spitzer_merged_vac_opt3as_irac4as_all_hpx_public.fits'
        
        ftout = '/beegfs/lofar/wwilliams/lofar_surveys/deep/agnfitter/LH_final_cross_match_catalogue-v0.7_outlier.fits'
        
    elif field == 'Bootes':
        ft = '/beegfs/lofar/wwilliams/lofar_surveys/deep/agnfitter/bootes_final_cross_match_catalogue-v0.7.fits'
        ftoutlier = '/beegfs/lofar/deepfields/science_ready_catalogs/Bootes_opt_spitzer_merged_vac_opt3as_irac4as_all_hpx_public.fits'
        
        ftout = '/beegfs/lofar/wwilliams/lofar_surveys/deep/agnfitter/Bootes_final_cross_match_catalogue-v0.7_outlier.fits'
        
    else:
        print('not handled:',field)
        sys.exit(1)
        
        
        
    ftmatch = 'temp_match.fits'



    os.system('stilts tmatch2 in1={in1} in2={in2} out={out} values1=ID values2=id join=all1 find=best suffix2=outlier matcher=exact'.format(in1=ft, in2=ftoutlier, out=ftmatch))


    t = Table.read(ft)
    toutlier = Table.read(ftmatch)


    for col in t.colnames:
        if 'flux_corr' in col:
            col1 = col.replace('_corr','')
            if col1 not in toutlier.colnames:  # handle the new FUV/NUV that aren't in Ken's files
                print(col1)
                continue
            Table([t[col], toutlier[col1]])
            print( col1, col, np.sum(toutlier[col1] == -90),np.sum(np.isnan(t[col])),)
            t[col][toutlier[col1] == -90] = np.nan
            print(np.sum(np.isnan(t[col])))

    t.write(ftout,overwrite=True)


