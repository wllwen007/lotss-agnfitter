from astropy.table import Table, Column, join
import numpy as np

import sys
field = sys.argv[1]

if field=='all':
    fields = ['Bootes','EN1','LH']
else:
    fields = [field]
    
    
    
for field in fields:    

    if field == 'EN1':

        tnew = Table.read('/beegfs/lofar/deepfields/ELAIS_N1_FIR_prelim/XID+_lofar_ELAIS-N1_v0.5_20200113.fits')
        #t = Table.read('/beegfs/lofar/deepfields/data_release/en1/final_cross_match_catalogue-v0.5.fits')
        #t = Table.read('/beegfs/lofar/deepfields/lgz/en1/final-v0.6.fits')
        t = Table.read('/beegfs/lofar/deepfields/lgz/en1/final-v0.7.fits')
        
        #ftout ='/beegfs/lofar/wwilliams/lofar_surveys/deep/agnfitter/EN1_final_cross_match_catalogue-v0.5_mergenew.fits'
        #ftout ='/beegfs/lofar/wwilliams/lofar_surveys/deep/agnfitter/EN1_final_cross_match_catalogue-v0.6_mergenew.fits'
        ftout ='/beegfs/lofar/wwilliams/lofar_surveys/deep/agnfitter/EN1_final_cross_match_catalogue-v0.7_mergenew.fits'


        tmoc = Table.read('/beegfs/lofar/deepfields/ELAIS_N1_optical/catalogues/correct_merging/add_uncat/EN1_MASTER_opt_spitzer_merged_cedit_apcorr_adduncat_extra_run2_ecorr_dr-1_lite_fservs.fits')

    elif field == 'LH':
        
        tnew = Table.read('/beegfs/lofar/deepfields/Lockman_FIR/XID+_lofar_Lockman_v0.5_20200303.fits')
        #t = Table.read('/beegfs/lofar/deepfields/data_release/lockman/final_cross_match_catalogue-v0.5.fits')
        #t = Table.read('/beegfs/lofar/deepfields/lgz/lockman/final-v0.6.fits')
        t = Table.read('/beegfs/lofar/deepfields/lgz/lockman/final-v0.7.fits')
        
        #ftout ='/beegfs/lofar/wwilliams/lofar_surveys/deep/agnfitter/LH_final_cross_match_catalogue-v0.5_mergenew.fits'
        #ftout ='/beegfs/lofar/wwilliams/lofar_surveys/deep/agnfitter/LH_final_cross_match_catalogue-v0.6_mergenew.fits'
        ftout ='/beegfs/lofar/wwilliams/lofar_surveys/deep/agnfitter/LH_final_cross_match_catalogue-v0.7_mergenew.fits'
        
    elif field == 'Bootes':
        
        tnew = Table.read('/beegfs/lofar/deepfields/Bootes_FIR/XID+_lofar_Bootes_v0.5_20200209.fits')
        #t = Table.read('/beegfs/lofar/deepfields/data_release/bootes/final_cross_match_catalogue-v0.5.fits')
        #t = Table.read('/beegfs/lofar/deepfields/lgz/bootes/final-v0.6.fits')
        t = Table.read('/beegfs/lofar/deepfields/lgz/bootes/final-v0.7.fits')
        
        #ftout ='/beegfs/lofar/wwilliams/lofar_surveys/deep/agnfitter/Bootes_final_cross_match_catalogue-v0.5_mergenew.fits'
        #ftout ='/beegfs/lofar/wwilliams/lofar_surveys/deep/agnfitter/Bootes_final_cross_match_catalogue-v0.6_mergenew.fits'
        ftout ='/beegfs/lofar/wwilliams/lofar_surveys/deep/agnfitter/Bootes_final_cross_match_catalogue-v0.7_mergenew.fits'
        
    else:
        
        print('not handled:',field)
        sys.exit(1)
        
        
        
    #toutlier = Table.read('/beegfs/lofar/deepfields/science_ready_catalogs/EN1_opt_spitzer_merged_vac_opt3as_irac4as_all_hpx_public.fits')

    colsnew = tnew.colnames
    cols = t.colnames

    copycols = []
    for cc in colsnew:
        if cc in cols:
            copycols.append(cc)

    # for the sources in tnew - replace the columns
    minds = []
    for ti in tnew:
        mind = np.where(t['Source_Name'] == ti['Source_Name'])[0][0]
        minds.append(mind)
        print ti['Source_Name'],mind
    for cc in copycols:
        t[cc][minds] = tnew[cc]


    ## add a flag for these sources...
    newflag = np.zeros(len(t), dtype=bool)
    newflag[minds] = 1
    t.add_column(Column(name='newXIDp', data=newflag))


    # add SERVS moc info for field EN1 only
    if field == 'EN1':
        # for the sources in tmoc - add the column
        minds = []
        mindst = []
        for ti in t:
            if (ti['ID'] in tmoc['ID']):
                mind = np.where(ti['ID'] == tmoc['ID'])[0][0]
                minds.append(mind)
                mindst.append(True)
                #print ti['Source_Name'],mind
            else:
                mindst.append(False)
            
        t.add_column(Column(name='FLAG_SERVS', data=np.zeros(len(t),dtype=bool)))
        t['FLAG_SERVS'] = tmoc['FLAG_SERVS'][minds]


    # save
    t.write(ftout,overwrite=True)


