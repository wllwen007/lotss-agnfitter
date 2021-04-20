from astropy.table import Table, Column, join
import numpy as np

tnew = Table.read('/beegfs/lofar/deepfields/ELAIS_N1_FIR_prelim/XID+_lofar_ELAIS-N1_v0.5_20200113.fits')
t = Table.read('/beegfs/lofar/deepfields/data_release/en1/final_cross_match_catalogue-v0.5.fits')


tmoc = Table.read('/beegfs/lofar/deepfields/ELAIS_N1_optical/catalogues/correct_merging/add_uncat/EN1_MASTER_opt_spitzer_merged_cedit_apcorr_adduncat_extra_run2_ecorr_dr-1_lite_fservs.fits')

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

#t['ID'][t['ID'].mask] =-1
#tmoc['ID'][tmoc['ID'].mask] =-2 ## there shouldn't be any in here...
#tt = join(t, tmoc, keys=['ID','ALPHA_J2000','DELTA_J2000'])
#tt['ID'][tt['ID']==-1] = np.nan


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
t.write('/beegfs/lofar/wwilliams/lofar_surveys/deep/agnfitter/EN1_final_cross_match_catalogue-v0.5_mergenew.fits',overwrite=True)


