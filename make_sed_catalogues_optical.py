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

zpapply = 'atlas'

for field in fields:
    path = '/beegfs/lofar/wwilliams/lofar_surveys/deep/science_ready_catalogs/filter_information/'
    filterpathFIR = 'FIR_filters_filters/'
    if field == 'EN1':
        finfo = path+'EN1_filters.res.info'
        filterpath = 'EN1_filters_filters/'
        zpfile = path+'../zeropoint_offsets/en1_merged_zeropoint_offsets.txt'
        outfilterpathherschel = 'filters/Herschel/'
        ftrans = path+'EN1.filter.translate'
        datfile = 'agnfitter_cols_en1_wl.csv'
        catname = '/beegfs/lofar/deepfields/science_ready_catalogs/EN1_opt_spitzer_merged_vac_opt3as_irac4as_all_hpx_public.fits'
    elif field == 'Bootes':
        finfo = path+'filter.bootes_mbrown_2014a.res.info'
        filterpath = 'filter.bootes_mbrown_2014a_filters/'
        zpfile = path+'../zeropoint_offsets/bootes_merged_zeropoint_offsets.txt'
        outfilterpathherschel = 'filters/Herschel/'
        ftrans = path+'brown.zphot.2014.translate'
        datfile = 'agnfitter_cols_bootes_wl.csv'
        catname = '/beegfs/lofar/deepfields/science_ready_catalogs/LH_opt_spitzer_merged_vac_opt3as_irac4as_all_hpx_public.fits'
    elif field == 'LH':
        finfo = path+'Lockman-SWIRE_filters.res.info'
        filterpath = 'Lockman-SWIRE_filters_filters/'
        outfilterpathherschel = 'filters/Herschel/'
        zpfile = path+'../zeropoint_offsets/lh_merged_zeropoint_offsets.txt'
        ftrans = path+'LH.filter.translate'
        datfile = 'agnfitter_cols_lockman_wl.csv'
        catname = '/beegfs/lofar/deepfields/science_ready_catalogs/Bootes_opt_spitzer_merged_vac_opt3as_irac4as_all_hpx_public.fits'
    else:
        print 'not implemented'
        sys.exit(1)
        
    field = field+'_opt'
        
    outfilterpath = 'filters/'+field+'/'


    filterdat = Table.read('filters/'+datfile)
    

    if not os.path.exists(outfilterpath):
        os.mkdir(outfilterpath)

    ### read the catalogue ###
    cat = Table.read(catname)
    cat.add_column(Column(np.arange(len(cat)),'radioID'))
   
    
    # apply zeropoint offsets - multiply fluxes
    zpoffsets = Table.read(zpfile, format='ascii')
    for zpi in zpoffsets:
        zfilt = zpi['filter'] #+'_corr' ## cat actually named filt_flux_corr - not opt cat
        if zfilt in cat.colnames:
            cat[zfilt] = cat[zfilt]*zpi[zpapply]
            
            print 'scaling ',zfilt,'by', zpi[zpapply]
        else:
            print 'no column:',zfilt

    #keep_cols = ['id','z1_median']
    #if 'ID' in cat.colnames:
        #cat.rename_column('ID','id')
    #if 'Z_BEST' in cat.colnames:
        #cat.rename_column('Z_BEST','z_best')
    #if field == 'EN1':
        #outcols = [cat['Source_Name'], cat['radioID'], cat['ID'], cat['z1_median'], cat['Z_BEST'], cat['newXIDp'], cat['FLAG_SERVS']]
    #else:
    
    #make Z_BEST column
    zbest = cat['z1_median'].copy()
    hasspec = np.isfinite(cat['Z_SPEC'])
    zbest[hasspec] = cat['Z_SPEC'][hasspec]
    cat.add_column(Column(data=zbest, name='Z_BEST'))
    
    if 'id' in cat.colnames:
        cat.rename_column('id','ID')
    outcols = [ cat['ID'], cat['z1_median'], cat['Z_BEST']]
        
    filterlist = []
    for fi in range(len(filterdat)):

        fluxcol = filterdat['influx_colname'][fi]
        fluxcolname = filterdat['outflux_colname'][fi]
        errcol = filterdat['influxerr_colname'][fi]
        errcolu = filterdat['influxerru_colname'][fi]
        errcoll = filterdat['influxerrl_colname'][fi]
        errcolname = filterdat['outfluxerr_colname'][fi]
        filt =  filterdat['id'][fi]
        
        fluxcol = fluxcol.replace('_corr','')
        errcol = errcol.replace('_corr','')
        
        filterlist.append(filterdat['id'][fi])
        fluxes = cat[fluxcol]
        fluxes[fluxes==-99] = np.nan
        fluxes[fluxes==-90] = np.nan
        if fluxes.unit=='muJy':
            fluxes.unit='microJansky'
        ### mips is unitless but is in mJy (maybe)
        if fluxes.unit is None:
            
            if ('PACS' in fluxcolname) or ('SPIRE' in fluxcolname):
                fluxes = fluxes *au.mJy
                print 'using mJy for',filt
            #elif 'MIPS' in fluxcolname:
                #fluxes = fluxes *au.microJansky
                #print 'using microJansky for',filt
            else:
                #print 'not implemented', fluxcolname
                fluxes = fluxes *au.microJansky
                print 'using muJy for',filt
                #sys.exit()
        fluxes = fluxes.to(au.Jy)
        if errcol != '':
            errfluxes = cat[errcol]
            errfluxes[errfluxes==-99] = np.nan
            errfluxes[errfluxes==-90] = np.nan
            
            if errfluxes.unit=='muJy':
                errfluxes.unit='microJansky'
            ### mips is unitless but is in mJy (maybe)
            if errfluxes.unit is None:
                
                if ('PACS' in fluxcolname) or ('SPIRE' in fluxcolname):
                    errfluxes = errfluxes *au.mJy
                    print 'using mJy for',filt
                #elif 'MIPS' in fluxcolname:
                    #errfluxes = errfluxes *au.microJansky
                    #print 'using microJansky for',filt
                else:
                    #print 'not implemented', fluxcolname
                    errfluxes = errfluxes *au.microJansky
                    print 'using muJy for',filt
                    #sys.exit()
            errfluxes = errfluxes.to(au.Jy)
        else:
            errfluxesl = cat[errcolu]
            errfluxesu = cat[errcoll]
            errfluxesl[errfluxesl==-99] = np.nan
            errfluxesu[errfluxesu==-99] = np.nan
            errfluxesl[errfluxesl==-90] = np.nan
            errfluxesu[errfluxesu==-90] = np.nan
            if errfluxesl.unit=='muJy':
                errfluxesl.unit='microJansky'
            if errfluxesl.unit is None:
                if ('PACS' in fluxcolname) or ('SPIRE' in fluxcolname):
                    errfluxesl = errfluxesl *au.mJy
                    print 'using mJy for',filt
                #elif 'MIPS' in fluxcolname:
                    #errfluxesl = errfluxesl *au.microJansky
                    #print 'using microJansky for',filt
                else:
                    #print 'not implemented', fluxcolname
                    errfluxesl = errfluxesl *au.microJansky
                    print 'using muJansky for',filt
            if errfluxesu.unit=='muJy':
                errfluxesu.unit='microJansky'
            if errfluxesu.unit is None:
                if ('PACS' in fluxcolname) or ('SPIRE' in fluxcolname):
                    errfluxesu = errfluxesu *au.mJy
                    print 'using mJy for',filt
                #elif 'MIPS' in fluxcolname:
                    #errfluxesu = errfluxesu *au.microJansky
                    #print 'using microJansky for',filt
                else:
                    #print 'not implemented', fluxcolname
                    errfluxesu = errfluxesu *au.microJansky
                    print 'using microJansky for',filt
            errfluxesl = errfluxesl.to(au.Jy)
            errfluxesu = errfluxesu.to(au.Jy)
            errfluxesl = np.abs(fluxes-errfluxesl)
            errfluxesu = np.abs(errfluxesu-fluxes)
            errfluxes = np.nanmax((errfluxesl, errfluxesu),axis=0) *au.Jy
            
        # add 10% flux in quadrature to flux errors
        errfluxes_raw = errfluxes.copy()
        errfluxes = np.sqrt(errfluxes**2. + (0.1*fluxes)**2. )
            
        # save columns
        outcols.append(Column(fluxes, fluxcolname, unit=au.Jy))
        outcols.append(Column(errfluxes, errcolname, unit=au.Jy))
        outcols.append(Column(errfluxes_raw, errcolname+'_raw', unit=au.Jy))
        
    outcat = Table(outcols)
    outcat.write(field+'_sedfit.fits', overwrite=True)
    
    
    
    ## add _wl cols for agnfitter
    for fi in range(len(filterdat)):
        #fcol = filterdat['fluxcol'][fi]
        #if fcol == '': continue

        #fluxcol = filterdat['fluxcol'][fi]
        lamcol = filterdat['outflux_colname'][fi].replace('_f','_wl')
        lam = filterdat['af_wavelength'][fi]
        outcols.append(Column(lam*np.ones(len(cat)), lamcol, unit=au.angstrom))
        
    outcat = Table(outcols)
    
        
    ## remove duplicate bands for agnfitter
    print 'unique wavelengths', len(np.unique(filterdat['af_wavelength'])) 
    print 'wavelengths', len(filterdat['af_wavelength']>0)
    if  len(np.unique(filterdat['af_wavelength']))  != len(filterdat['af_wavelength']):
        print 'duplicate wavelengths'

    waves, counts = np.unique(filterdat['af_wavelength'], return_counts=True)
    removecols = []
    for w in waves[counts>1]:
        fiw = []
        if w == 0:
            continue
        for fi in range(len(filterdat)):
            if filterdat['af_wavelength'][fi] == w:
                fiw.append(fi)
        if len(fiw) > 1:
            print 'not implemented'
            sys.exit()
            
        #should be dealt with already
            
        ##fcol = [filterdat['outflux_colname'][ff] for ff in fiw]
        #fcol = [filterdat['outflux_colname'][ff] for ff in fiw]
        #print w, fcol
        #fecol = [filterdat['outfluxerr_colname'][ff] for ff in fiw]
        
        #for ff in fcol:
            #if 'swire' in ff:
                #fswirecol = ff
            #elif 'servs' in ff:
                #fservscol = ff
            #else:
                #print 'not implemented'
                #sys.exit()
        #for ff in fecol:
            #if 'swire' in ff:
                #feswirecol = ff
            #elif 'servs' in ff:
                #feservscol = ff
            #else:
                #print 'not implemented'
        
        
        #seler1 = outcat[fecol[1]] < outcat[fecol[0]]  # er1 is smaller so should use flux 1
        #cat[seler1] = 
        #outcat[fcol[0]][seler1] = outcat[fcol[1]][seler1]
        #outcat[fecol[0]][seler1] = outcat[fecol[1]][seler1]
        #selservs = np.where(outcat['FLAG_SERVS']==1)[0]   # we should use servs
        #outcat[fswirecol][selservs] = outcat[fservscol][selservs]
        #outcat[feswirecol][selservs] = outcat[feservscol][selservs]
        
        ## set extra row to null so don't use the other flux
        #filterdat['select'][fiw[1]] = 0 
        ## and remove them from the output columns
        #removecols.append(fservscol)
        #removecols.append(feservscol)
        #removecols.append(fservscol.replace('_f','_wl'))
        
    #filterdat_ex = filterdat[filterdat['select'] == 0]
    #filterdat = filterdat[filterdat['select'] == 1]
        
        

    for t in outcat.colnames:
        if 'raw' in t:
            removecols.append(t)
    outcat.remove_columns(removecols)
    #keep only those with ID and z
    #outcat = outcat[np.isfinite(outcat['Z_BEST'])]
    #outcat = outcat[np.isfinite(outcat['ID'])]
    outcat.write(field+'_agnfitter.fits', overwrite=True)

    #Hfilternames = {1: 'PACS_100mu.txt', 2:  'PACS_160mu.txt',  3:'SPIRE_250mu.txt', 4:'SPIRE_350mu.txt', 5:'SPIRE_500mu.txt'}
    with open(field+'_agnfitter_filters.lst','w') as f:
        for fi in range(len(filterdat)):
            #fcol = filterdat['fluxcol'][fi]
            #if fcol == '': continue
            ffilter = filterdat['af_filtername'][fi]
            f.write(outfilterpath+ffilter+'.filter\n')
        #for hfi in Hfilternames.keys():
            #ffilter = Hfilternames[hfi]
            #f.write(outfilterpathherschel+ffilter+'\n')
        


    # create new settings file
    settingsin = '/beegfs/lofar/wwilliams/lofar_surveys/deep/agnfitter/SETTINGS_AGNfitter_deep_template.py'
    settingsout = settingsin.replace('template',field)
    with open(settingsin,'r') as fin:
        with open(settingsout,'w') as fout:
            for line in fin.readlines():
                if '##FIELD##' in line:
                    line = line.replace('##FIELD##',field)
                fout.write(line)


