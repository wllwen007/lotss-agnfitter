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

for field in fields:
    path = '/beegfs/lofar/wwilliams/lofar_surveys/deep/science_ready_catalogs/filter_information/'
    filterpathFIR = 'FIR_filters_filters/'
    if field == 'EN1':
        finfo = path+'EN1_filters.res.info'
        filterpath = 'EN1_filters_filters/'
        outfilterpathherschel = 'filters/Herschel/'
        ftrans = path+'EN1.filter.translate'
        #catname = '/beegfs/lofar/wwilliams/lofar_surveys/deep/science_ready_catalogs/EN1_select.fits'
        catname = '/beegfs/lofar/deepfields/data_release/en1/final_cross_match_catalogue-v0.5.fits'
    elif field == 'Bootes':
        finfo = path+'filter.bootes_mbrown_2014a.res.info'
        filterpath = 'filter.bootes_mbrown_2014a_filters/'
        outfilterpathherschel = 'filters/Herschel/'
        ftrans = path+'brown.zphot.2014.translate'
        #catname = '/beegfs/lofar/wwilliams/lofar_surveys/deep/science_ready_catalogs/Bootes_select.fits'
        catname = '/beegfs/lofar/deepfields/data_release/bootes/final_cross_match_catalogue-v0.5.fits'
    elif field == 'LH':
        finfo = path+'Lockman-SWIRE_filters.res.info'
        filterpath = 'Lockman-SWIRE_filters_filters/'
        outfilterpathherschel = 'filters/Herschel/'
        ftrans = path+'LH.filter.translate'
        #catname = '/beegfs/lofar/wwilliams/lofar_surveys/deep/science_ready_catalogs/LH_select.fits'
        catname = '/beegfs/lofar/deepfields/data_release/lockman/final_cross_match_catalogue-v0.5.fits'
    else:
        print 'not implemented'
        sys.exit(1)
        
    outfilterpath = 'filters/'+field+'/'

    if not os.path.exists(outfilterpath):
        os.mkdir(outfilterpath)

    filterdata = Table(names=('id', 'filtername', 'filterpath', 'fluxcol', 'fluxerrcol', 'fluxerrcoll', 'fluxerrcolu', 'fluxcolname','fluxerrcolname', 'wavelength'), dtype=('i4', 'S35', 'S200', 'S25',  'S25', 'S25', 'S25', 'S25', 'S25', 'f16'))
    def add_blankrow(t, tid):
        t.add_row((tid, '', '', '', '', '', '', '', '', 0.))
        return t

    #filternames = {}
    #filterpaths = {}
    # read *res.info file
    with open(finfo, 'r') as f:
        tt = f.readlines()
        for t in tt:
            C = t.strip().split()
            if len(C) == 0 : continue
            fi = int(C[0])
            if fi not in filterdata['id']:
                filterdata = add_blankrow(filterdata, fi)
                
            filterdata['filtername'][filterdata['id']==fi] = C[1]
            filterdata['filterpath'][filterdata['id']==fi] = filterpath
            #filternames[fi] = C[1]
            #filterpaths[fi] = filterpath



    # read the translate file 
    # this matches colname in catalog to filter number in *.res.info (to filter name)
    #fluxcols = {}
    #fluxerrcols = {}
    # the output column names (for most these are the same as the inputs)
    #fluxerrcolnames = {}
    #fluxcolnames = {}
    with open(ftrans, 'r') as f:
        tt = f.readlines()
        for t in tt:
            ### TEMP - for now for TESTING - EN1 - remove the servs flux
            if 'servs' in t: continue
        
            if '#' in t: continue
        
            C = t.strip().split()
            if len(C) == 0 : continue

        
            fi = int(C[1][1:])
            
            if fi not in filterdata['id']:
                filterdata = add_blankrow(filterdata, fi)
                print 'warning: not in res.info file???', C
                
            if C[1][0] == 'F':
                #fluxcols[fi] = C[0]
                #fluxcolnames[fi] = C[0]
                filterdata['fluxcol'][filterdata['id']==fi] = C[0]+ '_corr'
                filterdata['fluxcolname'][filterdata['id']==fi] = C[0] 
            elif C[1][0] == 'E':
                #fluxerrcols[fi] = [C[0]]
                #fluxerrcolnames[fi] = C[0]
                filterdata['fluxerrcol'][filterdata['id']==fi] = C[0]+ '_corr'
                filterdata['fluxerrcolname'][filterdata['id']==fi] = C[0]
            else:
                print 'something is wrong'
                sys.exit()

    #sys.exit()
    #need to add PACS and SPIRE
    # manually add the fir columns
    fi = max(filterdata['id']) +1 
    for addf,addfname in zip(['MIPS_24' , 'PACS_100', 'PACS_160', 'SPIRE_250', 'SPIRE_350', 'SPIRE_500'],
                            ['mips24.filter','PACS_100mu.txt','PACS_160mu.txt','SPIRE_250mu.txt','SPIRE_350mu.txt','SPIRE_500mu.txt']):

        if fi not in filterdata['id']:
            filterdata = add_blankrow(filterdata, fi)
        fii = filterdata['id']==fi
        filterdata['fluxcol'][fii] = 'F_'+addf
        filterdata['fluxcolname'][fii] = addf+'_flux'
        filterdata['fluxerrcoll'][fii] = 'FErr_'+addf+'_l'
        filterdata['fluxerrcolu'][fii] = 'FErr_'+addf+'_u'
        filterdata['fluxerrcolname'][fii] = addf+'_fluxerr'
        filterdata['filtername'][fii] = addfname
        filterdata['filterpath'][fii] = filterpathFIR
        fi += 1



    # filter files should be in Angstrom
    c=    2.997e8
    Angstrom = 1e10

    #wavelengths = {}
    for fi in range(len(filterdata)):
        fcol = filterdata['fluxcol'][fi]
        if fcol == '': continue
        ffilter = filterdata['filtername'][fi]
        band_lambda, band_factor =  np.loadtxt(path+filterdata['filterpath'][fi]+ffilter+'.filter', usecols=(0,1),unpack= True)
        os.system('cp '+path+filterdata['filterpath'][fi]+ffilter+'.filter '+outfilterpath)
        print('cp '+path+filterdata['filterpath'][fi]+ffilter+'.filter '+outfilterpath)

        central_lamb = np.sum(band_lambda*band_factor)/np.sum(band_factor)
        central_nu = float(np.log10((Angstrom*c)/central_lamb))
        filterdata['wavelength'][fi] = central_lamb
        
    #for fi in fluxcols.keys():
    #for fi in range(len(filterdata)):
        #ffilter = filternames[fi]
        #print filterdata[fi], fluxcols[fi],fluxerrcols[fi],filternames[fi],wavelengths[fi]

    ### read the catalogue ###
    cat = Table.read(catname)
    # temp for TESTING
    # TODO remove
    #cat = cat[0:2000]

    #keep_cols = ['id','z1_median']
    #if 'ID' in cat.colnames:
        #cat.rename_column('ID','id')
    #if 'Z_BEST' in cat.colnames:
        #cat.rename_column('Z_BEST','z_best')
    outcols = [cat['ID'], cat['z1_median'], cat['Z_BEST']]
    filterlist = []
    for fi in range(len(filterdata)):
        fcol = filterdata['fluxcol'][fi]
        if fcol == '': continue

        fluxcol = filterdata['fluxcol'][fi]
        fluxcolname = filterdata['fluxcolname'][fi]
        errcol = filterdata['fluxerrcol'][fi]
        errcolu = filterdata['fluxerrcolu'][fi]
        errcoll = filterdata['fluxerrcoll'][fi]
        errcolname = filterdata['fluxerrcolname'][fi]
        filt =  filterdata['filtername'][fi]
        
        filterlist.append(filterdata['filtername'][fi])
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
        outcols.append(Column(fluxes, fluxcolname, unit=au.Jy))
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
            errfluxesl = fluxes-errfluxesl
            errfluxesu = errfluxesu-fluxes
            errfluxes = np.nanmax((errfluxesl, errfluxesu),axis=0)
            
        outcols.append(Column(errfluxes, errcolname, unit=au.Jy))
        
        
    print 'unique wavelengths', len(np.unique(filterdata['wavelength'])) 
    print 'wavelengths', len(filterdata['wavelength']>0)
    if  len(np.unique(filterdata['wavelength']))  != len(filterdata['wavelength']):
        print 'duplicate wavelengths'

    waves, counts = np.unique(filterdata['wavelength'], return_counts=True)
    removecols = []
    for w in waves[counts>1]:
        fiw = []
        if w == 0:
            continue
        for fi in range(len(filterdata)):
            if filterdata['wavelength'][fi] == w:
                fiw.append(fi)
        if len(fiw) > 2:
            print 'not implemented'
            sys.exit()
            
        fcol = [filterdata['fluxcolname'][ff] for ff in fiw]
        print w, fcol
        
        fecol = [filterdata['fluxerrcolname'][ff] for ff in fiw]
        seler1 = cat[fecol[1]] < cat[fecol[0]]  # er1 is smaller so should use flux 1
        #cat[seler1] = 
        cat[fcol[0]][seler1] = cat[fcol[1]][seler1]
        cat[fecol[0]][seler1] = cat[fecol[1]][seler1]
        
        # set extra row to null so don't use the other flux
        filterdata['fluxcol'][fiw[1]] = ''  
        # and remove them from the output columns
        removecols.append(fcol[1])
        removecols.append(fecol[1])
        
    filterdata_ex = filterdata[filterdata['fluxcol'] == '']
    filterdata = filterdata[filterdata['fluxcol'] != '']
        
    for fi in range(len(filterdata)):
        fcol = filterdata['fluxcol'][fi]
        if fcol == '': continue

        fluxcol = filterdata['fluxcol'][fi]
        lamcol = filterdata['fluxcolname'][fi].replace('_flux','_wl')
        lam = filterdata['wavelength'][fi]
        outcols.append(Column(lam*np.ones(len(cat)), lamcol, unit=au.angstrom))
        

    outcat = Table(outcols)
    outcat.remove_columns(removecols)
    #keep only those with ID and z
    outcat = outcat[np.isfinite(outcat['Z_BEST'])]
    outcat = outcat[np.isfinite(outcat['ID'])]
    outcat.write(field+'_agnfitter.fits', overwrite=True)

    #Hfilternames = {1: 'PACS_100mu.txt', 2:  'PACS_160mu.txt',  3:'SPIRE_250mu.txt', 4:'SPIRE_350mu.txt', 5:'SPIRE_500mu.txt'}
    with open(field+'_agnfitter_filters.lst','w') as f:
        for fi in range(len(filterdata)):
            fcol = filterdata['fluxcol'][fi]
            if fcol == '': continue
            ffilter = filterdata['filtername'][fi]
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


