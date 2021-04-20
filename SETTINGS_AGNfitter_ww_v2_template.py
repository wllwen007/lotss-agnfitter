def CATALOG_settings():

    """==================================
    ***USER INPUT NEEDED***

    Set the right values to be able to read your catalog's format.
    FITS option is not available yet.
    =================================="""


    cat = dict()


    ##GENERAL
    cat['path'] ='/car-data/wwilliams/bootes/multi_bootes/agnfitter/git/AGNfitter/'

    cat['outpath'] ='/car-data/wwilliams/bootes/multi_bootes/agnfitter/'
    
    cat['filename'] = cat['outpath']+'/__CATFILE__'


    ##OTHERS (user input ***not*** needed)
    cat['output_folder'] =  cat['outpath']+'/__OUTDIR__/'#if no special OUTPUT folder, leave default
    #cat['dict_path'] = cat['outpath']+'models/MODELSDICT_bootesz' 
    cat['dict_path'] = cat['outpath']+'/__MODEL__'


    cat['filetype'] = 'FITS' ## catalog file type: 'ASCII' or 'FITS'. 
    cat['name'] = 'id' #'ID'            ## If ASCII: Column index (int) of source IDs
    cat['redshift'] = 'z' #'z'              ## If ASCII:  Column index(int) of redshift
 
    ##FREQUENCIES/WAVELENGTHS 
    cat['freq/wl_suffix'] = '_wl'
    cat['freq/wl_format'] = 'wavelength' ## Gives your catalog *observed*
    cat['freq/wl_unit'] = u.Angstrom       ## Astropy unit of freq or wavelength

    ##FLUXES  
    cat['flux_unit'] = u.Jy             ## Astropy unit of *flux* (astropy-units)
    cat['flux_suffix'] = '_f'
    cat['fluxerr_suffix'] = '_e'
    cat['ndflag_bool'] = True          ## Does you catalog has columns with flags 1(0) for 
    cat['ndflag_suffix'] = '_nd'         ## If ASCII: List of column indexes (int)

    return cat


def FILTERS_settings():

    """==================================
    Set the photometric bands included in your catalog,
    in order to integrate the models over their response curves.
    =================================="""

    filters = dict()
    
    #filters['dict_zarray'] = zz  # The grid of redshifts needed to fit your catalog
    #filters['dict_zarray'] = np.arange(0,7.00001,0.01)  # The grid of redshifts needed to fit your catalog
    filters['Bandset'] = 'BANDSET_settings' # OPTIONS: 
                                           # 'BANDSET_default' (for testing)
                                           # 'BANDSET_settings' (choosing relevant filters below, as given by your catalog)
                                           # if your filter is not included, go to DICTIONARIES_AGNfitter to add.

    filters['SPIRE500']= True
    filters['SPIRE350']= True
    filters['SPIRE250']= True
    filters['PACS160']=False
    filters['PACS100']=False

    filters['MIPS160']=False      
    filters['MIPS70']=False    
    filters['MIPS24']=True

    filters['IRAC4']=True       
    filters['IRAC3']=True
    filters['IRAC2']=True
    filters['IRAC1']=True

    filters['WISE4']=False
    filters['WISE3']=False
    filters['WISE2']=False
    filters['WISE1']=False

    filters['Ks_2mass']=False
    filters['H_2mass']=False
    filters['J_2mass']=False

    filters['H_VISTA']=False
    filters['J_VISTA']=False
    filters['K_VISTA']=False
    filters['Y_VISTA']=False
    filters['Z_VISTA']=False

    filters['u_SDSS']=False  
    filters['g_SDSS']=False
    filters['r_SDSS']=False
    filters['i_SDSS']=False  
    filters['z_SDSS']=False

    filters['g_SUBARU']=False
    filters['r_SUBARU']=False
    filters['i_SUBARU']=False  
    filters['z_SUBARU']=False
    filters['B_SUBARU']=False
    filters['V_SUBARU']=False

    filters['u_CHFT']=False  
    filters['g_CHFT']=False
    filters['r_CHFT']=False
    filters['i_CHFT']=False  
    filters['z_CHFT']=False

    filters['GALEX_2500']=True
    filters['GALEX_1500']=False
    
    
    filters['Y_NDWFS']=True
    filters['U_NDWFS']=True
    filters['B_NDWFS']=True
    filters['R_NDWFS']=True
    filters['I_NDWFS']=True
    filters['J_NDWFS']=True
    filters['H_NDWFS']=True
    filters['K_NDWFS']=True
    filters['z_NDWFS']=True
    
    
    return filters

def MCMC_settings():

    """==================================
    Set your preferences for the MCMC sampling.
    =================================="""

    mc = dict()

    mc['Nwalkers'] = 100  ## number of walkers 
    mc['Nburnsets']= 2   ## number of burn-in sets
    mc['Nburn'] = 4000 ## length of each burn-in sets
    mc['Nmcmc'] = 10000  ## length of each burn-in sets
    mc['iprint'] = 1000 ## show progress in terminal in steps of this many samples

    return mc

def OUTPUT_settings():

    """==================================
    Set your preferences for the production of OUTPUT files. 
    =================================="""

    out = dict()

    out['plot_format'] = 'png'
    out['plot_residuals'] = False

    #CHAIN TRACES
    out['plot_tracesburn-in'] = False    
    out['plot_tracesmcmc'] = False

    #BASIC OUTPUT
    out['Nsample'] = 1000 ## out['Nsample'] * out['Nthinning'] <= out['Nmcmc']
    out['Nthinning'] = 10 ## This describes thinning of the chain to sample
    out['writepar_meanwitherrors'] = True ##Write output values for all parameters in a file.
    out['plot_posteriortriangle'] = False ##Plot triangle with all parameters' PDFs?
    
    out['save_chains'] = True

    #INTEGRATED LUMINOSITIES
    out['calc_intlum'] = True  
    out['realizations2int'] = 1000 #This process is very time consuming.
                                #Around 100-1000 is recomendend for computational reasons.
                                #If you want to plot posterior triangles of 
                                #the integrated luminosities, should be > 1000.
    out['plot_posteriortrianglewithluminosities'] = False  # requires out['calc_intlum']=True

    #INTEGRATION RANGES
    out['intlum_models'] = ['sb','bbb', 'bbbdered', 'gal', 'gal', 'tor','sb','bbb', 'bbbdered', 'gal', 'tor','sb', 'gal', 'tor','sb']  #leave 'sb' always 
                                                                        #as first element
    out['intlum_freqranges_unit'] = u.micron   #Astropy unit 
    out['intlum_freqranges'] = np.array([[8.,1000.],[0.1,1.],[0.1,1.],[0.1,1.],[1.,100.],[1.,100.],[1.,100.],[1.,30.],[1.,30.],[1.,30.],[1.,30.],[1.,30.],[1.,8.],[1.,8.],[1.,8.]])
    out['intlum_names'] = ['LIR_8_1000','Lbb_0p1_1', 'Lbbdered_0p1_1', 'Lga_0p1_1', 'Lga_1_100', 'Ltor_1_100','Lsb_1_100','Lbb_1_30', 'Lbbdered_1_30', 'Lga_1_30', 'Ltor_1_30','Lsb_1_30', 'Lga_1_8', 'Ltor_1_8','Lsb_1_8']

    #SED PLOTTING
    out['realizations2plot'] = 10

    out['plotSEDrealizations'] = True

    return out
