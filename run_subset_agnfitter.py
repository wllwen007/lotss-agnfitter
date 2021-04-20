import sys
import os
import subprocess as sub
import datetime
import select
import numpy as np
from astropy.table import Table
from astropy import units as u
import argparse
import linecache
import os

def parse_output(af_outfile):
    fits_out = af_outfile.replace('.txt', '.ascii')
    fits_outt = af_outfile.replace('.txt', '.ascii.tmp')
    
    idn = linecache.getline(af_outfile, 1).strip().split()[-1]
    linehead = linecache.getline(af_outfile, 4).strip() + ' '
    linehead = linehead.replace('log Mstar','log_Mstar')
    linehead = linehead.replace('-ln_like','neg_ln_like')
    header = '# radioID '+linehead.replace(' ','_med ')+linehead.replace(' ','_p16 ')+linehead.replace(' ','_p84 ')+linehead.replace(' ','_ml ')

    lines0 = header+'\n'

    errlowline = linecache.getline(af_outfile, 6).strip()
    medianline = linecache.getline(af_outfile, 7).strip()
    errhighline = linecache.getline(af_outfile, 8).strip()
    mlline = linecache.getline(af_outfile, 10).strip()

    newline = idn+' '+medianline+' '+errlowline+' '+errhighline+' '+mlline

    lines0 += newline+'\n'
    cat = Table.read(lines0, format='ascii')
    cat.write(fits_outt, format='ascii', overwrite=True)
    
    # strip to single line, and remove ascii file
    os.system("tail -n1 {inf} > {out}".format(inf=fits_outt, out=fits_out))
    os.system('rm -rf {inf}'.format(inf=fits_outt))

    
    return 

def run_log(cmd,logfile,quiet=False):
    logfile = open(logfile, 'w')
    logfile.write('Running process with command: '+cmd+'\n')
    proc=sub.Popen(cmd, shell=True, stdout=sub.PIPE, stderr=sub.STDOUT,universal_newlines=True)
    while True:
        try:
            select.select([proc.stdout],[],[proc.stdout])
        except select.error:
            pass
        line=proc.stdout.readline()
        if line=='':
            break
        if not quiet:
            sys.stdout.write(line)
        ts='{:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now())
        logfile.write(ts+': '+line)
        logfile.flush()
    retval=proc.wait()
    logfile.write('Process terminated with return value %i\n' % retval)
    return retval



parser = argparse.ArgumentParser()

parser.add_argument("field", type=str, help="field to run")
parser.add_argument("--subset_ind", type=int, default=0, help="subset index to start at")
parser.add_argument("--npersubset", type=int, default=1, help="number of sources to process in this job")
parser.add_argument("-c","--clobber", action="store_true", default=False, help="overwrite any existing outputs")

args = parser.parse_args()

field = args.field
subset_ind = args.subset_ind
npersubset = args.npersubset
clobber = args.clobber


exec(compile(open('SETTINGS_AGNfitter_deep_{field}.py'.format(field=field), "rb").read(), 'SETTINGS_AGNfitter_deep_{field}.py'.format(field=field), 'exec'))
cat = CATALOG_settings()
sourcecat = Table.read(cat['filename'])

print(len(sourcecat), 'sources in source catalogue',cat['filename'])
print("will run jobs", list(range(subset_ind*npersubset, (subset_ind+1)*npersubset)))

for job_ind in  range(subset_ind*npersubset, (subset_ind+1)*npersubset):
    
    print("running job {job_ind}".format(job_ind=job_ind))


    jobdir = '/beegfs/lofar/wwilliams/lofar_surveys/deep/agnfitter/{field}/agnfitter_out/{job_ind}'.format(job_ind=job_ind, field=field)
    if os.path.isdir(jobdir) and clobber:
        os.system('rm -rf '+jobdir)
    if not os.path.isdir(jobdir):
        os.system('mkdir -p '+jobdir)    

    cmd = 'RUN_AGNfitter_multi.py  --independent --sourcenumber {job_ind} SETTINGS_AGNfitter_deep_{field}.py'.format(job_ind=job_ind, field=field)

    logfile = os.path.join(jobdir, 'agnfitter_{job_ind}.log'.format(job_ind=job_ind))


    # we need a redshift
    if np.isfinite(sourcecat['Z_BEST'][job_ind]):
        
        # don't run again if this already exists
        # it will have been removed if clobber is on
        if not os.path.isfile('{jobdir}/SED_manyrealizations_{job_ind}.png'.format(job_ind=job_ind, jobdir=jobdir)):
            print(cmd) 
            run_log(cmd, logfile)


        if os.path.isfile('{jobdir}/SED_manyrealizations_{job_ind}.png'.format(job_ind=job_ind, jobdir=jobdir)):
            
            parse_output('{jobdir}/parameter_outvalues_{job_ind}.txt'.format(job_ind=job_ind, jobdir=jobdir))
            
            with open('{jobdir}/done_{job_ind}'.format(job_ind=job_ind,jobdir=jobdir),'w') as flog:
                flog.write('SUCCESS {:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()))
                
            print ('source done')
            # these take up too much disk space
            os.system('rm {jobdir}/samples*'.format(jobdir=jobdir))
            os.system('rm {jobdir}/MODELSDICT_{job_ind}'.format(jobdir=jobdir, job_ind=job_ind))
        else:
            with open('{jobdir}/fail_{job_ind}'.format(job_ind=job_ind,jobdir=jobdir),'w') as flog:
                flog.write('FAIL {:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()))
            print ('source failed {job_ind}'.format(job_ind=job_ind))
                
    else:
        # not possible to run
        with open('{jobdir}/na_{job_ind}'.format(job_ind=job_ind, jobdir=jobdir),'w') as flog:
            flog.write('NOT RUN {:%Y-%m-%d %H:%M:%S}'.format(datetime.datetime.now()))
        print ('source not applicable {job_ind}'.format(job_ind=job_ind))
        

    print("done job {job_ind}".format(job_ind=job_ind))
    print()
