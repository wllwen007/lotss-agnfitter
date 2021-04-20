qsub -t 100-299%10 -l walltime=24:00:00 -v FIELD=Bootes_opt,NPER=10 -N agnf-Bootes_opt run_agnfitter.qsub
qsub -t 100-299%10 -l walltime=24:00:00 -v FIELD=LH_opt,NPER=10 -N agnf-LH_opt run_agnfitter.qsub
qsub -t 100-299%10 -l walltime=24:00:00 -v FIELD=EN1_opt,NPER=10 -N agnf-EN1_opt run_agnfitter.qsub

qsub -t 100-299%10 -l walltime=24:00:00 -v FIELD=Bootes,NPER=10 -N agnf-Bootes run_agnfitter.qsub
qsub -t 100-299%10 -l walltime=24:00:00 -v FIELD=LH,NPER=10 -N agnf-LH run_agnfitter.qsub
qsub -t 100-299%10 -l walltime=24:00:00 -v FIELD=EN1,NPER=10 -N agnf-EN1 run_agnfitter.qsub
