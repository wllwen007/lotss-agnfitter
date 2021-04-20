qsub -t 0-317%20 -l walltime=24:00:00 -v FIELD=EN1,NPER=100,VERSION='v0.6' -N agnf-EN1 run_agnfitter.qsub
qsub -t 0-192%20 -l walltime=24:00:00 -v FIELD=Bootes,NPER=100,VERSION='v0.6' -N agnf-Bootes run_agnfitter.qsub
qsub -t 0-312%20 -l walltime=24:00:00 -v FIELD=LH,NPER=100,VERSION='v0.6' -N agnf-LH run_agnfitter.qsub

qsub -t 0-99%20 -l walltime=24:00:00 -v FIELD=EN1_opt,NPER=10,VERSION='none' -N agnf-EN1-opt run_agnfitter.qsub
qsub -t 0-99%20 -l walltime=24:00:00 -v FIELD=Bootes_opt,NPER=10,VERSION='none' -N agnf-Bootes-opt run_agnfitter.qsub
qsub -t 0-99%20 -l walltime=24:00:00 -v FIELD=LH_opt,NPER=10,VERSION='none' -N agnf-LH-opt run_agnfitter.qsub

#Bootes: radio - 19179
#Bootes: optical - 2214845
#EN1: radio - 31645
#EN1: optical - 2105993
#LH: radio - 31163
#LH: optical - 3041794

