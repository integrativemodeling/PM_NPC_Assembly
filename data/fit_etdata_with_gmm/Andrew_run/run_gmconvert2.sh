#!/bin/bash                         
#$ -S /bin/bash                     
#$ -o log 
#$ -e err 
#$ -cwd                            
#$ -pe smp 6 
#$ -r n                           
#$ -j y                           
#$ -l h_rt=48:00:00                
#$ -l h_vmem=24G

# add gmconvert to path
export PATH="$PATH:/.../gmconvert/"

# required to limit scikit learn's attempt to utilize multithreading 
export OMP_NUM_THREADS=1

MAXMEM=24000
NGAUSS=150
RVAL=0
MRCGW=10
NSLOTS=6
DATADIR="/.../data/fit_etdata_with_gmm/Andrew_run/mask_densities"
OUTDIR="/.../data/fit_etdata_with_gmm/Andrew_run/data"


KEYS=("5min" "6min" "8min" "10min" "15min" "mature" "mature_no_channel")
DFILES=("5min_masked_${RVAL}.mrc"\
        "6min_masked_${RVAL}.mrc"\
        "8min_masked_${RVAL}.mrc"\
        "10min_masked_${RVAL}.mrc"\
        "15min_masked_${RVAL}.mrc"\
        "mature_masked_${RVAL}.mrc"\
        "mature_masked_no_channel_${RVAL}.mrc")

# run gmconvert
/wynton/home/sali/aplatham/Programs/parallel-20221022/bin/bin/parallel --link --jobs $NSLOTS gmconvert -imap ${DATADIR}/{1} -ogmm $OUTDIR/{2}_${SGE_TASK_ID}.pdb -zth 0.0 -ng $NGAUSS -max_memory $MAXMEM -rseed {3} ::: ${DFILES[@]} ::: ${KEYS[@]} ::: $( for val in `seq 1 $NSLOTS`; do echo $RANDOM; done )

/wynton/home/sali/aplatham/Programs/parallel-20221022/bin/bin/parallel --jobs $NSLOTS gmconvert G2V -igmm $OUTDIR/{}_${SGE_TASK_ID}.pdb -omap $OUTDIR/{}_${SGE_TASK_ID}.mrc -gw $MRCGW -max_memory $MAXMEM ::: ${KEYS[@]}

# let parallel calls terminate
wait

date
hostname

qstat -j $JOB_ID
