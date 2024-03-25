#!/bin/bash                         
#$ -S /bin/bash                     
#$ -cwd                            
#$ -pe smp 1 
#$ -l mem_free=4G
#$ -r n                           
#$ -j y                           
#$ -m e
#$ -t 1-5
#$ -l h_rt=6:00:00                


# add gmconvert to path
export PATH="$PATH:~/local/bin"

# required to limit scikit learn's attempt to utilize multithreading 
export OMP_NUM_THREADS=1

NTRIALS=1
NGAUSS=2
MAXMEM=4000
MRCGW=10
NSLOTS=1
DATADIR="../../output/templates"

DFILES=("ir_channel_template.pdb"\
        "ir_core1_template.pdb"\
        "ir_core2_template.pdb"\
        "Nup155_template.pdb"\
        "yc_template.pdb")

let kinx=$SGE_TASK_ID-1

# boilerplate ensuring that TMPDIR is created; if not make new dir under username
if [[ -z "$TMPDIR" ]]; then
    if [[ -d /scratch ]]; then TMPDIR=/scratch/$USER; else TMPDIR=/tmp/$USER; fi
    mkdir -p "$TMPDIR"
    export TMPDIR
fi


# use this bash trick to sub output file name extension 
ofile="${infile%.*}_.gmm"

# uses GNU parallel to execute multiple fitting runs
# convert atomic template subcomplex using gmconvert
parallel --link --jobs $NSLOTS gmconvert -ipdb ${DATADIR}/${DFILES[$kinx]} -ogmm $TMPDIR/${DFILES[$kinx]%.*}_gmm_{1}.pdb -ng $NGAUSS -max_memory $MAXMEM -rseed {2} ::: $( seq 1 $NTRIALS ) ::: $( for val in `seq 1 $NTRIALS`; do echo $RANDOM; done )
wait

# convert output to mrc 
parallel --jobs $NSLOTS gmconvert G2V -igmm $TMPDIR/${DFILES[$kinx]%.*}_gmm_{1}.pdb -omap $TMPDIR/${DFILES[$kinx]%.*}_{1}.mrc -gw $MRCGW -max_memory $MAXMEM ::: $( seq 1 $NTRIALS )

# convert pdb files to isd format
#scl enable rh-python36 "parallel --jobs $NSLOTS python ../tools/pdb_to_isd_gmm.py $TMPDIR/${KEYS[$kinx]}_{1}.pdb $TMPDIR/${KEYS[$kinx]}_{1}.gmm ::: $( seq 1 $NTRIALS )" 
wait

# group output into tar files
#cd "$TMPDIR" 

#tar -cf ${DFILES[$kinx]%.*}_gmm_all.tar *.pdb
#rm *.pdb

#tar -cf ${DFILES[$kinx]%.*}_mrc_all.tar ./*.mrc
#rm *.mrc

#tar -cf ${KEYS[$kinx]}_gmm_all.tar ./*.gmm 
#rm *.gmm

# back to calling dir
#cd -

# copy files from TMP to HOME output
cp $TMPDIR/* ${DATADIR}

#python fit_gmm_gmconvert.py ${KEYS[$kinx]} $i
#gmconvert -imap ${DATA_FILE} -ogmm ${OUTFILE}.pdb -zth ${THRESHOLD} -ng $i -max_memory 3000
#gmconvert -igmm ${OUTFILE}.pdb -omap ${OUTFILE}.mrc -gw 10	

date
hostname

qstat -j $JOB_ID
