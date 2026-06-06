# batch shell script
LSF_DOCKER_VOLUMES="/storage1/fs1/michael.landis/Active:/storage1/fs1/michael.landis/Active"
JOBDIR="/storage1/fs1/michael.landis/Active/Sean/BePhyNE_pletho"

RUN_LIST=$(seq 1 20)

for run_id in ${RUN_LIST[@]}
do

  NAME="Pletho_full_${run_id}"

  bsub -G compute-michael.landis \
    -cwd "$JOBDIR" \
    -o "$JOBDIR/stdout/$NAME.out" \
    -J "$NAME" \
    -q general \
    -g /m.seanwmchugh/BePhyNE_Kappa \
    -n 1 -M 3GB -R "rusage[mem=3GB] span[hosts=1]" \
    -a 'docker(sswiston/phylo_docker:full_amd64)' \
    /bin/bash "$JOBDIR/BePhyNE_pletho_full.sh" "$run_id"

done
