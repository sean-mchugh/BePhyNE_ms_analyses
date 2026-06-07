LSF_DOCKER_VOLUMES="/storage1/fs1/michael.landis/Active:/storage1/fs1/michael.landis/Active"
JOBDIR="/storage1/fs1/michael.landis/Active/Sean/BePhyNE_sims"

RUN_LIST=$(seq 1 50)

HEIGHT_MEANS="0.5 0.7 0.95"
HEIGHT_SDS="0.15 0.5 1.0"

for r in $RUN_LIST
do
  for hmean in $HEIGHT_MEANS
  do
    for hsd in $HEIGHT_SDS
    do

      NAME="Height_runs_${r}_mean_${hmean}_sd_${hsd}"

      bsub -G compute-michael.landis \
        -cwd "$JOBDIR" \
        -o "$JOBDIR/stdout/$NAME.out" \
        -J "$NAME" \
        -q general \
        -g /m.seanwmchugh/BePhyNE_Kappa \
        -n 1 -M 3GB -R "rusage[mem=3GB] span[hosts=1]" \
        -a 'docker(sswiston/phylo_docker:full_amd64)' \
        /bin/bash "$JOBDIR/BePhyNE_Height_priors.sh" "$r" "$hmean" "$hsd"
    

    done
  done
done
