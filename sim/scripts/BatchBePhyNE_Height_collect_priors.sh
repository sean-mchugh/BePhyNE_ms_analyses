LSF_DOCKER_VOLUMES="/storage1/fs1/michael.landis/Active:/storage1/fs1/michael.landis/Active"
JOBDIR="/storage1/fs1/michael.landis/Active/Sean/BePhyNE_sims"



HEIGHT_MEANS="0.5 0.7 0.95"
HEIGHT_SDS="0.15 0.5 1.0"

for hmean in $HEIGHT_MEANS
do
  for hsd in $HEIGHT_SDS
  do

    NAME="Collect_mean_${hmean}_sd_${hsd}"

    bsub -G compute-michael.landis \
      -cwd "$JOBDIR" \
      -o "$JOBDIR/stdout/$NAME.out" \
      -J "$NAME" \
      -q general \
      -g /m.seanwmchugh/BePhyNE_collect \
      -n 1 -M 3GB -R "rusage[mem=3GB] span[hosts=1]" \
      -a 'docker(sswiston/phylo_docker:full_amd64)' \
      /bin/bash "$JOBDIR/BePhyNE_Height_collect_priors.sh" "$hmean" "$hsd"

  done
done

