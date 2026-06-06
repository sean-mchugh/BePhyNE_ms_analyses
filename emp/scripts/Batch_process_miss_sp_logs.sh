LSF_DOCKER_VOLUMES="/storage1/fs1/michael.landis/Active:/storage1/fs1/michael.landis/Active"
JOBDIR="/storage1/fs1/michael.landis/Active/Sean/BePhyNE_pletho"



    NAME="Collect_miss_logs"

    bsub -G compute-michael.landis \
      -cwd "$JOBDIR" \
      -o "$JOBDIR/stdout/$NAME.out" \
      -J "$NAME" \
      -q general \
      -g /m.seanwmchugh/BePhyNE_process \
      -n 1 -M 3GB -R "rusage[mem=3GB] span[hosts=1]" \
      -a 'docker(sswiston/phylo_docker:full_amd64)' \
      /bin/bash "$JOBDIR/process_missing_sp_logs.sh" 


