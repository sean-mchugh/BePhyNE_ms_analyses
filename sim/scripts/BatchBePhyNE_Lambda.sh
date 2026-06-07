LSF_DOCKER_VOLUMES="/storage1/fs1/michael.landis/Active:/storage1/fs1/michael.landis/Active"
JOBDIR="/storage1/fs1/michael.landis/Active/BePhyNE_sims"


RUN_LIST=$(seq 1 100)


for r in ${RUN_LIST[@]}
do



NAME="Lambda_runs_$r"

bsub -G compute-michael.landis \
-cwd /storage1/fs1/michael.landis/Active/BePhyNE_sims/ \
-o /storage1/fs1/michael.landis/Active/BePhyNE_sims/stdout/$NAME  \
-J $NAME \
-q general \
-g /m.seanwmchugh/BePhyNE_Lambda \
-n 1 -M 3GB -R "rusage [mem=3GB] span[hosts=1]" \
-a 'docker(sswiston/rb_tp:6)' /bin/bash /storage1/fs1/michael.landis/Active/BePhyNE_sims/BePhyNE_Lambda.sh

done