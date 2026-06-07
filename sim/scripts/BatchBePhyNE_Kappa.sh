LSF_DOCKER_VOLUMES="/storage1/fs1/michael.landis/Active:/storage1/fs1/michael.landis/Active"
JOBDIR="/storage1/fs1/michael.landis/Active/Sean/BePhyNE_sims"


RUN_LIST=$(seq 1 200)


for r in ${RUN_LIST[@]}
do



NAME="New_New_Kappa_runs_$r"

bsub -G compute-michael.landis \
-cwd /storage1/fs1/michael.landis/Active/Sean/BePhyNE_sims/ \
-o /storage1/fs1/michael.landis/Active/Sean/BePhyNE_sims/stdout/$NAME  \
-J $NAME \
-q general \
-g /m.seanwmchugh/BePhyNE_Kappa \
-n 1 -M 3GB -R "rusage [mem=3GB] span[hosts=1]" \
-a 'docker(sswiston/rb_tp:test)' /bin/bash /storage1/fs1/michael.landis/Active/Sean/BePhyNE_sims/BePhyNE_Kappa.sh

done
