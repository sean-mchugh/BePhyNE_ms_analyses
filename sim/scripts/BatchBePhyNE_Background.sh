LSF_DOCKER_VOLUMES="/storage1/fs1/michael.landis/Active:/storage1/fs1/michael.landis/Active"
JOBDIR="/storage1/fs1/michael.landis/Active/Sean/BePhyNE_sims"


#RUN_LIST=$(seq 201 400)
#
#
#for r in ${RUN_LIST[@]}
#do
#
#
#
#NAME="Background_runs_$r"
#
#bsub -G compute-michael.landis \
#-cwd /storage1/fs1/michael.landis/Active/Sean/BePhyNE_sims/ \
#-o /storage1/fs1/michael.landis/Active/Sean/BePhyNE_sims/stdout/$NAME  \
#-J $NAME \
#-q general \
#-g /m.seanwmchugh/BePhyNE_Kappa \
#-n 1 -M 3GB -R "rusage [mem=3GB] span[hosts=1]" \
#-a 'docker(sswiston/rb_tp:test)' /bin/bash /storage1/fs1/michael.landis/Active/Sean/BePhyNE_sims/BePhyNE_Background.sh
#
#done


RUN_LIST=$(seq 0 200)

for r in ${RUN_LIST[@]}; do
  for bool in TRUE FALSE; do

    NAME="${bool}_${r}"

    bsub -G compute-michael.landis \
      -cwd /storage1/fs1/michael.landis/Active/Sean/BePhyNE_sims/ \
      -o /storage1/fs1/michael.landis/Active/Sean/BePhyNE_sims/stdout/$NAME \
      -J "$NAME" \
      -q general \
      -g /m.seanwmchugh/BePhyNE_Kappa \
      -n 1 -M 3GB -R "rusage [mem=3GB] span[hosts=1]" \
      -a 'docker(sswiston/rb_tp:test)' \
      /bin/bash /storage1/fs1/michael.landis/Active/Sean/BePhyNE_sims/BePhyNE_Background.sh

  done
done