#!/bin/bash
round=7
operation=heating
operation2=heating

mkdir dumpoutputs/round${round}
mkdir dumpseriesoutputs/round${round}

sbatch GAP_MD_farm_submission.sh ${operation} amorph round${round}
#sbatch GAP_MD_farm_submission.sh ${operation} divacancy round${round}
#sbatch GAP_MD_farm_submission.sh ${operation} interstitial round${round}
#sbatch GAP_MD_farm_submission.sh ${operation2} liquid round${round}
#sbatch GAP_MD_farm_submission.sh ${operation} vacancy round${round}

echo "All batches submitted."
exit
