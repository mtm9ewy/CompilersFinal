#!/bin/bash

start="$1"
end="$2"

# Arbitrary choice, but it evenly divides total
step=91
# 16 * 15 * 14 * 13
total=43680

end=$((end > total ? total : end))

while [[ $start -le $end ]]
do
    last=$((start+step-1 < end ? start+step-1 : end))
    echo "Submitting a job for permutations ${start} through ${last}"
    sbatch --exclude=adriatic01,adriatic02,adriatic03,adriatic04,adriatic05,adriatic06,affogato01,affogato02,affogato03,affogato04,affogato05,affogato11,affogato12,affogato13,affogato14,affogato15,ai01,ai02,ai03,ai04,ai05,ai06,ai07,ai08,ai09,ai10,cheetah01,cheetah02,cheetah03,doppio01,doppio02,doppio03,doppio04,doppio05,epona,heartpiece,hydro,jaguar01,jaguar02,jaguar03,jaguar04,lotus,lynx01,lynx02,lynx03,lynx04,lynx05,lynx06,lynx07,lynx08,lynx09,lynx10,lynx11,lynx12,optane01,panther01,pegasusboots,puma01,ristretto01,ristretto02,ristretto03,ristretto04,sds01,sds02,slurm1,slurm2,slurm3,slurm4,slurm5,titanx01,titanx02,titanx03,titanx04,titanx05,titanx06 runtests.sh "$start" "$last" 
    ((start+=step))
done

