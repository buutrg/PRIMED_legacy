
########################################
# Compute sscore files
# scorefile: the score file from TERRA
# bed: the bed file prefix
########################################


ncols=$(zcat ${scorefile} | head -n1 | awk "{print NF}")

zcat ${scorefile} | awk '{$1="chr"$1; print $0}' OFS="\t" | \
awk '{
    if (FNR==1) { print $0; next}
    split($1, a, ":")
    comp['A']='T'; comp['T']='A'; comp['C']='G'; comp['G']='C'
    if (a[4] != $2 && a[4] != comp[$2]) {
        for (i=3; i<=NF; i++) {
            $i=$i*-1
        }
    }
    $2=a[4]
    print $0
}' OFS="\t" > ${scorefile}.tmp

plink2 --bfile ${bed} --score ${scorefile}.tmp no-mean-imputation header-read cols=+scoresums --score-col-nums 3-${ncols} --out ${OUT_SCORE_PATH}/chr${chr}_part${part_i}_sub${part_sub}


