
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


#####################################
# Ancestry adjustment
#####################################

# Convert format of scores from reference population
pgscatalog-aggregate -s ref_batch1.sscore ref_batch3.sscore ref_batch3.sscore \
    --no-split \
    -o aggregated_ref \
    --verbose

# Convert format of scores from target population
pgscatalog-aggregate -s target_batch1.sscore target_batch3.sscore target_batch3.sscore \
    --no-split \
    -o aggregated_target \
    --verbose

zcat aggregated_ref/aggregated_scores.txt.gz > aggregated_ref/aggregated_scores.txt
zcat aggregated_target/aggregated_scores.txt.gz  | tail -n+2 > aggregated_target/aggregated_scores.txt
cat aggregated_ref/aggregated_scores.txt aggregated_target/aggregated_scores.txt > all_aggregated_scores.txt

# Read more about file format at https://pygscatalog.readthedocs.io/en/latest/how-to/guides/ancestry.html
pgscatalog-ancestry-adjust --agg_scores all_aggregated_scores.txt.gz \
    --psam target.psam \
    --ref_pcs ref.pcs \
    --target_pcs target.pcs \
    --reference_related ref.king.cutoff.id \
    --dataset hgdp \
    --reference reference \
    --outdir ancestry_results
