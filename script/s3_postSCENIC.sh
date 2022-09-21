
echo "run postSCENIC for avg20.rep1"
f_loom_path_scenic=s2/s2_avg20_rep1.pyscenic.loom
ctx_output=s2/s2_avg20_rep1.reg.tsv
sample_name=s3/s3_avg20_rep1
threads=10
min_regulon_size=10

python s3_postSCENIC.py $f_loom_path_scenic $ctx_output $sample_name $threads $min_regulon_size

echo "run postSCENIC for avg20.rep2"
f_loom_path_scenic=s2/s2_avg20_rep2.pyscenic.loom
ctx_output=s2/s2_avg20_rep2.reg.tsv
sample_name=s3/s3_avg20_rep2

python s3_postSCENIC.py $f_loom_path_scenic $ctx_output $sample_name $threads $min_regulon_size

echo "run postSCENIC for avg20.rep3"
f_loom_path_scenic=s2/s2_avg20_rep3.pyscenic.loom
ctx_output=s2/s2_avg20_rep3.reg.tsv
sample_name=s3/s3_avg20_rep3

python s3_postSCENIC.py $f_loom_path_scenic $ctx_output $sample_name $threads $min_regulon_size

echo "run postSCENIC for avg2.rep1"
f_loom_path_scenic=s2/s2_avg2_rep1.pyscenic.loom
ctx_output=s2/s2_avg2_rep1.reg.tsv
sample_name=s3/s3_avg2_rep1

python s3_postSCENIC.py $f_loom_path_scenic $ctx_output $sample_name $threads $min_regulon_size

echo "run postSCENIC for avg2.rep2"
f_loom_path_scenic=s2/s2_avg2_rep2.pyscenic.loom
ctx_output=s2/s2_avg2_rep2.reg.tsv
sample_name=s3/s3_avg2_rep2

python s3_postSCENIC.py $f_loom_path_scenic $ctx_output $sample_name $threads $min_regulon_size

echo "run postSCENIC for avg2.rep3"
f_loom_path_scenic=s2/s2_avg2_rep3.pyscenic.loom
ctx_output=s2/s2_avg2_rep3.reg.tsv
sample_name=s3/s3_avg2_rep3

python s3_postSCENIC.py $f_loom_path_scenic $ctx_output $sample_name $threads $min_regulon_size
