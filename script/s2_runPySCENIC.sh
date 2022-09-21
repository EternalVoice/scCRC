## inputs
sn=$1
f_loom_path_crc=s1/s1_expr.mat.loom
f_loom_path_scenic=s1/s1_${sn}.loom

## outputs
grn_output=s2/s2_${sn}.adj.tsv
ctx_output=s2/s2_${sn}.reg.tsv
f_pyscenic_output=s2/s2_${sn}.pyscenic.loom

## reference
f_tfs=/home/tongqiang/reference/cisTarget/TFs/hs_hgnc_tfs.txt
f_motif_path=/home/tongqiang/reference/cisTarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl
f_db_names=`find /home/tongqiang/reference/cisTarget/mc9nr -name "hg38*.feather"`

arboreto_with_multiprocessing.py \
    $f_loom_path_scenic \
    $f_tfs \
    --method grnboost2 \
    --output $grn_output \
    --num_workers 10 \
    --seed 777

pyscenic ctx \
    $grn_output \
    $f_db_names \
    --annotations_fname $f_motif_path \
    --expression_mtx_fname $f_loom_path_scenic \
    --output $ctx_output \
    --num_workers 10
    ## Note:
    ## The reason why I didn't use `--mask_dropouts` parameters:
    ## In R version of SCENIC, mask dropouts is False (default when pySCENIC>0.9.6, current version: 0.10.1). 
    ## Here we match the behavior of R version.


# ## set TMPDIR to current path, in case of no enough disk space on /tmp/
#export TMPDIR=`pwd` 

pyscenic aucell \
    $f_loom_path_crc \
    $ctx_output \
    --output $f_pyscenic_output \
    --num_workers 10 \
    --seed 777
