cd ./pySCENIC/scripts
python3 arboreto_with_multiprocessing.py ~/data.loom ~/mm_tfs.txt --method grnboost2 --output ~/tmp/adj.tsv --num_workers 28 --seed 777
pyscenic ctx ~/tmp/adj.tsv   ~/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather  --annotations_fname  ~/motifs-v9-nr.mgi-m0.001-o0.0.tbl --expression_mtx_fname ~/data.loom --output ~/tmp/reg.csv --mask_dropouts --num_workers 28
pyscenic aucell  --num_workers 28 ~/data.loom ~/tmp/reg.csv  -o ~/tmp/AUC.csv