
## Fetch data from FIMM cluster
## scRNAseq-data
me=$(whoami)
scrna_files=outs/filtered_feature_bc_matrix
molecule_file=outs/molecule_info.h5
local_folder=/Users/$me/Dropbox/NordCML/data/scRNAseq


scrna_path=/fas/NGS/pipes/cellranger/fimm_sca_mustjoki/TIM3/Batch3-7_100619-120619/count_190701_A00464_0079_BHKLL7DMXX
names=$(echo 2260_C22D1 2414_BM_C3D1 2474_BM_scr 2298_BM_C3D1 2414_BM_scr 2515_BM_C3D1 2298_BM_scr 2474_BM_C3D1 2515_BM_scr)

scp -r jhuuhtan@atlas.genome.helsinki.fi:/fas/NGS/pipes/cellranger/fimm_sca_mustjoki/TIM3/Batch3-7_100619-120619/count_190701_A00464_0079_BHKLL7DMXX/2474_BM_C3D1/outs/molecule_info.h5 $local_folder/


for name in $names; do
  mkdir $local_folder/$name/;
  # scp -r jhuuhtan@atlas.genome.helsinki.fi:$scrna_path/$name/$scrna_files $local_folder/$name/
  scp -r jhuuhtan@atlas.genome.helsinki.fi:$scrna_path/$name/$molecule_file $local_folder/$name/
done















## ============= TCRab data
file1=all_contig_annotations.csv
file2=clonotypes.csv
file3=consensus_annotations.csv
file4=filtered_contig_annotations.csv
file5=metrics_summary.csv
tcr_files=$(echo $file1 $file2 $file3 $file4 $file5)

local_folder=/Users/$me/Dropbox/NordCML/data/scRNAseq+TCRseq/
mkdir $local_folder


## Batch1
remote_folder=/fas/NGS/pipes/cellranger/fimm_sca_mustjoki/NordCML007/Batch1-5_151019-121119/vdj_191204_A00464_0135_BHHVTCDRXX
names=$(echo 706_dg_TCR 706_3mo_TCR 706_12mo_TCR
             716_dg_TCR 716_3mo_TCR 716_12mo_TCR
             720_dg_TCR 720_3mo_TCR 720_12mo_TCR
             730_dg_TCR 730_3mo_TCR 730_12mo_TCR)

names=$(echo  730_3mo_TCR 730_12mo_TCR)

for name in $names; do
  echo $name
  mkdir $local_folder/$name/;
  for file in $tcr_files; do
    echo $file
    scp jhuuhtan@atlas.genome.helsinki.fi:$remote_folder/$name/outs/$file $local_folder/$name;
  done;
done
