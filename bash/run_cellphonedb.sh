
## Activate virtual env
me=$(whoami)
# source /Users/$me/Dropbox/lag3/applications/cellphone_db_venv/bin/activate
source /Users/janihuuh/Dropbox/AML_TIM3/applications/cpdb-venv/bin/activate
cd /Users/$me/Dropbox/NordCML/results/cellphonedb/

## Use cellphonedb

######## tot_0m
cellphonedb method statistical_analysis --iterations=1000 --threads=20 \
  --counts-data hgnc_symbol \
  --project-name tot_0m \
  input_files/tot_0m_meta.txt input_files/tot_0m_counts.txt

######## tot_1m
cellphonedb method statistical_analysis --iterations=1000 --threads=20 \
  --counts-data hgnc_symbol \
  --project-name tot_1m \
  input_files/tot_1m_meta.txt input_files/tot_1m_counts.txt

######## tot_3m
cellphonedb method statistical_analysis --iterations=1000 --threads=20 \
  --counts-data hgnc_symbol \
  --project-name tot_3m \
  input_files/tot_3m_meta.txt input_files/tot_3m_counts.txt




######## r_0m
cellphonedb method statistical_analysis --iterations=1000 --threads=20 \
  --counts-data hgnc_symbol \
  --project-name r_0m \
  input_files/r_0m_meta.txt input_files/r_0m_counts.txt

######## r_1m
cellphonedb method statistical_analysis --iterations=1000 --threads=20 \
  --counts-data hgnc_symbol \
  --project-name r_1m \
  input_files/r_1m_meta.txt input_files/r_1m_counts.txt

######## r_3m
cellphonedb method statistical_analysis --iterations=1000 --threads=20 \
    --counts-data hgnc_symbol \
    --project-name r_3m \
    input_files/r_3m_meta.txt input_files/r_3m_counts.txt






######## n_0m
cellphonedb method statistical_analysis --iterations=1000 --threads=20 \
  --counts-data hgnc_symbol \
  --project-name n_0m \
  input_files/n_0m_meta.txt input_files/n_0m_counts.txt

######## n_1m
cellphonedb method statistical_analysis --iterations=1000 --threads=20 \
  --counts-data hgnc_symbol \
  --project-name n_1m \
  input_files/n_1m_meta.txt input_files/n_1m_counts.txt

######## n_3m
cellphonedb method statistical_analysis --iterations=1000 --threads=20 \
    --counts-data hgnc_symbol \
    --project-name n_3m \
    input_files/n_3m_meta.txt input_files/n_3m_counts.txt
