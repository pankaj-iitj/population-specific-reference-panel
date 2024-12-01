# install conda
wget https://repo.anaconda.com/archive/Anaconda3-2024.06-1-Linux-x86_64.sh
bash Anaconda3-2024.06-1-Linux-x86_64.sh

conda init
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda

conda install bcftools samtools bwe gatk4 fastqc multiqc fastp

# LDtect
virtualenv ldetect
source ldetect/bin/activate
pip install ldetect
deactivate

# snpflip
virtualenv snpflip
source snpflip/bin/activate
pip install snpflip
deactivate
