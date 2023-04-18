# basecall_nanopore

## Description
This pipeline automates a series of task to convert Nanopore's fast5 file to fastq (`Guppy GPU`), perform basic QC (`pycoQC`), adapter trimming (`Porechop`) and low quality read filtering (drop bottom 5%; `Filtlong`).

## Requirements
* Linux computer equipped with a capable NVIDA graphics card and a working NVIDIA driver.
* GPU version of Guppy v6.3.8. It could work with other versions, but has not been tested.
* Requirement are listed in the `requirements.txt` file.

## Installation
Guppy v6.3.8:
```commandline
# Go to my home folder
cd ~

# Create folder for program
mkdir prog

# Change to "prog" directory
cd prog  

# Dowload Guppy GPU version 6.3.8
wget https://cdn.oxfordnanoportal.com/software/analysis/ont-guppy_6.3.8_linux64.tar.gz

# Uncompress
tar zxvf ont-guppy_6.3.8_linux64.tar.gz

# Find "guppy_basecaller" absolute path
cd ont-guppy/bin
pwd

# Add it to your PATH
# Say the "pwd" command results was /home/bioinfo/prog/ont-guppy/bin
echo "export PATH=\$PATH:/home/bioinfo/prog/ont-guppy/bin" | tee -a ~/.bashrc

# Apply changes so your bash terminal knows where to find guppy
source ~/.bashrc
```
Next we need to create and activate a virtual environment to install all the tools we need for the pipeline. We're going to use `miniconda` for this purpose (installation instruction can be found [here](https://conda.io/projects/conda/en/stable/user-guide/install/linux.html) if not already installed).

The dependencies:
```commandline
# Create environment
conda create -n nanopore -y -c bioconda porechop filtlong psutil pandas pycoqc pysam
# Activate environment
conda activate nanopore
```
The pipeline itself:
```commandline
# Clone GitHub repository
git clone https://github.com/duceppemo/basecall_nanopore

# Test (don't forgfet to avitvate the "nanopore" environment first if not already done)
python basecall_nanopore.py -h  # Should not display any error and show help message
```

## Usage
There are two ways to use Guppy:
1- Using a configuration files. Preferred method.
2- Using a combination of sequencing kit, flowcell and sequencer.

The allowed values for the configuration file, the library preparation kit, barcode kit are located in `kits.py` file.

## About barcodes
1- If samples were barcoded, providing the specific barcode kit used will speed up the basecalling/demultiplexing.
2- If the run contained barcodes, but you don't know which kit was used, then just use "unknown" for `--barcode-kit`.
3- Providing a barcode description file (meaning using `--description`) will result in deleting any sequence assigned to a barcode no present in the barcode description file.
```commandline
usage: python basecall_nanopore.py [-h] -i /path/to/input_folder/ -o /path/to/output_folder/ [-s {minion,promethion}] [-c dna_r9.4.1_450bps_sup.cfg] [-f FLO-MIN106] [-l SQK-LSK109] [-b EXP-NBD104]
                                   [-d /path/to/barcode_description.tsv] [-r] [-t 16] -g "cuda:0" [-p 2] [-m 57] [-v]

Basecall Nanopore raw data to fastq.

options:
  -h, --help            show this help message and exit
  -i /path/to/input_folder/, --input /path/to/input_folder/
                        Folder that contains the fast5 files. Mandatory.
  -o /path/to/output_folder/, --output /path/to/output_folder/
                        Folder to hold the result files. Mandatory.
  -s {minion,promethion}, --sequencer {minion,promethion}
                        Sequencer used. "minion" includes all sequencers except "promethion". Optional. Default is "minion".
  -c dna_r9.4.1_450bps_sup.cfg, --config dna_r9.4.1_450bps_sup.cfg
                        Guppy config file. Typically found in "/opt/ont/guppy/data" ending with ".cfg". This is the prefered methods over choosing a "library kit/flowcell" combination. Both methods are uncompatible.
                        Optional.
  -f FLO-MIN106, --flowcell FLO-MIN106
                        Flowcell type used for sequencing. Optional.
  -l SQK-LSK109, --library-kit SQK-LSK109
                        Library kit used. Optional.
  -b EXP-NBD104, --barcode-kit EXP-NBD104
                        Barcoding kit used. Use "unknown" if you know barcodes were used, but do not know which kit. Not using this option will not perform barcode splitting. Optional
  -d /path/to/barcode_description.tsv, --description /path/to/barcode_description.tsv
                        Tab-separated file with two columns with barcode assignments. First column contains barcode names [barcode01, barcode02, etc.]. Second column contains sample name. Avoid using special character.
                        Sample file in data folder. Optional.
  -r, --recursive       Look for fast5 recursively. Useful if fast5 are in multiple sub-folders. Optional
  -t 16, --threads 16   Number of threads. Default is maximum available(16). Optional.
  -g "cuda:0", --gpu "cuda:0"
                        GPU device tp use. Typically "cuda:0" is one compatible graphics card is installed. Use "cuda:0 cuda:1" (including the quotes) to use two graphics cards. Default is "cuda:0". Mandatory.
  -p 2, --parallel 2    Number of samples to process in parallel for trimming and filtering. Default is 2. Optional.
  -m 57, --memory 57    Memory in GB. Default is 85% of total memory (57). Optional.
  -v, --version         show program's version number and exit
```

## Examples
Different scenario:
1- No barcodes, R9.4.1 flowcell, Super Accuracy basecalling using config file.
```commandline
python basecall_nanopore.py \
    -i /data/fast5 \
    -o /analyses/basecalled \
    --recursive \
    -c dna_r9.4.1_450bps_sup.cfg
    -t 4 \
    -g "cuda:0"
```
2- Barcodes, barcode description, R10.3 flowcell, Fast basecalling using config file.
```commandline
python basecall_nanopore.py \
    --input /data/fast5 \
    --output /analyses/basecalled \
    --recursive \
    --config dna_r10.3_450bps_fast.cfg \
    --barcode-kit EXP-NBD104 \
    --description /data/bc.txt \
    --threads 4 \
    --gpu "cuda:0"
```
3- Barcodes unknown, no barcode description, R9.4.1 flowcell. Using 2 GPUs.
```commandline
python basecall_nanopore.py \
    --input /data/fast5 \
    --output /analyses/basecalled \
    --recursive \
    --flowcell FLO-MIN106 \
    --library-kit SQK-LSK109
    --barcode-kit unknown \
    --threads 4 \
    --gpu "cuda:0 cuda:1"
```
4-  2 Barcode kits, barcode description,  R9.4.1 flowcell, Super Accuracy basecalling using config file.
```commandline
python basecall_nanopore.py \
    -i /data/fast5 \
    -o /analyses/basecalled \
    --recursive \
    -c dna_r9.4.1_450bps_sup.cfg
    --barcode-kit "EXP-NBD104 EXP-NBD114" \
    --description /data/bc.txt \
    -t 4 \
    --gpu "cuda:0"
    
```