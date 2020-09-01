### Prerequisites

#### Mac

``` bash
brew install openssl
export C_INCLUDE_PATH=${C_INCLUDE_PATH}:/usr/local/Cellar/openssl/your_version/include
export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/usr/local/Cellar/openssl/your_version/lib/"
```

#### Linux
``` bash
sudo apt-get install openssl
```

### Install the module using:

``` bash
sudo python setup.py install
# or 
python setup.py install --user
```

#### Install conda python 3.6

Follow the instruction: https://conda.io/miniconda.html to install python 3.6

``` bash
wget -c https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
```

#### Install the dependencies

``` bash
bash install.sh
```


#### Start to use

Activate the environment

``` bash
source activate lisa
```

##### Input a gene list

The size should be at least 10 genes, no more than 3000, it is too many.

without l1 regularization selected DNase-seq sample peak region to limit motif region,

```bash
time THEANO_FLAGS='mode=FAST_RUN,device=cpu,floatX=float32,openmp=True' OMP_NUM_THREADS=8 lisa logit --gene AR.symbol --tf no --name AR --species hg38 -O ${HOME}/public_html/lisa/ARfolder

```

with l1 regularization selected DNase-seq sample peak region to limit motif region,

```bash
time THEANO_FLAGS='mode=FAST_RUN,device=cpu,floatX=float32,openmp=True' OMP_NUM_THREADS=8 lisa logit --gene AR.symbol --tf no --name AR --species hg38 -O ${HOME}/public_html/lisa/ARfolder --DNase
```

with ChIP-seq peak to validate the cis-element prediction.
`AR_ChIP-seq_peak` is made through intersection between genome wide 1kb windows and AR peak bed, return the index of the genome window, then use `np.load` and `np.save` to save as numpy binary array file.

```bash
time THEANO_FLAGS='mode=FAST_RUN,device=cpu,floatX=float32,openmp=True' OMP_NUM_THREADS=8 lisa logit --gene AR.symbol --tf AR_ChIP-seq_peak --name AR --species hg38 -O ${HOME}/public_html/lisa/ARfolder --DNase
```

By default, LISA will run all 10 marks, for selecting marks, use below:

```bash
time THEANO_FLAGS='mode=FAST_RUN,device=cpu,floatX=float32,openmp=True' OMP_NUM_THREADS=8 lisa logit --gene AR.symbol --tf AR_ChIP-seq_peak --name AR --species hg38 -O ${HOME}/public_html/lisa/ARfolder --DNase --histone "H3K27ac,DNase,H3K4me3,H3K27me3,ATAC-seq"

time THEANO_FLAGS='mode=FAST_RUN,device=cpu,floatX=float32,openmp=True' OMP_NUM_THREADS=8 lisa logit --gene AR.symbol --tf AR_ChIP-seq_peak --name AR --species hg38 -O ${HOME}/public_html/lisa/ARfolder --DNase --histone DNase
```

##### Input fastq and gene list

`lisa fastq` mode needs to fill the configuration files with bwa genome index `bwa_index`.

Remember, the `--name` is the prefix for new-generated HDF5 sample id, avoid using numbers and try to make them as unique as possible, such as "LnCaP_Study_2017".

``` bash
mkdir -p input_fastq_folder
mkdir -p input_fastq_folder/H3K4me3
mkdir -p input_fastq_folder/H3K27me3
mkdir -p input_fastq_folder/H3K27ac
# then cp the corresponding fastq files with suffix .fastq, .fastq.gz, .fq or .fq.gz into the sub-directory
```

The `input_fastq_folder` contains sub-directories, such as: H3K27me3,H3K4me3,H3K27ac. Under each of the sub-directory, there are single-ended ChIP-seq fastq files. Then, specify "--histone H3K27me3,H3K4me3,H3K27ac". The folder names and `--histone` option should be `H3K27ac,DNase,H3K4me1,H3K4me3,H3K4me2,H3K27me3,H3K36me3,H3K9me3,ATAC-seq,H3K9ac` or its subset. Other factors are not supported yet.

Then, run the `fastq` mode:
``` bash
lisa fastq --fastq input_fastq_folder --name LnCaP_Study_2017 --gene genes.txt --tf no -O output_html_folder --histone H3K4me3,H3K27me3,H3K27ac
snakemake -j 4 --use-conda
```


##### with in-house HDF5 to run other gene list

``` bash
lisa logit --histone H3K4me3 -O output_html --name test --gene gene.list --tf no --species hg38 --additional_h5 H3K4me3.h5 2>>{log}"
```

##### External filter regions (under experiment)

For `fastq` and `logit` mode, optionally, specify a region to filter chip-seq peak or motif regions.

``` bash
--EXPANDDNase new_regions.bed
```

###### TroubleShooting
Sometimes pandas would throw out HDFstore error, try:

``` bash
conda install --upgrade pandas pytables h5py
or 
pip install --upgrade pandas tables h5py
```

One possible problem is about mismatch between hdf5 header and library when using pandas.read_hdf, which calls pytables. Set the env of *HDF5_DIR* and the *LD_LIBRARY_PATH* to the same hdf5 in *.bashrc* as follows, or just set nothing, let the system decide..

``` bash
export HDF5_DIR=hdf5-x.x.x/hdf5
export LD_LIBRARY_PATH=hdf5-x.x.x/local/lib:$LD_LIBRARY_PATH
```
