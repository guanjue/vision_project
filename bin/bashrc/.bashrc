# .bashrc                                  

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi
#module load python/2.7.11
module load python/2.7.8.ucs4
# Get the aliases and functions                  
# User specific environment and startup programs 
#PATH=$PATH:/gpfs/home/gzx103/scratch/guanjue/ROTATION/Peak-Shape_Calling/wd/pc/DNAshape/bin/
PATH=$PATH:~/work/bin/homer/bin
PATH=$PATH:$HOME/bin
PATH=$PATH:~/work/bin/fastx_toolkit/bin
#PATH=$PATH:/gpfs/home/gzx103/work/0606hg19data
#PATH=$PATH:/gpfs/home/gzx103/work/bin
PATH=$PATH:/gpfs/home/gzx103/work/bin/wavelets/bin
PATH=$PATH:~/work/bin/bedtools2/bin
#PATH=$PATH:/gpfs/home/gzx103/work/bin/Cmake/cmake-3.0.0/bin
PATH=$PATH:/gpfs/home/gzx103/work/bin/bin1
PATH=$PATH:/gpfs/home/gzx103/meme/bin
#PATH=$PATH:/gpfs/home/gzx103/group/software/bowtie-1.0.1
#PATH=$PATH:/gpfs/home/gzx103/group/software/samtools-0.1.18
PATH=$PATH:/gpfs/home/gzx103/work/bin/cluster_gene/bin
#PATH=$PATH:~/group/software/bowtie-1.0.1
PATH=$PATH:/gpfs/home/gzx103/work/bin/R-3.1.2/bin
#PATH=$PATH:~/group/projects/guanjue/bin/homer/bin
#PATH=$PATH:~/group/projects/guanjue/bin/blat/userApps/bin
#PATH=$PATH:~/group/projects/guanjue/bin/Ghostscript/ghostscript-9.16-linux-x86_64
#PATH=$PATH:~/group/projects/guanjue/bin/weblogo/weblogo
#PATH=$PATH:~/group/projects/guanjue/bin/samtools/samtools-1.2
#PATH=$PATH:~/group/projects/guanjue/bin/epigram-pipeline-0.004/bin
#PATH=$PATH:~/group/projects/guanjue/bin/epigram-pipeline-0.004/bin/epigram-0.002
PATH=$PATH:~/group/software/samtools
PATH=$PATH:~/work/bin/Cython-0.22
PATH=$PATH:~/work/bin/HTSeq-0.6.1/scripts
#PATH=$PATH:~/work/bin/SRA/sratoolkit.2.3.5-2-ubuntu64/bin
PATH=$PATH:~/src/bowtie2-2.2.4/
PATH=$PATH:~/work/bin/dreg/dREG/
PATH=$PATH:~/src/sratoolkit.2.3.5-2-ubuntu64/bin
PATH=$PATH:~/work/bin/bwa/bwa-0.7.10
PATH=$PATH:~/work/bin/macs/bin
PATH=$PATH:~/work/bin/PeakSplitter_Cpp/PeakSplitter_Linux64
PATH=$PATH:~/src/FastQC
PATH=$PATH:~/scratch/cegr_tools/cegr-tools/
#PATH=$PATH:~/work/bin/glibc
THEANO_FLAGS='cuda.root=/gpfs/apps/cuda-rhel6/cuda/7.0/,device=gpu,floatX=float32'
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/gpfs/apps/cuda-rhel6/cuda/7.0/lib/
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/gpfs/apps/cuda-rhel6/cuda/7.0/lib64/
#CUDA_ROOT=/gpfs/apps/cuda-rhel6/cuda/7.0/bin/
#THEANO_FLAGS='cuda.root=/usr/global/cuda/7.0/bin/,device=gpu,floatX=float32'

export PERL5LIB=/gpfs/home/gzx103/meme/lib/perl
export R_LIBS="/gpfs/work/g/gzx103/bin/R-3.1.2/lib64/R/library"
export PATH="$HOME/.linuxbrew/bin:$PATH"
export MANPATH="$HOME/.linuxbrew/share/man:$MANPATH"
export INFOPATH="$HOME/.linuxbrew/share/info:$INFOPATH"

export PATH

PS1='\w$ '

#PYTHONPATH=${PYTHONPATH}:~/work/bin/python2.7/Python2.7.8/lib/python2.7/site-packages
#PYTHONPATH=${PYTHONPATH}:/usr/global/python/2.7.8/lib/python2.7/site-packages
#PYTHONPATH=${PYTHONPATH}:~/work/bin/macs/lib/python2.7/site-packages
#PYTHONPATH=${PYTHONPATH}:/usr/global/python/2.7.11/bin/
#PYTHONPATH=${PYTHONPATH}:/gpfs/apps/cuda-rhel6/python/2.7.11/bin/lib/python2.7/site-packages
#PYTHONPATH=${PYTHONPATH}:~/work/bin/python2.7/Python2.7.11/lib
export PYTHONPATH

export SWIG_LIB=~/work/bin/swig/swig-3.0.8/Lib


#alias python="~/work/bin/python2.7/Python2.7.8/bin/python2.7"
alias Rscript="~/work/bin/R-3.1.2/bin/Rscript"
alias R="~/work/bin/R-3.1.2/bin/R"
alias awk="awk -F '\t' -v OFS='\t'"  # split and separate by tabs
alias c="clear"
alias qstatgj="qstat -u gzx103"

alias bridge="scp gxiang@bridges.psc.xsede.org"
#alias fseq="~/group/projects/guanjue/bin/Fseq/F-seq/dist~/fseq/bin/fseq"
#alias pip="~/work/bin/python2.7/Python2.7.8/lib/python2.7/site-packages/pip"
if [ -x /usr/bin/dircolors ]; then
    eval "`dircolors -b /gpfs/group/mahony/scripts/env/dir_colors`"
    alias ls='ls --color=auto'
fi


#if [ -f /gpfs/group/mahony/scripts/env/bashrc ]; then
#	source /gpfs/group/mahony/scripts/env/bashrc
#fi


cphere ~/code
export CLASSPATH=${CLASSPATH}:~/code/passwords
export CLASSPATH=${CLASSPATH}:~/code/seqcode-core/build/classes
export CLASSPATH=${CLASSPATH}:~/code/seqcode-dev/build/classes


#chmod -R o-rwx /gpfs/group/mahony/ 2>/dev/null
