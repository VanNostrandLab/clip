# clip

A shortcut for using eCLIP pipeline to process CLIP data on HPC Cluster.

## Installation
Clone this repository, cd into clip directory and then run `./install.sh`:

```shell script
$ git clone https://github.com/VanNostrandLab/clip.git
$ cd clip
$ ./install.sh
```

Note: in case you got `permission denied: ./install.sh` change the mode of `install.sh` 
then try it again:

```shell script
$ chmod +x install.sh
$ ./install.sh
```
If everything went well, the `install.sh` will take care of the virtual environment 
set up and install all required packages and components automatically. The script 
`install.sh` will be moved to directory `source` and three new directories: venv, eclip and 
merge_peak will be created. Also, an executable script `clip` will be created inside 
current directory. 

## Usage:

### For Van Nonstrand Lab users
The package has been installed in the following directory:
`/storage/vannostrand/software/clip`

Please `cd` into the directory and call `./clip -h`
or call `/storage/vannostrand/software/clip/clip -h` 
to see its usage.

### For other users
Once `clip` was successfully installed, it will tell you how to display its
usage at the very end. In case you did not noticed that, you can always check
its usage by yourself.

```shell script
$ ./clip -h
usage: clip [-h] [-o OUTDIR] [-j JOBNAME] [-e EMAIL] [-s SCHEDULER] [-t TIME] [-m MEMORY] 
            [-n CORES] [--debug] MANIFEST

A shortcut for using eCLIP pipeline to process CLIP data on HPC Cluster.

positional arguments:
  MANIFEST      Manifest YAML or JSON file describing paths of your dataset.

optional arguments:
  -h, --help    show this help message and exit
  -o OUTDIR     Path of the base output directory, default: the manifest file's parent directory.
  -j JOBNAME    Name of your job, default: eCLIP
  -e EMAIL      Email address for notify you the status (start, end, and abort) of you job.
  -s SCHEDULER  Name of the scheduler on your cluster, e.g., PBS (or qsub) or SBATCH (or slurm).
  -t TIME       Time (in integer hours) needed for your job, default: 32.
  -m MEMORY     Amount of memory (in GB) for all cores needed for your job, default: 32.
  -n CORES      Number of cores needed for your job, default: 8.
  --debug       Run the analysis in debug mode and keep cache files.
```

Note: in case you are no logger in the installation directory, you need to call `clip` using 
its relative path or absolute path.
