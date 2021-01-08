# clip

A wrapper for running eCLIP and merge_peak (IDR) pipelines from command line. 

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
$ ./install
```
If everything went well, the `install.sh` will take care of the virtual environment 
set up and install all required packages and components automatically. The source code 
scripts (`install.sh` and `clip.py`) will be moved to directory `source` and a executable 
script `clip` will be created inside current directory. 

## Usage:
Once `clip` was successfully installed, it will tell you how to display its 
usage at the very end. In case you did not noticed that, you can always check 
its usage by yourself.

```shell script
$ ./clip -h
```

Note: in case you are no logger in the installation directory, you need to call `clip` using 
its relative path or absolute path.