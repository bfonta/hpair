The bulk of this repository is taken from Michael Spira's [HPAIR code](http://tiger.web.psi.ch/hpair/). 

My contribution was to create a wrapper script to perform parameter scans. Additionally I added a very basic I/O functionality to enable passing the names of input and out files to HPAIR's executable.

# HPAIR installation
 + Make sure you have the Fortran compiler ```gfortran``` installed
 + Clone this repository and ```cd hpair```
     + Alternatively, to get the original version of HPAIR, download the compressed file stored [here](http://tiger.web.psi.ch/hpair/), do a ```mkdir hpair```, copy the downloaded file to this folder and decompress it: ```tar xzvf hpair.tar.gz```.
 + Install LHAPDF following [these instructions](https://lhapdf.hepforge.org/install.html). I personanly used:
```shell
wget https://lhapdf.hepforge.org/downloads/?f=LHAPDF-6.5.4.tar.gz -O LHAPDF-6.5.4.tar.gz
tar xf LHAPDF-6.5.4.tar.gz
cd LHAPDF-6.5.4
sudo mkdir /opt/lhapdf
sudo ./configure --prefix=/opt/lhapdf
make
make install
```
 + Add the required LHA PDF set with ```sudo wget http://lhapdfsets.web.cern.ch/lhapdfsets/current/PDF4LHC15_nlo_100.tar.gz -O- | sudo tar xz -C /opt/lhapdf/share/LHAPDF``` (possibly adapting the paths if you installed ```lhapdf``` somewhere else)
 + Compile HPAIR: ```make && make install```
 + Check the configuration stored in ```hpair.template.in```:
     + For LO-only, set ```LOOP = 1```
     + For SM (```MODEL = 0```) without new resonances, use ```PROCESS = 0```
     + To change the trilinear coupling, one can modify ```C_3 =  1.0D0```

# Run this wrapper
You'll need to install the language ```go```. 
In order to vary the "C_3" parameter (Higgs trilinear coupling) and produce a ```results.txt``` files with the resulting cross-sections, run:

```shell
go run parameter_scan.go m1 0 1 2 3
```

where ```go```'s native parallelism is exploited and the following happens:

+ 5 HPAIR input files are stored under ```inputs/``` (the suffix ```m``` denotes negative numbers)
+ HPAIR is run over those inputs
+ the standard HPAIR output files (also 5 for this example) will be stored under ```outputs/```
+ a summary of the results is stored in ```outputs/results.csv``` for easier access (in the same order as provided by the user) with the following structure:
    + 1st line: values of the born cross-section
    + 2nd line: errors of the born cross-section
    + 3rd line: values of the NLO cross-section
    + 4th line: errors of the NLO cross-section
+ the values above will have different meanings depending on the ```hpair.template.in``` file
