The bulk of this repository is taken from Michael Spira's [HPAIR code](http://tiger.web.psi.ch/hpair/). My contribution was simply to write some code around it to perform parameter scans.

# HPAIR installation
 + Make sure you have the Fortran compiler ```gfortran``` installed
 + Download the compressed file stored [here](http://tiger.web.psi.ch/hpair/).
 + ```mkdir hpair``` and copy the above file to this folder
 + Decompress the file: ```tar xzvf hpair.tar.gz```
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
 + Check the configuration stored in ```hpair.in```:
     + For LO-only, set ```LOOP = 1```
     + For SM (```MODEL = 0```) without new resonances, use ```PROCESS = 0```
     + To change the trilinear coupling, one can modify ```C_3 =  1.0D0```

# Run this wrapper
If you want to vary the "C_3" parameter and produce a ```results.txt``` files with the resulting cross-sections, run:

```shell
go run create_files.go m1 0 1 2 3
```

where the following happens:

+ HPAIR input files are stored under ```inputs/``` (the suffix ```m``` denotes negative numbers)
+ HPAIR is run over those inputs
+ the standard HPAIR output files will be stored under ```outputs/```.
