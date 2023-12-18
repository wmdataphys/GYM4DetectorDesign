# Getting Started

In order to get started, It is better to have a linux based operating system since the framework that is being used, natively runs on most of the linux flavours. However, procedures mentioned in this documentation must apply to all operating system flavours through the mode of Virtual Machines and singularity containers. 

The portable simulator makes uses of Singularity container that was used by ECCE protocolaboration. The Singularity has been modified such that it is relative lighter to better make the singularity more mobile. The singularity has built in Fun4All framwork. A `Geant4` based framework. The portable simulator has 2 mains components, The ECCE Singularity and the macros that steers the design. 

Please report any issues or bugs to the [GitHub repository](https://github.com/karthik18495/BoTorch4EIC/issues) or click on the `GitHub` tab on the top right corner of the page.

# Singularity and source code

Follow the steps below to install the singularity container:

1. Make sure to have singularity installed on your system. More information can be found [here](https://sylabs.io/guides/3.5/user-guide/quick_start.html)
2. Download and install `ECCE` singularity environment. 
```tcsh
git clone https://github.com/ECCE-EIC/Singularity.git
cd Singularity/
./updatebuild.sh
```
3. Once downloaded, one can start the singularity container using the following command.
```tcsh
singularity shell -B /path/to/Singularity/cvmfs:/cvmfs,/any/other/path/to/bind:/work /path/to/Singularity/cvmfs/eic.opensciencegrid.org/singularity/rhic_sl7_ext.sif
```
4. Load the environment `source /cvmfs/eic.opensciencegrid.org/ecce/gcc-8.3/opt/fun4all/core/bin/ecce_setup.sh`. If this does not work try with `source /cvmfs/eic.opensciencegrid.org/ecce/gcc-8.3/opt/fun4all/core/bin/ecce_setup.sh - n new`.

5. Try to run the following command to check if the singularity is working.
```tcsh
root-config --version
python -c "import ROOT"
```

```{note}
The last two steps to load singularity is only to check if the singularity is working. Once the singularity is working, please exit singularity and follow the next steps outlined here.
```

6. Download the `src code` for the portable simulator.
```tcsh
git clone https://github.com/karthik18495/BoTorch4EIC.git
```
7. Navigate to `cd BoTorch4EIC` and check out the branch `git checkout dev`


Great! Now you are ready to run the simulation.



