# BoTorch4EIC

A repo to perform baysian based pareto optimization using BoTorch for designing detectors for the Electron Ion Collider

# Portable simulator
## Installation of the singularity container

Steps to install the singularity container:

1. Make sure to have singularity installed on your system. More information can be found [here](https://sylabs.io/guides/3.5/user-guide/quick_start.html)
2. Down and install `ECCE` singularity environment. 
```tcsh
git clone https://github.com/ECCE-EIC/Singularity.git
cd Singularity/
./updatebuild.sh
```
3. Start the singularity container.
```tcsh
singularity shell -B /path/to/Singularity/cvmfs:/cvmfs,/any/other/path/to/bind:/work /path/to/Singularity/cvmfs/eic.opensciencegrid.org/singularity/rhic_sl7_ext.sif
```
4. Load the environment `source /cvmfs/eic.opensciencegrid.org/ecce/gcc-8.3/opt/fun4all/core/bin/ecce_setup.sh`

## Downloading necessary `src code`

1. Download the `src code` as
```tcsh
git clone https://github.com/karthik18495/BoTorch4EIC.git
```
2. navigate to `cd BoTorch4EIC`

3. Run the script `python runRandomPoint.py`. Note that the basically, creates a random design point (Note it does not check for "all" constraints) and visualizes the design.

4. `>> .L DisplayOn.C` -- to activate the visualization drivers required by Geant4

5. `>> PHG4Reco *g4 = QTGui();` --  to visualize the design in `QT`

6. Should be able to see the tracking detector overview

## Structure of the Design code 

1. We have two main files that drive the design simulation. The `.config` file which are necessary design parameter file that will be written dynamically by a wrapper script (`BO` suggested design points). Next is the `.setting` file which defines all the hyper parameters for the optimization. `HolyGrail.setting` is the one corresponding to the maximum number of design paramters. Find the details of it below. 

## The design problem -- Holy Grail 

We are adapting to solve Tracker for ECCE. We have a set of cylinderical detectors and a set of disk type detectors. 

* `nCylinderical layers` - `7`
* `nElectronDisks` - `7`
* `nHadronDisks` - `7`

|    Design Variable   	|                              Remarks                             	| Number of Design Parameters 	|
|:--------------------:	|:----------------------------------------------------------------:	|:---------------------------:	|
|    `Si_Cyl_Radius`   	|                 Radius of  cylindrical detectors                 	|             `7`             	|
|  `Si_Cyl_Thickness`  	|                Thickness of  cylindrical detectors               	|             `7`             	|
|   `Si_Cyl_eLength`   	|          Electron-side length of  cylindrical detectors          	|             `7`             	|
|   `Si_Cyl_hLength`   	|            Hadron-side length of cylindrical detectors           	|             `7`             	|
|    `Si_Cyl_pitch`    	|  pitch used for  resolution calculation of cylindrical detectors 	|             `7`             	|
|     `Si_EDisk_Z`     	|                  `z` position of  Electron disk                  	|             `7`             	|
| `Si_EDisk_thickness` 	|                      Electron disk thickness                     	|             `7`             	|
|   `Si_EDisk_pitch`   	| pitch used for  resolution calculation  of cylindrical detectors 	|             `7`             	|
|     `Si_HDisk_Z`     	|                   `z` position of  Hadron disk                   	|             `7`             	|
| `Si_HDisk_thickness` 	|                       Hadron disk thickness                      	|             `7`             	|
|   `Si_HDisk_pitch`   	| pitch used for  resolution calculation  of cylindrical detectors 	|             `7`             	|
|   `Magnetic Field`   	|            Central magnetic field (Using sPHENIX Map)            	|             `1`             	|
|  `Total Parameters`  	|                  Total number of  design params                  	|             `78`            	|

![Sample display of a random design point](Images/Detector-example.png "Sample display of a random design point")

## The Support struture

The tracking support is made of mainly `5` materials, this can also, be as a design parameter. The materials are 
* Aluminimum 
* Polythene
* Water
* PEEK
* Iron 

There are two cylindrical supports made of carbon fibre 
* One for the "vertex" layers
* One for the "sagitta" layers



## Status
* ~~Implement tracking support~~
* ~~Having Overlaps between tracking support and Cylinders/Disks when changing thickness~~
* Need to write down constraints and embed them
    * ~~The barrel layer lengths have to be in multiple of 15. THis is due to technological choice, ITS3~~
    * The FST Disks have to have a minimum of 10cm distance between them. (But double layered design not possible)
* Need to run a 1M event sample of `pi-` and `e-` to extract resolutions
    * This is to ensure "statistical error" on simulation. 
    * Also, perform a resource / time estimation on simulation of 1M points across various parallel cores.
* Need to write down documentation on these sections
    * ~~Introduction to Singularity container and how to use it~~
    * The tracker at EIC
        * Explain about the tracker in general
        * The Silicon technology
        * What are the scenarios not considered in this simulation
        * What is momentum resolution and how it is computed
    * The design problem
        * The various configurations
        * Design dimensions for these configurations
        * Caveats for the design problem
    * How to use the scripts in this repository
    * Resource / time estimation
* Start a `jupyter-book` for documentation and as a logbook

