# Running first simulations

This chapter will show you how to run the first simulation. There are multiple ways to run a simulation. Here is a brief description of how the running simulations are envisioned. 
1. Running within singularity - This is the simplest way to run the simulation. Initially, one start the `singularity` container[^WhatIsSig]. This makes debugging and code developement easier. Also, it may be suited when running the simulation on a personal machine. One can also parallelize the simulation using multiple threads. 
[^WhatIsSig]: The `singularity` container is a virtual machine that contains all the required software and libraries to run the simulation.
2. Running outside of singularity - This is suited when, one wants to run the simulation on multiple clusters. The MOBO algorithm suggests points, these design points along with the location of the singularity container is wrapped in a `SLURM` script to run on a `HPC` node. Again, One can also parallelize the simulation using multiple threads either inside or outside the singularity container. It is though more practical to run simulations outside of singularity and parallelize the simulation using multiple threads within the singularity when submitting to `SLURM` node. 

More information about parallelization can be found in [](content:headers:Parallelization)

Hence, to getting started with simulations, there are three main modes one can run the simulations currently supporting both `batch` mode and `interactive` mode. In the `batch` mode, the simulations can be run on a single thread or can be distributed into multiple threads. While in the `interactive` mode, the simulation is run on a single thread and a visualization of the detector geometry is saved. 

```{note}
`interactive` mode is useful for debugging and understanding the geometry alongside getting familiar to `Geant4`. It is **NOT** used for running simulations.
```

A script named `runRandomPoint.py` is provided in the `BoTorch4EIC` directory. This script is meant to be used as a tutorial script for running simulations in various modes. 

```{dropdown} runRandomPoint.py usage
    
    python runRandomPoint.py --config HolyGrail_mobo_emulator_cfg.json --HowToRun 1
    args : 
        --config -c   : Configuration JSON file to initialize variables for optimization
        --HowToRun -r : 0, [1], 2 
                    : 0 - Run simulations in batch mode within in singularity. User when running optimization on a node with parallel simulations
                    : 1 - Build the geometry and visualize it and save the .root file. This option used only when outside of singularity.
                    : 2 - Run simulations in single thread. But this runs the simulations outside of singularity
    
```

# Run interactive visualization

To better understand the geometry and the problem in hand, let us start with running an interactive visualization of a random design point. In order to run this interactive visualization, follow the set of instructions below.

1. Navigate to `BoTorch4EIC` directory and checkout to `dev` branch. You should be able to see `HolyGrail_mobo_emulator_cfg.json` file along with other files.

2. Open the file named `HolyGrail_mobo_emulator_cfg.json` and change the following variables accordingly.
    * `RUNNING_DIR` -- path where the simulation is to be run.
    * `SINGULARITY_IMAGE` -- path where the singularity image is located.
    * `MOUNT_PATH` -- the first element in the path is where the `cvmfs` is mounted. 
    * `src_dir` -- path to `BoTorch4EIC` directory.
More information about the configuration file can be found in the section [The Configuration `JSON` file](../wraperScript/ConfigurationJSONFile.md).

3. Run the script `python runRandomPoint.py -c HolyGrail_mobo_cfg.json -r 1`. Note that the basically, creates a random design point with all the constraints and visualizes the design. The output file `Support_Structure.root` file gets saved in the `RUNNING_DIR`. 

4. The file can be viewed in [JSRoot](https://jsroot.gsi.de/latest), which can used to open the `Support_Structure.root` file and open the geometry named `World`. Hide the Blackhole geometery. Select the background to open `controls` and select `clipping` on `x`. 

Here is an example of a geometry made from the script 

<div class="video-container">
    <iframe src="https://karthik18495.github.io/ECCE_GeomFiles/JsRoot633/?nobrowser&file=../geom/Support_Structure.root&item=World&opt=all;clipx;roty317,rotz351,zoom25" height="315" width="100%" allowfullscreen="1" frameborder="0">
    </iframe>
</div>


# Run simulaton in `batch` mode

Until now, the geometry has been visualization. Simulations have not yet run. So in order to run simulation of tracks, it is better to run them in quiet mode rather than visualizating and simulating the tracks. In order to run `10000` events of $\pi^{+}$ tracks, 

## Outside `singularity` simualtion

This is to illustrate, the capabilities of running a python script which should load singularity and run a bunch of commands and exit. Therefore, run the following command outside `singularity` container.

```tcsh
python runRandomPoint.py --config HolyGrail_mobo_emulator_cfg.json --HowToRun 1
```

The output of the file is then analyzed and a bunch of pdf files are generated. These files can be found in the `RUNNING_DIR` directory with a random `uuid` name. The PDFs are the `resolutions` extracted for detector response variables such as the momentum resolution ($\frac{\sigma_{p}}{p}$) and angular resolution ($\sigma_{\theta}$, $\sigma_{\phi}$).

## Inside `singularity` simulation

This is useful when one is using the singularity container and running simulations from within the singularity on a single machine. This is how the whole framework was designed initially back in 2021. 

```tcsh
python runRandomPoint.py --config HolyGrail_mobo_emulator_cfg.json --HowToRun 0
```

Once when the script is executed, the resolution values along with Costproxies are printed. 

