# How are `NEvents` run in parallel

The logic behind how parallel computations are run

(content:headers:Parallelization)=
# Types of parallelizations 

There are two level of parallelizations. `SLURM` submission can be way too time consuming depending on the queue time. So should discuss if this should be considered.

# Drawbacks of parallelizations

Mainly, bottleneck from `Geant4` simulations. The libraries have to load, and this takes a fixed amount of time. Also, if large number of Events are simulated at once, the probablitly of job failing will increase and hence, risk of loosing the entire simulation.