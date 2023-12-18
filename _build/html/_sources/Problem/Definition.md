# Gentle introduction to the tracker

In this section, a brief non technical introduction about the tracker is outlined in detail. Trackers in Nuclear and Particle Physcics experiments measure the properties of charged particles. When charged particle passes through a material, it ionizes the enviroment around it triggering a signal. Stacking a bunch of layers of such material can be used to measure the trajectory of the particle.

The role of a tracking detector is briefly summarized below
1. It helps in reconstructing vertices. Vertex is a point where more than one particle emerge from the same point. 
2. In order to measure the momentum of the charged particle, often a magnetic field is used. The charged particle then bends in the magnetic field and the radius of curvature is used to measure the momentum of the particle. Techniques such as kalman fitering is used to measure the trajectory of the particle and hence deduce its momentum. 
3. Tracking detectors also helps in measuring finite particle lifetimes. 

## Caveats of tracking
1. Each layer of tracker can have an hit efficiency ($\epsilon$), which is the probablity of having a trigger when a charged particle passes through it. Hence, more number of layers for tracking will increase the accuracy of the measurement.
2. The material of the tracker can cause scattering of the charged particle. This can cause the trajectory of the particle to change and hence, the momentum of the particle to change. Hence, it is important to reduce the amount of material in the tracker.
3. The magnetic field can cause the charged particle to bend. For small values of momentum the radius of curvature is larger and for large values of momentum, the radius of curvature is smaller (more straight tracks). Hence, the radius of curvature have to be maximized in order to measure accurately the small values of momentum.

## Tracker coverage for the toy detector

Inspired by EIC generic tracker, similar definitions for phase space is adapted here. Quantity $\eta$ is called as pseudo-rapibity which is a measure of the angle $\theta$ in the lab frame. The pseudo-rapidity is defined as 

$$\eta = -\ln\Big[{\tan{(\frac{\theta}{2})}}\Big]$$

Also, particles originates from the vertex with momentum ranging as low as 0.1 GeV/c to as high as 100 GeV/c. Hence, the tracker is divided into three regions based on the pseudo-rapidity $\eta$ as shown in {numref}`eic-detector-coverage`. The Hadron endcap region is also called as the forward region, and the Electron endcap region is also called as the backward region. The region where there is a transition from the central barrel region to the end cap region is called as transition region.

```{figure} ../images/eic-detector-coverage.png
---
height: 300px
name: eic-detector-coverage
---

EIC Detector 1 $\eta$ coverage. Figure from {cite:p}`EIC_YellowReport`
```

A summary of phase space that will be used in this simulation is summarized in the table below 

|                            |  **Range** | **Bins** |
|:--------------------------:|:----------:|:--------:|
| **Momentum ($p$) [GeV/c]** |   1 - 30   |    30    |
|      **Eta ($\eta$)**      | -3.4 - 3.4 |     5    |

## References
```{bibliography}
```