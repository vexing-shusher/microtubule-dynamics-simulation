# Simulation of stochastic dynamics of tubulin microtubules
My Bachelor's degree project (Moscow State University, Dept. of Physics): microtubule dynamics simulation.

# Disclaimer
1. This work is not (and was not) intended for publication in any peer-reviewed scientific media. The results presented here may turn out to be not scientifically significant. In essense, I just gave it a go.
2. At the time of working on this project, my programming skills were pretty poor. The code should not contain any bugs, however, I cannot guarantee that :)

# Introduction

## Importance
The construction of a detailed computational model of a microtubule will make it possible to carry out computer simulations of complex experiments (including optical ones) and, even more significantly, will facilitate their interpretation from a theoretical point of view.

## Microtubule properties
1. Dynamic instability
2. Polarity of the ends
3. Ability to generate force
4. Aging

## Model description
Two main coordinate systems are used in this project:

![Molecular coordinate system (black)](https://github.com/vexing-shusher/microtubule-dynamics-simulation/blob/main/materials/mtpic1.png)
![Monomer coordinate system (blue)](https://github.com/vexing-shusher/microtubule-dynamics-simulation/blob/main/materials/mtpic2.png)

The model pipeline can be described by the following diagram:

![Model diagram](https://github.com/vexing-shusher/microtubule-dynamics-simulation/blob/main/materials/prog.png)

In essense, it solves the following system of equations (using 4th order Runge-Kutta algorithm):

![Model equations](https://github.com/vexing-shusher/microtubule-dynamics-simulation/blob/main/materials/eq1.png)

The parameters of the model can be found inside the .cpp file, as well as in the following table:

![Model parameters](https://github.com/vexing-shusher/microtubule-dynamics-simulation/blob/main/materials/pars.png)

# Experiments

## Experiment 1: Einstein's formula for Brownian motion

It is known from the Einstein's formula for Brownian motion that a linear dependence of the dispersion on time should be observed for both translational and rotational coordinates (see equations below). 

![Einstein's formula](https://github.com/vexing-shusher/microtubule-dynamics-simulation/blob/main/materials/eq2.png)

The conducted experiment demonstrates that this is indeed true for the model presented here.

![Einstein's formula experiment](https://github.com/vexing-shusher/microtubule-dynamics-simulation/blob/main/materials/tabD.png)

