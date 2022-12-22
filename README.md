# Simulation of stochastic dynamics of tubulin microtubules
My Bachelor's degree project (Moscow State University, Dept. of Physics): microtubule dynamics simulation.
Advisor: Y.M. Romanovsky
Co-Advisor: N.B. Gudimchuk

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

![Molecular coordinate system (black)](https://github.com/vexing-shusher/microtubule-dynamics-simulation/blob/main/materials/mtpic1.png | width = 200)
![Monomer coordinate system (blue)](https://github.com/vexing-shusher/microtubule-dynamics-simulation/blob/main/materials/mtpic2.png){ width: 200px; }

The model pipeline can be described by the following diagram:

![Model diagram](https://github.com/vexing-shusher/microtubule-dynamics-simulation/blob/main/materials/prog.png){ width: 200px; }

In essense, it solves the following system of equations (using 4th order Runge-Kutta algorithm):

![Model equations](https://github.com/vexing-shusher/microtubule-dynamics-simulation/blob/main/materials/eq1.png){ width: 200px; }

The parameters of the model can be found inside the .cpp file, as well as in the following table:

![Model parameters](https://github.com/vexing-shusher/microtubule-dynamics-simulation/blob/main/materials/pars.png){ width: 200px; }

# Experiments

## Experiment 1: Einstein's formula for Brownian motion

It is known from the Einstein's formula for Brownian motion that a linear dependence of the dispersion on time should be observed for both translational and rotational coordinates (see equations below). 

![Einstein's formula](https://github.com/vexing-shusher/microtubule-dynamics-simulation/blob/main/materials/eq2.png){ width: 200px; }

The conducted experiment demonstrates that this is indeed true for the model presented here.

![Einstein's formula experiment](https://github.com/vexing-shusher/microtubule-dynamics-simulation/blob/main/materials/tabD.png){ width: 200px; }

## Experiment 2: Tetramer -- free bottom monomer

In this experiment, as well as in the following ones, the dynamics of a tetramer (a structure of 4 tubulin monomers) is studied. Two types of tetramers are considered: 3j6e (all monomers are bound to GTP) and 3j6f (2nd and 4th monomers are bound to GDP, simulating the microtubule protofilament after hydrolysis). Additionally, two boundry conditions are considered: free bottom monomer (i.e. freely floating tetramer) and fixed bottom monomer (i.e. the tetramer at the end of the microtubule).

Here are the equilibrium coordinates of the monomers in the tetramer (bottom monomer is free).

![Free ends1](https://github.com/vexing-shusher/microtubule-dynamics-simulation/blob/main/materials/fe_xz.png){ width: 200px; }
![Free ends2](https://github.com/vexing-shusher/microtubule-dynamics-simulation/blob/main/materials/fe_yz.png){ width: 200px; }

## Experiment 3: Tetramer -- fixed bottom monomer

Here is how the coordinates of the tetramer change with time if the bottom monomer is fixed.

### 3j6e (all-GTP monomers)
![Fixed ends1](https://github.com/vexing-shusher/microtubule-dynamics-simulation/blob/main/materials/3j6e_xz.png){ width: 200px; }
![Fixed ends2](https://github.com/vexing-shusher/microtubule-dynamics-simulation/blob/main/materials/3j6e_yz.png){ width: 200px; }

### 3j6f (GTP-GDP monomers)
![Fixed ends3](https://github.com/vexing-shusher/microtubule-dynamics-simulation/blob/main/materials/3j6f_xz.png){ width: 200px; }
![Fixed ends4](https://github.com/vexing-shusher/microtubule-dynamics-simulation/blob/main/materials/3j6f_yz.png){ width: 200px; }

### A microtubule fragment consisting of tetramers without lateral bounds
![Fixed ends5](https://github.com/vexing-shusher/microtubule-dynamics-simulation/blob/main/materials/mtnb.png){ width: 200px; }


## Experiment 4: Tetramer -- lateral bounds are active
Finally, here are the equilibrium conformational angles for the tetramers when the lateral bounds are active.

### 3j6e (all-GTP monomers)
![Fixed ends1](https://github.com/vexing-shusher/microtubule-dynamics-simulation/blob/main/materials/3j6e_bnb_xz.png){ width: 200px; }
![Fixed ends2](https://github.com/vexing-shusher/microtubule-dynamics-simulation/blob/main/materials/3j6e_bnb_yz.png){ width: 200px; }

### 3j6f (GTP-GDP monomers)
![Fixed ends3](https://github.com/vexing-shusher/microtubule-dynamics-simulation/blob/main/materials/3j6f_bnb_xz.png){ width: 200px; }
![Fixed ends4](https://github.com/vexing-shusher/microtubule-dynamics-simulation/blob/main/materials/3j6f_bnb_yz.png){ width: 200px; }

# Conclusion
1. The characteristic equilibrium time for the tetramer was determined to be 0.7 μs;
2. It was found that the equilibrium values of the conformational angles predicted by the model slightly differ from those included in it;
3. It was found that with the parametrization chosen in the work, the lateral bonds formed by four monomers are not enough to compensate for the bend in the radial plane;
4. It was found that the bending of the protofilament in the (y, z) plane remains virtually unchanged in the presence of lateral bonds.

# Final remarks
Even though the real scientific significance of this work is (in my opinion) questionable, it was an exciting journey I was happy to make. I would like to thank my advisor, Y.M. Romanovsky, as well as N.B. Gudimchuk, for their support and insights into this fascinating topic.


