# the-Runge-Kutta-method-of-the-4th-order
This repository considers one of the molecular dynamics problems in the following statement:
A cube with a side of 5 nanometers, n1 particles of gas of the 1st type and n2 particles of gas of the second type are given, their initial linear velocities and positions are determined. It is necessary to simulate the evolution of the system in 1 nanosecond. In this problem, we took 16 xenon atoms and 16 helium atoms.

 A two-parameter potential, the Lennard-Jones potential, was chosen as the interatomic interaction potential. Note that the following assumptions were made when solving this problem: 
#### 1. Collisions of gas molecules are absolutely elastic. 
#### 2. Particle rotation can be neglected.

Also note that our goal is to obtain the temperature distribution of this mixture at each moment of the simulation. The calculation error is controlled using the energy conservation law (in the code it is denoted by eps)
Gases of different types in the code are created through the class "Hu". This repository contains the header file "Hu.h" and "main.cpp"
