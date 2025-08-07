#This repository is still under development. The classes "Fullerene", "Fullerene_new", and "Schemes" are implemented.
instances of the "Fullerene" class have all the physical characteristics of a fullerene: 
the speed of the center of mass, 
coordinates
kinetic momentum
rotation speed

the class "Fullerene_new" is needed to obtain shifted coordinates 
for calculating the coefficients of the R-K method.

The "Schemes" class stores attributes of the numerical method
for solving the system of equations of motion.
Such attributes include all coefficients of this method, the time step.

#Statement of the problem
Let there be n1 molecules of uncharged fullerene C60, consisting of n2 = n1·60 carbon atoms, in a representative cube with a side of 10 nm. Due to their large molecular radius, fullerenes have both translational and rotational motion components. We also assume that the bond length in fullerene between carbon atoms is fixed and is equal to 0.144±0.001 nm or 0.139±0.001, depending on the type of this bond. Van der Waals forces are present between the molecules, which are taken into account in the interaction potential.

<img width="176" height="50" alt="image" src="https://github.com/user-attachments/assets/2c6f6710-4aef-4445-8b50-97cf48c567b2" />

Due to its simplicity and a fairly accurate description of Van der Waltz interactions between neutral molecules, the Lennard-Jones potential with the following parameters was chosen as the interaction potential between carbon atoms composing fullerenes.:
σ=0.335 nm,ε=17.25*10^(-23) J.

The data on the initial velocities of translational and rotational motion, the position of the center of mass and the atoms of each fullerene are presented below.:

<img width="205" height="140" alt="image" src="https://github.com/user-attachments/assets/6e98c739-3f13-40a8-9141-88941f6c8a35" />

equations of motion that will be solved by the Runge-Kutta method

<img width="408" height="49" alt="image" src="https://github.com/user-attachments/assets/988b73fa-61cd-402f-88ff-1ccd104c6468" />


<img width="482" height="46" alt="image" src="https://github.com/user-attachments/assets/bc6a6137-f25f-40f5-8e82-9f9bedb7f417" />


<img width="336" height="47" alt="image" src="https://github.com/user-attachments/assets/85a66f09-bb29-4344-9704-1c357f777ca3" />

<img width="391" height="106" alt="image" src="https://github.com/user-attachments/assets/55a2e8df-78fe-4650-a25a-011094db872d" />



<img width="414" height="316" alt="image" src="https://github.com/user-attachments/assets/3d0e542e-f6d9-4613-ad42-cb9e2c8c4d9a" />

<img width="396" height="87" alt="image" src="https://github.com/user-attachments/assets/7c5942e7-a547-4d34-a837-b16249a5c2eb" />


<img width="566" height="180" alt="image" src="https://github.com/user-attachments/assets/ac3dee71-4590-4495-bf8c-4b4f75de9803" />


<img width="615" height="92" alt="image" src="https://github.com/user-attachments/assets/4de508bf-1a43-4f77-b6c9-eb5d511919c5" />


Ultimately, the fullerene equations of motion are solved in the following order:
1. Calculation of projections of external forces (6) in terms of the derivative of the interaction potential.
2. Finding the components of the velocity vector of the center of mass from equations (5).
3. Calculating the coordinates of the center of mass from equations (7) taking into account (8).
4. Finding the coordinates of all fullerene atoms using the equations (10)–(12).
5. Determination of the components of the inertia tensor by formulas (14)-(19).
6. Finding projections of the kinetic moment on the Oxyz axis from (24)-(26).
7. The solution of SLOUGH (20)-(22).

#Runge-Kutta method

<img width="250" height="44" alt="image" src="https://github.com/user-attachments/assets/775bd142-97f1-4636-99ed-06de04642980" />

<img width="310" height="158" alt="image" src="https://github.com/user-attachments/assets/a6a712b1-28cf-4994-a2f5-d4f2895e900a" />

<img width="328" height="57" alt="image" src="https://github.com/user-attachments/assets/d79aaf13-0bec-4d84-95a9-da037beb8c6c" />











