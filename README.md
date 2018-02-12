# Infinite Horizon CLQR
C++ code for solving the infinite horizon constrained LQR problem with the dual proximal gradient method. Comparison with linear MPC.

## What it does
The regulation of a 2-state toy-system to the origin while respecting state and input constraints is solved using (i) the infinite horizon constrained LQR approach (**CLQR**) as presented in the corresponding papers (see References) and (ii) using conventional model predictive control (**MPC**). The problem is solved in closed loop and a perturbation to the initial state is applied at the beginning of each solve.

## How to use it
There are three folders, namely *data*, *main_CLQR* and *main_MPC*. 
The *data* folder contains the system's data (matrices & vectors) as well as the initial state denoted by 'xinit_data.dat'. 
The *main_CLQR* folder contains the compiled code for the CLQR solution. Execute `./clqr` from a terminal, located in `main_CLQR/build/src/CLQR`
The *main_MPC* folder contains the compiled code for the MPC solution. Execute `./mpc` from a terminal, located in `main_MPC/build/src/MPC`
After having performed any change to the codes (e.g., to the CLQR code), run the following commands:
```
cd main_CLQR/build
make
cd main_CLQR/build/src/CLQR
./clqr
```
## References
1. G. Stathopoulos, M. Korda and C. N. Jones, "Solving the infinite-horizon constrained {LQR} problem using splitting techniques", 19th IFAC World Congress 2014, (https://infoscience.epfl.ch/record/197368?ln=en)
2. G. Stathopoulos, M. Korda and C. N. Jones, "Solving the Infinite-Horizon Constrained LQR Problem Using Accelerated Dual Proximal Methods ", IEEE Transactions on Automatic Control, 62, 4, 1752-1767, 2016 (https://infoscience.epfl.ch/record/228793?ln=en)