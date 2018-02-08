# InfHorCLQR
C++ code for solving the infinite horizon constrained LQR problem with the dual proximal gradient method. Comparison with linear MPC.

## Further details
The regulation of a 2-state toy-system to the origin while respecting state and input constraints is solved using (i) the infinite horizon constrained LQR approach (**CLQR**) as presented in[]() and (ii) using conventional model predictive control (**MPC**). The problem is solved in closed loop and a perturbation in the initial state is applied in the beginning of each solve.

## Instructions
There are three folders, namely *data*, *main_CLQR* and *main_MPC*. 
The *data* folder contains the system's data (matrices & vectors) as well as the initial state denoted by 'xinit_data.dat'. 
The *main_CLQR* folder contains the compiled code for the CLQR solution. Execute `./clqr` from a terminal, located in `main_CLQR/build/src/CLQR`
The *main_MPC* folder contains the compiled code for the MPC solution. Execute `./mpc` from a terminal, located in `main_MPC/build/src/MPC`
After performing any change in the codes (e.g., to the CLQR code), run the following commands:
```
cd to main_CLQR/build
make
cd to main_CLQR/build/src/CLQR
./clqr
```
