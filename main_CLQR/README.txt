{\rtf1\ansi\ansicpg1252\cocoartf1187\cocoasubrtf400
{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
\paperw11900\paperh16840\margl1440\margr1440\vieww19940\viewh17460\viewkind0
\pard\tx566\tx1133\tx1700\tx2267\tx2834\tx3401\tx3968\tx4535\tx5102\tx5669\tx6236\tx6803\pardirnatural

\f0\fs24 \cf0 README for:\
PCSC Project: Constrained linear quadratic regulator\
===========================================================\
Quickstart:\
There are two folders, the main one 'clqr' and a 'Testsuite'\
\
Testsuite contains small examples that demonstrate the functionality of the methods.\
To run the examples go to build/opt/\
\
make\
\
then go to build/opt/Testsuite\
\
make\
\
and execute with \
\
./testing\
\
clqr contains the main algorithm.\
To run it go to build/opt/\
\
make\
\
then go to build/opt/Sources/ProjectCLQR/\
\
make\
\
and execute with \
\
./clqr\
\
IMPORTANT: The paths to read the data in main.c, function ReadInMatVecs (lines 149-224) has to be altered before execution. The data exist in the Matlab_files folder.\
\
The algorithm should terminate and return the optimal control sequence u and a hitting time t.\
\
At convergence, the control sequence elements should be almost inside the unit box [-1,1].\
Decreasing the tolerance from line 19 of the code should result to more iterations and 'more constrained' u values.\
}