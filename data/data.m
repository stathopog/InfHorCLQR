%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Infinite horizon CLQR using FBS
%
%  Data generation for toy system
%
%  Author: Giorgos Stathopoulos - 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Defining the problem parameters

% generate random system with two states and one input
[dat,model] = small_random;
dat.px = size(dat.Cx,1);
dat.pu = size(dat.Cu,1);
dat.pf = size(model.x.terminalSet.H,1);
load('X_INIT_toy.mat');
dat.x_init = X_INIT(:,323);

% horizon, perturbation on the initial state and stepsize
dat.N = 100; 
perturbation = [zeros(length(X_INIT),1) 0.05*(1-2*rand(length(X_INIT),dat.N-1))]';
perturbation = perturbation(:,1);
[dat] = upper_bound_beta(dat.N,model,dat); 
dat.rho = 1 / dat.beta;

% condense
con = [];
[bar,con,H,F1,G1,F2,G2] = dense_form_generation(dat);
bar.cxx = bar.cx-bar.Cx*bar.A*dat.x_init;
con.c = [];
for i = 1:dat.N
    con.c = [con.c; bar.cu((i-1)*dat.no.con.u+1:i*dat.no.con.u,:); bar.cxx((i-1)*dat.no.con.x+1:i*dat.no.con.x,:)];
end
con.C = full(con.C);

%% write matrix to file
delete sizes_data;
fi1 = fopen('sizes_data','w');
fprintf(fi1,'%u ',dat.nx); fprintf(fi1,'%u ',dat.nu);
fprintf(fi1,'%u ',dat.px); fprintf(fi1,'%u ',dat.pu); fprintf(fi1,'%u ',dat.pf); fprintf(fi1,'%u ',dat.N);
fprintf(fi1,'%6.6f ',dat.beta); fprintf(fi1,'%6.6f ',dat.w);
fclose(fi1);

delete xinit_data;
csvwrite('xinit_data.dat',dat.x_init)

delete A_data;
csvwrite('A_data.dat',dat.A) 

delete B_data;
csvwrite('B_data.dat',dat.B)

delete Q_data;
csvwrite('Q_data.dat',dat.Q)

delete R_data;
csvwrite('R_data.dat',dat.R)
 
delete Qbar;
csvwrite('Qbar.dat',bar.Q)

delete Rbar;
csvwrite('Rbar.dat',bar.R)

delete S_data;
csvwrite('S_data.dat',dat.S)

dat.M = dat.R+dat.B'*dat.S*dat.B;
dat.L = chol(dat.M,'lower');
delete M_data;
csvwrite('M_data.dat',dat.M)
delete L_data;
csvwrite('L_data.dat',dat.L)

delete K_data;
csvwrite('K_data.dat',dat.K)

delete H_data;
csvwrite('H_data.dat',dat.Hf)
delete hf_data;
csvwrite('hf_data.dat',dat.hf)

delete Cbar;
csvwrite('Cbar.dat',con.C)
delete dbar;
csvwrite('dbar.dat',con.c)

delete Bbar;
csvwrite('Bbar.dat',bar.B)

delete Abar;
csvwrite('Abar.dat',bar.A)

delete Gbar;
csvwrite('Gbar.dat',G1)

delete D_data;
csvwrite('D_data.dat',dat.Cu)

delete C_data;
csvwrite('C_data.dat',dat.Cx)

delete cu_data;
csvwrite('cu_data.dat',dat.cu)

delete cx_data;
csvwrite('cx_data.dat',dat.cx)


for iii = 1:dat.N
    dat.W((iii-1)*dat.no.con.all+1:iii*dat.no.con.all) = ones(dat.no.con.all,1)*dat.w^(iii-1);
end
delete W;
csvwrite('W.dat',dat.W')

delete perturb;
csvwrite('perturb.dat',perturbation);