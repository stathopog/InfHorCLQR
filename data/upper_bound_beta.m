function [dat] = upper_bound_beta(T,model,dat)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Computes an upper bound of the stepsize rho = 2 / beta
%
%  beta <= |C|^2 |H^-1|. 
%
%  Author: Giorgos Stathopoulos - 2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% To compute |C| = sigma(Cu) + Hinf(sys), choose weight w such that A
% is stable
if ( max(eig(model.A))>1 )
    dat.w = 1/(1.1*max(eig(model.A)));
    tempA = dat.w * model.A;
else
    dat.w = 1;
    tempA = model.A;
end
m = size(model.u.penalty.H,1);
sys = ss(tempA,model.B,dat.Cx,zeros(size(dat.Cx,1),m),1);
sigmas = sigma(sys);
for i = 1:m
    sigmas = max(sigmas);
end
Hinf = sigmas;

Cnorm = norm(dat.Cu,2) + Hinf;

dat.beta = Cnorm^2 * (1 / min(eig(model.u.penalty.H)));

for iii = 1:T
    dat.W((iii-1)*dat.no.con.all+1:iii*dat.no.con.all) = ones(dat.no.con.all,1)*dat.w^(iii-1);
end
            
            
    




