function [dat,model] = small_random()

% Generating a random system
dat.A = [1.1 2; 0 0.95]; 
dat.B = [0; 0.0787];
C = [-1 1];
dat.nx = size(dat.A,1); dat.nu = size(dat.B,2);
dat.no.o = 1;

dat.Q = 2*(C'*C); dat.R = 2*eye(dat.nu); dat.N = zeros(dat.nx,dat.nu);

dat.Cu = [eye(dat.nu); -eye(dat.nu)];
dat.Cx = [eye(dat.nx); -eye(dat.nx)];
dat.no.con.u = size(dat.Cu,1);
dat.no.con.x = size(dat.Cx,1);
dat.no.con.all = dat.no.con.u+dat.no.con.x;  

% Constraints and matrices
umax = 1; % small
umin = -umax;
xmax = 10; % small
xmin = -xmax;
dat.cu = [ones(dat.nu,1)*umax; -ones(dat.nu,1)*umin];
dat.cx = [ones(dat.nx,1)*xmax; -ones(dat.nx,1)*xmin];
dat.umax = umax;
dat.umin = umin;
dat.xmax = xmax;
dat.xmin = xmin;

%% Compute invariant set
[dat.K,~,~] = dlqr(dat.A,dat.B,0.5*dat.Q,0.5*dat.R,zeros(dat.nx,dat.nu));
H = [dat.Cx; -dat.Cu*dat.K]; h = [dat.cx; dat.cu];
Omega = Polyhedron(H,h);
% Compute the maximal positively invariant set
iii = 1;
while 1
	Omegaprev = Omega;
	F = Omega.H(:,1:end-1);
    f = Omega.H(:,end);
	% Compute the pre-set
	Omega = Polyhedron([F;F*(dat.A-dat.B*dat.K)],[f;f]);
    Omega.minHRep();
	if Omega == Omegaprev, break; end
	iii = iii + 1;
end


% MPT3 form and terminal cost and set computation
model = LTISystem('A', dat.A, 'B', dat.B);
model.u.min = umin;
model.u.max = umax;
model.setDomain('x',Polyhedron(dat.Cx,dat.cx));
model.x.penalty = QuadFunction(0.5*dat.Q,zeros(1,dat.nx),0);
model.u.penalty = QuadFunction(0.5*dat.R,zeros(1,dat.nu),0);

P = model.LQRPenalty;
% Omega.plot();
model.x.with('terminalPenalty');
model.x.with('terminalSet');
model.x.terminalPenalty = P;
model.x.terminalSet = Omega;
dat.S = 2 * P.H;
dat.Hf = model.x.terminalSet.H(:,1:end-1);
dat.hf = model.x.terminalSet.H(:,end);
[dat.K,~,~] = dlqr(dat.A,dat.B,dat.Q,dat.R,dat.N);