function [bar,con,H,F1,G1,F2,G2] = dense_form_generation(dat)

% Generated the dense form MPC problem for a given horizon
bar.Cx = sparse(blkdiag(kron(eye(dat.N-1),dat.Cx),dat.Hf));  bar.Cu = sparse(kron(eye(dat.N),dat.Cu));
bar.Q =  blkdiag(kron(eye(dat.N-1),dat.Q), dat.S);   bar.R = kron(eye(dat.N),dat.R);
bar.A = dat.A;
for i = 2:dat.N
    bar.A = [bar.A; dat.A^i];
end
bar.B_coln = dat.B;
for i = 1:dat.N-1
    bar.B_coln = [bar.B_coln; dat.A^i*dat.B];
end
bar.B = [bar.B_coln zeros(dat.N*dat.nx,(dat.N-1)*dat.nu)];
for i = 1:dat.N-1
    bar.B(:,i*dat.nu+1:(i+1)*dat.nu) = [zeros(i*dat.nx,dat.nu); bar.B_coln(1:end-i*dat.nx,:)];
end
noise_bar_coln = 0*eye(dat.nx);
for i = 1:dat.N-1
    noise_bar_coln = [noise_bar_coln; dat.A^i];
end
bar.noise = [noise_bar_coln zeros(dat.N*dat.nx,(dat.N-1)*dat.nx)];
for i = 1:dat.N-1
    bar.noise(:,i*dat.nx+1:(i+1)*dat.nx) = [zeros(i*dat.nx,dat.nx); noise_bar_coln(1:end-i*dat.nx,:)];
end
H = bar.B'*bar.Q*bar.B + bar.R;
G1 = bar.A'*bar.Q*bar.B; G2 = bar.noise'*bar.Q*bar.B;
F1 = bar.A'*bar.Q*bar.A + dat.Q; F2 = bar.noise'*bar.Q*bar.noise;
bar.cx = [repmat(dat.cx,dat.N-1,1); dat.hf];  bar.cu = repmat(dat.cu,dat.N,1);
bar.Cxx = bar.Cx*bar.B;
con.C = []; 
for i = 1:dat.N
    con.C = [con.C; bar.Cu((i-1)*dat.no.con.u+1:i*dat.no.con.u,:); bar.Cxx((i-1)*dat.no.con.x+1:i*dat.no.con.x,:)];
end