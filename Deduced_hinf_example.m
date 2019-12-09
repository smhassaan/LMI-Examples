
clc
clear
close all

% Set YALMIP Options

verbosity_level = 0;
debug_level = 0;
use_solver = 'sedumi';

%% User specified variable

% Select an arbitrary value to check if it bounds the H-inf norm of the
% system
gamma = 1;

%% System

A =[-1     1     0     1     0     1;
    -1    -2    -1     0     0     1;
     1     0    -2    -1     1     1;
    -1     1    -1    -2     0     0;
    -1    -1     1     1    -2    -1;
     0    -1     0     0    -1    -3];

B =[ 0 -1 -1; 
     0  0  0;
    -1  1  1;
    -1  0  0;
     0  0  1;
    -1  1  1];

C =[0     1     0    -1    -1    -1;
    0     0     0    -1     0     0;
    1     0     0     0    -1     0];

D=[0 0 0;
   0 0 0;
   0 0 0];

num_x = size(A,2);
num_w1 = size(B,2);
num_w2 = size(C,1);
num_u = num_w1;
num_y = num_w2;
num_z1 = num_y;
num_z2 = num_u;

B1 = [B zeros(size(B,1),num_w2)];
B2 = B;
C1 = [C; zeros(num_z2,size(C,2))];
C2 = C;

D11 = [D zeros(size(D,1),num_w2);
       zeros(num_z2,num_w1), zeros(num_z2,num_w2)];
D12 = [D; eye(num_z2)];
D21 = [D eye(size(D,1))];
D22 = D;

% Required matrices:
A = A;
B = B1;
C = C2;
D = D21;

n = size(A,1);
m = size(C,1);
r = size(B,2);

% Define required decision variables

X = sdpvar(num_x);
Omega = sdpvar(n+r,n+r,'full');

% Define terms for LMI problem

Theta = [zeros(n), X, zeros(n,m+n+r);
         X, -X, zeros(n,m+n+r);
         zeros(m,n+n), -gamma*eye(m), zeros(m,n+r);
         zeros(n,n+n+m), -X, zeros(n,r);
         zeros(r,n+n+m+n), gamma*eye(r)];
     
Phi = [-eye(n), A', C', eye(n), zeros(n,r);
       zeros(r,n) B', D', zeros(r,n) -gamma*eye(r)];

Psi = [eye(n), zeros(n,15+r);
       zeros(r,n+15), -gamma*eye(r)];
   
constr = [Theta+(Phi')*Omega*Psi+((Phi')*Omega*Psi)' <= -eps];

% Set YALMIP Options
ops = sdpsettings('verbose',verbosity_level,'debug',debug_level, ...
                                             'solver',use_solver);

% Solve feasibility problem
optim = optimize(constr, gamma, ops);

if optim.problem ~= 0 && optim.problem ~= 4
    if optim.problem == 1
        fprintf('Infeasible Problem, H-inf gain is greater than ');
        fprintf('the arbitrary value\n\n');
    else
        error(['Error solving problem: ' optim.info ])
    end
else
    gamma = (value(gamma));
    fprintf('H-inf bounds the system\n');
end
