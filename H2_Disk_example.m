clc
clear
close all

set(0,'DefaultFigureWindowStyle','docked');
addpath(genpath('..\..\Toolkit'));

% Set YALMIP Options

verbosity_level = 0;
debug_level = 0;
use_solver = {'sedumi'; 'gurobi'; 'cplex'; 'scip'};

%% Define the system matrices

A = [-1 3 2; 0 1 0; 1 2 -1];
B = [1 0; 2 3; 1 1];
C = [1 1 0];

n = size(A,1); 
m = size(B,2);
p = size(C,1);

D = zeros(p,m);

% Define desired Disk region:

r0 = 4; q = 4;

%% Define decision variables:

gamma = sdpvar(1,1);
K = sdpvar(m,p,'full');

% Define constraints:

Disk_constr = [-r0*eye(n) A+B*K*C+q*eye(n);
          (A+B*K*C+q*eye(n))', -r0*eye(n)] <= 0;
      
H2_constr = [-gamma*eye(m), K; K', -gamma*eye(p)] <= 0;

tot_constr = Disk_constr + H2_constr;

ops = sdpsettings('verbose',verbosity_level,'debug',debug_level, ...
                                             'solver',use_solver{1});

optim = optimize(tot_constr, gamma, ops);
if optim.problem ~= 0 && optim.problem ~= 4
    if optim.problem == 1
        error(['Infeasible Problem']);
    else
        error(['Error solving problem: ' optim.info ])
    end
else
    fprintf('\n\t\tController found: ');
    K = value(K)
    
    fprintf('\n\tH2 norm bound: ');
    gamma = value(gamma)
end