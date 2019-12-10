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

A = [-0.5, 0, 1; 0 -2 0; 0 0 -5];
B = [5 3; 6 5; 7 6];
C = [2 -2 1];

n = size(A,1); 
m = size(B,2);
p = size(C,1);

D = zeros(p,m);

% Define desired strip region:

r1 = -5; r2 = -1;

%% Define decision variables:

gamma = sdpvar(1,1);
K = sdpvar(m,p,'full');

% Define constraints:

Strip_constr_low = A+B*K*C+(A+B*K*C)' >= 2*r1*eye(n);
Strip_constr_high = A+B*K*C+(A+B*K*C)' <= 2*r2*eye(n);
      
H2_constr = [-gamma*eye(m), K; K', -gamma*eye(p)] <= 0;

tot_constr = Strip_constr_low + Strip_constr_high + H2_constr;

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