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

n = size(A,1); 
m = size(B,2);

%% Define decision variables:

Y = sdpvar(m,n,'full');
Q = sdpvar(n,n);

% Define constraints:

Posi_constr = Q >= 0;
      
Orthant_constr = Q*(A') + A*Q + B*Y + (B*Y)' <= 0;

tot_constr = Posi_constr+Orthant_constr;

ops = sdpsettings('verbose',verbosity_level,'debug',debug_level, ...
                                             'solver',use_solver{1});

optim = optimize(tot_constr, [], ops);
if optim.problem ~= 0 && optim.problem ~= 4
    if optim.problem == 1
        error(['Infeasible Problem']);
    else
        error(['Error solving problem: ' optim.info ])
    end
else
    fprintf('\n\t\tController found: ');
    K = value(Y)/value(Q)
end