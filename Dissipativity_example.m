clc
clear
close all

set(0,'DefaultFigureWindowStyle','docked');
addpath(genpath('..\Toolkit'));

% Set YALMIP Options

verbosity_level = 0;
debug_level = 0;
use_solver = {'sedumi'; 'gurobi'; 'cplex'; 'scip'};

%% Define the system matrices

A = [-1.3410  0.9933 0  -0.1689  -0.2518
     43.2230 -0.8693 0 -17.2510  -1.5766
      1.3410  0.0067 0   0.1689   0.2518
           0       0 0 -20.0000        0
           0       0 0        0 -20.0000 ];
       
B = [0 0; 0 0; 0 0; 20 0; 0 20];

C = [eye(3), zeros(3,2)];

n = size(A,1); 
m = size(B,2);
p = size(C,1);

D = zeros(p,m);

%% Define variables and feasibility constraint

P = sdpvar(n,n);
Q = sdpvar(m+p,m+p);

positive_constr = [P>=eps];

dissipativity_constr = [(A')*P+P*A, P*B; (B')*P, zeros(m)] ...
         - ([C, D; zeros(m,n), eye(m)]')*Q*[C, D; zeros(m,n), eye(m)] <= 0;
                    
tot_constr = positive_constr + dissipativity_constr;

%% Solve for feasibility
ops = sdpsettings('verbose',verbosity_level,'debug',debug_level, ...
                                             'solver',use_solver{1});

optim = optimize(tot_constr, [], ops);
if optim.problem ~= 0 && optim.problem ~= 4
    if optim.problem == 1
        fprintf('Infeasible Problem\n\n');
    else
        error(['Error solving problem: ' optim.info ])
    end
else
    fprintf('System Dissipative with Q:\n')
    Q = value(Q)
end