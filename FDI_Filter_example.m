% This example is based on the scenario given in the referenced paper

clc
clear
close all

set(0,'DefaultFigureWindowStyle','docked');
addpath(genpath('..\..\Toolkit'));

% Set YALMIP Options

verbosity_level = 0;
debug_level = 0;
use_solver = {'sedumi'; 'gurobi'; 'cplex'; 'scip'};

%% Given Physical Parameters

Ms = 6e3 * eye(4); Ms(4,4) = 360;

n_p_1 = size(Ms,2);

Tgs = -Ms*ones(n_p_1,1);

Tu = [zeros(n_p_1-2,1); 1];

Tn = Tu';

Tus = [-Tu; 1];

C = 1e3 * [ 12.4, -5.16,     0;
           -5.16,  12.4, -4.59;
               0, -4.59,   7.2];
           
K = 1e6 * [ 3.4, -1.8,    0;
           -1.8,  3.4, -1.6;
              0, -1.6,  1.6];
          
ka = 18819; ca = 365.3910;

Ks = [K+Tu*ka*Tn, -Tu*ka;
         -ka*Tn,     ka];
          
Cs = [C+Tu*ca*Tn, -Tu*ca;
         -ca*Tn,     ca];

%% System matrices:

As = [zeros(n_p_1), eye(n_p_1);
      -Ms\Ks, -Ms\Cs];
  
B1s = [zeros(n_p_1,1); Ms\Tgs]; B2s = [zeros(n_p_1,1); Ms\Tus];

Tau_s = eye(n_p_1) + [zeros(1,n_p_1); -eye(n_p_1-1), zeros(n_p_1-1,1)];

Tau = [Tau_s, zeros(n_p_1);zeros(n_p_1), Tau_s];

A = Tau*As/Tau; B1 = Tau*B1s; B2 = Tau*B2s; 

C = [1, 0, 0, 0, 0, 0, 0, 0;
     0, 1, 0, 0, 0, 0, 0, 0;
     0, 0, 1, 0, 0, 0, 0, 0];

%% Define decision variables and optimization problem:

gamma = sdpvar(1,1);
Phi = sdpvar(n_p_1*2,n_p_1*2);
Theta = sdpvar(n_p_1*2,n_p_1-1,'full');

m_11 = (A')*Phi + Phi*A + (Theta*C) + (Theta*C)';

m_12 = Phi*B1;
m_13 = Theta;
m_14 = eye(2*n_p_1);

m_22 = -gamma*eye(1); m_234 = zeros(1,n_p_1*2 + n_p_1-1);
m_33 = -gamma*eye(n_p_1-1); m_34 = zeros(n_p_1-1,2*n_p_1);
m_44 = -gamma*eye(n_p_1*2);

pos_constr = Phi >= 0;
gain_constr = [m_11, m_12, m_13, m_14; 
               m_12', m_22, m_234;
               [m_13'; m_14'], m_234', [m_33; m_34'], [m_34; m_44]] <= 0;
           
           
tot_constr = pos_constr + gain_constr;

%% Solve the problem
ops = sdpsettings('verbose',verbosity_level,'debug',debug_level, ...
                                             'solver',use_solver{1});

optim = optimize(tot_constr, gamma, ops);
if optim.problem ~= 0 && optim.problem ~= 4
    if optim.problem == 1
        fprintf('Infeasible Problem\n\n');
    else
        error(['Error solving problem: ' optim.info ])
    end
else
    fprintf('FDI Filter gains designed:');
    K = value(Phi)\value(Theta)
end
