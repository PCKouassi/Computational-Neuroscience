%% HW 5 Question 1 Script
clear
clc
% Initialize weight matrix 
M = zeros(3);

% Initialize number of neurons
N = 3;

% Initialize number of synapses
S = N^2-N;

%% Network Learning

% Initialize pattern vectors
v1 = [1;-1;1];
v2 = [-1;-1;1];
v3 = [1;-1;-1];
 
% Calcualte the transpose multiplication
v1_mat_mul = v1*v1';
v2_mat_mul = v2*v2';
v3_mat_mul = v3*v3';

% adding the new new network weight matrices together 
sum_1_2 = v1_mat_mul+v2_mat_mul;
M_non = v1*v1'+v2*v2';

M = M_non/N;
gamma_stat = [0.6294; 0.8116; -0.7460];
gamma = -1 + (1+1)*rand(3,3);
epselon = 10E-10;

v2m = M*v2;
v3m = M*v3;

v_new_1 = M*v1+epselon*gamma_stat;

v_new_2 = M*v2+epselon*gamma_stat;

v_new_3 = M*v3+epselon*gamma_stat;