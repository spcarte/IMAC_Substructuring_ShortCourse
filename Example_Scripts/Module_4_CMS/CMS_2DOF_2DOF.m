%% Set-up

clc
clear all
close all
format bank

%% Setting Up System A

% Defining Mass and Stiffness Parameters
m1  = 10;
m2  = 10;
k1 = 100000;
k2 = 100000;

% Assembling Mass and Stiffness Matrices
M_A = diag([m1, m2]);
K_A = [k1+k2 -k2; -k2 k2];

% Solving For Modal Parameters
[phi_A,fn_A] = eig(K_A,M_A);
fn_A = sqrt(diag(fn_A))/2/pi;
[fn_A, idx] = sort(fn_A);
phi_A = phi_A(:, idx);


%% Setting Up System B

% Defining Mass and Stiffness Parameters
m3  = 8;
m4  = 8;
k3 = 120000;
k4 = 120000;

% Assembling Mass and Stiffness Matrices
M_B = diag([m3, m4]);
K_B = [k3 -k3; -k3 k3+k4];

% Solving For Modal Parameters
[phi_B,fn_B] = eig(K_B,M_B);
fn_B = sqrt(diag(fn_B))/2/pi;
[fn_B, idx] = sort(fn_B);
phi_B = phi_B(:, idx);

%% Setting Up System Truth

% Assembling Truth Matrices
M_Truth = diag([m1, m2+m3, m4]);
K_Truth = [k1+k2, -k2 0;
       -k2 k2+k3 -k3;
       0 -k3 k3+k4];

% Solving for Truth Modal Parameters
[phi_Truth,fn_Truth] = eig(K_Truth,M_Truth);
fn_Truth = sqrt(diag(fn_Truth))/2/pi;
[fn_Truth, idx] = sort(fn_Truth);
phi_Truth = phi_Truth(:, idx);


%% Physical Substructuring

% Assembling uncoupled equations 
% in block diagnoal form
M_blk_Phys = blkdiag(M_A, M_B);
K_blk_Phys = blkdiag(K_A, K_B);

% Defining Constraints
B = [0 1 -1 0];
L= null(B);

% Synthesizing equations with 
% Constraints
M_Phys = L.'*M_blk_Phys*L;
K_Phys = L.'*K_blk_Phys*L;

% Solving for Synthesized Modal Parameters
[phi_Phys,fn_Phys] = eig(K_Phys,M_Phys);
fn_Phys = sqrt(diag(fn_Phys))/2/pi;

[fn_Phys, idx] = sort(fn_Phys);
phi_Phys = phi_Phys(:, idx);


%% CMS Modal Substructuring

% Assembling uncoupled equations in a MODAL 
% block diagonal form
M_blk_Modal = blkdiag(eye(size(M_A,1)),...
 eye(size(M_B,1)));
K_blk_Modal = blkdiag((diag(fn_A)*2*pi).^2,...
 (diag(fn_B)*2*pi).^2);

% Defining physical constraints
B = [0 1 -1 0];
% Casting constraints into the modal domain
Phi_stack = blkdiag(phi_A,phi_B);
B_bar =  B*Phi_stack;
L= null(B_bar);

% Synthesizing uncoupled modal equations
M_Modal = L.'*M_blk_Modal*L;
K_Modal = L.'*K_blk_Modal*L;

% Solving for synthesized modal parameters
[phi_Modal,fn_Modal] = eig(K_Modal,M_Modal);
fn_Modal = sqrt(diag(fn_Modal))/2/pi;
[fn_Modal, idx] = sort(fn_Modal);
phi_Modal = phi_Modal(:, idx);

%% Comparison of Results
disp(' ')
disp('Truth        Phys Substr  Modal Substr')
disp(num2str([fn_Truth fn_Phys  fn_Modal]))

 