%% Set-up

clc
clear all
close all
format bank

%% Setting Up System A

% Defining Mass and Stiffness Parameters
m_a  = 10;
k_a = 1e9;
numA = 10;

% Assembling Mass and Stiffness Matrices

M_A = zeros(numA,numA);
K_A = zeros(numA,numA);


M_A(1,1)= m_a;

for ii = 2:numA
   
    M_A(ii,ii)= m_a;
    K_A(ii, ii) =K_A(ii, ii)+k_a;
    K_A(ii-1, ii-1) =K_A(ii-1, ii-1)+k_a;
    K_A(ii, ii-1) =K_A(ii, ii-1)-k_a;
    K_A(ii-1, ii) =K_A(ii-1, ii)-k_a;
end
    
    

% Solving For Modal Parameters
[phi_A,fn_A] = eig(K_A,M_A);
fn_A = sqrt(diag(fn_A))/2/pi;
[fn_A, idx] = sort(fn_A);
phi_A = phi_A(:, idx);


%% Setting Up System B

% Defining Mass and Stiffness Parameters
m_b  = 8;
k_b = 1.2e9;
numB = 10;

% Assembling Mass and Stiffness Matrices

M_B = zeros(numB,numB);
K_B = zeros(numB,numB);


M_B(1,1)= m_b;

for ii = 2:numB
   
    M_B(ii,ii)= m_b;
    K_B(ii, ii) =K_B(ii, ii)+k_b;
    K_B(ii-1, ii-1) =K_B(ii-1, ii-1)+k_b;
    K_B(ii, ii-1) =K_B(ii, ii-1)-k_b;
    K_B(ii-1, ii) =K_B(ii-1, ii)-k_b;
end
    
    

% Solving For Modal Parameters
[phi_B,fn_B] = eig(K_B,M_B);
fn_B = sqrt(diag(fn_B))/2/pi;
[fn_B, idx] = sort(fn_B);
phi_B = phi_B(:, idx);


%% Setting Up System Truth

% Assembling Mass and Stiffness Matrices

M_Truth = zeros(numB+numA-1,numB+numA-1);
K_Truth = zeros(numB+numA-1,numB+numA-1);


M_Truth(1,1)= m_a;
for ii = 2:numA
    M_Truth(ii,ii)= m_a;
    K_Truth(ii, ii) =K_Truth(ii, ii)+k_a;
    K_Truth(ii-1, ii-1) =K_Truth(ii-1, ii-1)+k_a;
    K_Truth(ii, ii-1) =K_Truth(ii, ii-1)-k_a;
    K_Truth(ii-1, ii) =K_Truth(ii-1, ii)-k_a;
end

M_Truth(numA,numA)=M_Truth(numA,numA)+m_b;

for ii = 1:numB-1
    M_Truth(numA+ii,numA+ii)= m_b;
    K_Truth(numA+ii, numA+ii) =K_Truth(numA+ii, numA+ii)+k_b;
    K_Truth(numA+ii-1, numA+ii-1) =K_Truth(numA+ii-1, numA+ii-1)+k_b;
    K_Truth(numA+ii, numA+ii-1) =K_Truth(numA+ii, numA+ii-1)-k_b;
    K_Truth(numA+ii-1, numA+ii) =K_Truth(numA+ii-1, numA+ii)-k_b;
end

% Solving For Modal Parameters
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
B = [zeros(1,numA) zeros(1,numB)];
B(numA)=1;
B(numA+1)=-1;
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

num_modes_A = 5;
num_modes_B = 5;

% Assembling uncoupled equations in a MODAL 
% block diagonal form
M_blk_Modal = blkdiag(eye((num_modes_A)),...
 eye((num_modes_B)));
K_blk_Modal = blkdiag((diag(fn_A(1:num_modes_A))*2*pi).^2,...
 (diag(fn_B(1:num_modes_B)*2*pi).^2));

% Defining physical constraints
B = [zeros(1,numA) zeros(1,numB)];
B(numA)=1;
B(numA+1)=-1;
% Casting constraints into the modal domain
Phi_stack = blkdiag(phi_A(:,1:num_modes_A),phi_B(:,1:num_modes_B));
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
disp('        Truth      Phys Substr  Modal Substr')
disp((real([fn_Truth(1:9) fn_Phys(1:9)  fn_Modal(1:9)])))

 