% NUMERICAL EXAMPLE CASES FOR IMAC SHORT COUSR ON DYNAMIC SUBSTRUCTURING
%   This script goes through a case of utilizing the Transmission Simulator
%   method to improve the results from using CMS to join two numerical beam
%   models. This is extrapolated from the accompanying Numerical_Beam_Example.m 
%   script in the same folder. 
% 
% - Dimensions: height, width, total length, length ratio where the beam is
% split into the A and B subsystem, and total number of elements to use in
% modeling the full length beam. 
% 
% - Material Properties: Elastic Modulus, Density, Poissons Ratio, and 
% Rayleigh Damping, or Proportional Damping, parameters (gamma & beta),
% that are used to generate a C matrix via C=gamma*M+beta*K. The gamma and
% beta values are set by specifying the damption ratio at two frequencies.
% 
% - Boundary Conditions: BCr is a ratio of the diagonal elements of the
% stiffness matrices that is added to them. When BCr is small, this 
% effectively adds a soft boundary condition to the beam that, if it is
% free, raises the natural frequency of the rigid body modes to be not
% close to zero. This can help to ensure consistent substructuring results.
% The BC_fix_L & BC_fix_R flags are for setting the left side of Beam A
% and/or the right side of Beam B to be fixed and cantilevers. 
% 
% - Frequency Vectors: Set the frequency range and number of frequency
% lines for use in generating FRFs. This also sets the range over which 
% the FBS results are computed. 
% 
% - DVA Flag: Term to set if the FRFs are evaluated as displacement,
% velocity, or acceleration. 
% 

clc; clear; % Initialize Workspace

% % % % PARAMETERS TO DEFINE BEAM SUBSYSTEMS % % % % 

% % Geometry of Beams
Beam_h = 0.02; % height (m)
Beam_b = 0.02; % width  (m)
Beam_L = 1.00; % length (m) 
Beam_Ratio = .65; % Ratio of length where beams split
Beam_Elems = 100; % Total number of elements in full beam

% % % Steel Material Properties
% E = 200e9; % N/m^2
% rho = 7800; % kg/m^3
% nu = 0.3; % -

% % Aluminum Material Properties
E = 70e9; % N/m^2
rho = 2700; % kg/m^3
nu = 0.33; % -

% % Proportional Damping parameters
w1 =   10*2*pi; z1=.01; % 1st Specified wn & zeta, 10Hz & 1%
w2 = 1000*2*pi; z2=.01; % 2nd Specified wn & zeta, 1000Hz & 1%
gamma = 2*w1*w2*(z1*w2-z2*w1)/(w2^2-w1^2); % Mass Proportional Term
beta = 2*(z2*w2-z1*w1)/(w2^2-w1^2); % Stiffness Proportional Term

% % Boundary Condition
BCr = 1e-9; % Ratio of stiffness to add along K diagonal to have a soft BC
BC_fix_L = true; % Left side of Beam A free or fixed
BC_fix_R = false; % Right side of Beam B free or fixed

% % Frequency Vector for FRFs
fs = linspace(1,1000,1000)'; ws=fs*2*pi;

% Term to specify if FRFs are Displacement, Velocity or Acceleration
% H_DVA = 1; % Displacement
H_DVA = 1i*ws; % Velocity
% H_DVA = -ws.^2; % Acceleration

%% CONSTRUCT THE MODELS OF BEAM A, BEAM B, AND THE TRUTH ASSEMBLY BEAM C

% % % % Create a Model of Beam A % % % % 
N_A = round(Beam_Elems*Beam_Ratio)+1; % number of Nodes in Beam A
[K_A,M_A,node_A] = beamkm_ur(N_A,(N_A-1)*Beam_L/Beam_Elems,Beam_b,Beam_h,E,rho,nu,'2D');
DOF_A = reshape((1:size(K_A,1)/2)+[.1;.5],[],1); % DOF labels

% Boundary Conditions
K_A = K_A+diag(diag(K_A))*BCr; % Soft Distributed BC 
if BC_fix_L % Fix Left End
M_A=M_A(3:end,3:end); K_A=K_A(3:end,3:end); 
DOF_A=DOF_A(3:end); N_A=N_A-1; node_A=node_A(2:end);
end

% Generate Modal Model
[phi_A,fn_A] = eig(K_A,M_A,'vector'); % Generalized Eigen Solution
[fn_A, idx]=sort(fn_A); phi_A=phi_A(:,idx); % Make sure they are in Order
phi_A = phi_A./sqrt(diag(phi_A'*M_A*phi_A))'; % Mass normalize mode shapes
fn_A = sqrt(abs(fn_A))/2/pi; % Natural Frequencies

% Proportional Damping
C_A = beta.*K_A + gamma.*M_A; % Damping Matrix
zt_A = diag(phi_A'*C_A*phi_A)./(4*pi*fn_A); % Modal Damping Ratios

% FRF Synthesis
Hm_A = H_DVA./(-ws.^2 + 2i.*ws.*zt_A'.*fn_A'*2*pi + (fn_A*2*pi)'.^2); % Modal FRFs
H_A = pagemtimes(phi_A,phi_A'.*permute(Hm_A,[2 3 1]));  % Physical FRF matrix




% % % % Create a Model of Beam B % % % % 
N_B = Beam_Elems+1-round(Beam_Elems*Beam_Ratio); % number of Nodes in Beam B
[K_B,M_B,node_B] = beamkm_ur(N_B,(N_B-1)*Beam_L/Beam_Elems,Beam_b,Beam_h,E,rho,nu,'2D');
node_B = node_B+node_A(end); % Node Positions
DOF_B = reshape((1:size(K_B,1)/2)+[.1;.5],[],1)-1+floor(DOF_A(end-1)); % DOF labels

% Boundary Conditions
K_B = K_B+diag(diag(K_B))*BCr; % Soft Distributed BC 
if BC_fix_R % Fix Right End
M_B=M_B(1:end-2,1:end-2); K_B=K_B(1:end-2,1:end-2); 
DOF_B=DOF_B(1:end-2); N_B=N_B-1; node_B=node_B(1:end-1); 
end

% Generate Modal Model
[phi_B,fn_B] = eig(K_B,M_B,'vector'); % Generalized Eigen Solution
[fn_B, idx]=sort(fn_B); phi_B=phi_B(:,idx); % Make sure they are in Order
phi_B=phi_B./sqrt(diag(phi_B'*M_B*phi_B))'; % Mass normalize mode shapes
fn_B = sqrt(abs(fn_B))/2/pi; % Natural Frequencies

% Proportional Damping
C_B = beta.*K_B + gamma.*M_B; % Damping Matrix
zt_B = diag(phi_B'*C_B*phi_B)./(4*pi*fn_B); % Modal Damping Ratios

% FRF Synthesis
Hm_B = H_DVA./(-ws.^2 + 2i.*ws.*zt_B'.*fn_B'*2*pi + (fn_B*2*pi)'.^2); % Modal FRFs
H_B = pagemtimes(phi_B,phi_B'.*permute(Hm_B,[2 3 1]));  % Physical FRF matrix



% % % % Create Truth Model -> Beam C % % % % 
% Superimpose the MCK System Matrices of Beams A & B
node_C = [node_A ; node_B(2:end)]; % Node Locations
DOF_C = [DOF_A ; DOF_B(3:end)]; % DOF Labels

% Construct Global Mass Matrix
M_C = zeros(N_A*2+N_B*2-2);
M_C(1:N_A*2,1:N_A*2) = M_A;
M_C((N_A*2-1):end,(N_A*2-1):end) = M_C((N_A*2-1):end,(N_A*2-1):end) + M_B;

% Construct Global Damping Matrix
C_C = zeros(N_A*2+N_B*2-2);
C_C(1:N_A*2,1:N_A*2) = C_A;
C_C((N_A*2-1):end,(N_A*2-1):end) = C_C((N_A*2-1):end,(N_A*2-1):end) + C_B;

% Construct Global Stiffness Matrix
K_C = zeros(N_A*2+N_B*2-2);
K_C(1:N_A*2,1:N_A*2) = K_A;
K_C((N_A*2-1):end,(N_A*2-1):end) = K_C((N_A*2-1):end,(N_A*2-1):end) + K_B;

% Generate Modal Model
[phi_C,fn_C] = eig(K_C,M_C,'vector'); % Generalized Eigen Solution
phi_C = phi_C./sqrt(diag(phi_C'*M_C*phi_C))'; % Mass normalize mode shapes
fn_C = sqrt(abs(fn_C))/2/pi; % Natural Frequencies
zt_C = diag(phi_C'*C_C*phi_C)./(4*pi*fn_C); % Damping Ratios

% FRF Synthesis
Hm_C = H_DVA./(-ws.^2 + 2i.*ws.*zt_C'.*fn_C'*2*pi + (fn_C*2*pi)'.^2); % Modal FRFs
H_C = pagemtimes(phi_C,phi_C'.*permute(Hm_C,[2 3 1]));  % Physical FRF matrix

%% Component Mode Synthesis Substructuring

% % Keep All Modes
% num_A = true(length(fn_A),1);
% num_B = true(length(fn_B),1);

% Modal Truncation -> Try using modes under 1000Hz
num_A = abs(fn_A)<1000;
num_B = abs(fn_B)<1000;

M_Global_CMS = blkdiag(eye(nnz(num_A)), eye(nnz(num_B)));
C_Global_CMS = blkdiag(diag(2*zt_A(num_A).*fn_A(num_A)*2*pi), diag(2*zt_B(num_B).*fn_B(num_B)*2*pi));
K_Global_CMS = blkdiag(diag((fn_A(num_A)*2*pi).^2), diag((fn_B(num_B)*2*pi).^2));
Phi_Global_CMS = blkdiag(phi_A(:,num_A),phi_B(:,num_B));

% Signed Boolean Matrix
%   We have two constraint equations that join the translation and rotation
%   DOFs at the right end of Beam A and the left end of Beam B
B = zeros(2,size(K_A,1)+size(K_B,1)); 
B(1,N_A*2-1) =  1;
B(1,N_A*2+1) = -1;
B(2,N_A*2  ) =  1;
B(2,N_A*2+2) = -1;

B_CMS = B*Phi_Global_CMS;
L_CMS = null(B_CMS);

M_CMS = L_CMS.'*M_Global_CMS*L_CMS;
C_CMS = L_CMS.'*C_Global_CMS*L_CMS;
K_CMS = L_CMS.'*K_Global_CMS*L_CMS;

[phi_CMS,fn_CMS] = eig(K_CMS,M_CMS,'vector'); % Generalized Eigen Solution
[fn_CMS,idx]=sort(fn_CMS); phi_CMS=phi_CMS(:, idx); % Sort Result
phi_CMS = phi_CMS./sqrt(diag(phi_CMS'*M_CMS*phi_CMS))'; % Mass normalize mode shapes
fn_CMS = sqrt(abs(fn_CMS))/2/pi; % Natural Frequencies
zt_CMS = diag(phi_CMS'*C_CMS*phi_CMS)./(4*pi*fn_CMS); % Damping Ratios
phi_CMS = Phi_Global_CMS*L_CMS*phi_CMS; % Convert back to physical coordinates
phi_CMS = phi_CMS([1:size(K_A,1) (size(K_A,1)+3):end],:); % Remove repeated nodes

% FRF synthesis
Hm_CMS = H_DVA./(-ws.^2 + 2i.*ws.*zt_CMS'.*fn_CMS'*2*pi + (fn_CMS*2*pi)'.^2); % Modal FRFs
H_CMS = pagemtimes(phi_CMS,phi_CMS'.*permute(Hm_CMS,[2 3 1]));  % Physical FRF matrix

% CMS Results
fn_compare_CMS = [fn_C(1:length(fn_CMS)) fn_CMS (fn_CMS-fn_C(1:length(fn_CMS)))./fn_C(1:length(fn_CMS))*100]; 
zt_compare_CMS = [zt_C(1:length(zt_CMS)) zt_CMS (zt_CMS-zt_C(1:length(zt_CMS)))./zt_C(1:length(zt_CMS))*100]; 
MAC_compare(phi_C(:,1:min([size(phi_C,2) size(phi_CMS,2) 20])),...
    phi_CMS(:,1:min([size(phi_C,2) size(phi_CMS,2) 20])),2,4000,fn_C,fn_CMS);
xlabel('Truth'); ylabel('CMS'); 


%% Transmission Simulator Beam
% To improve the CMS results, we will now include a TS in the beam subsystems.
% The TS is overlaid on the right end of A and extends off the left side of B. 

% Define a model of the TS
N_TS = 4; % Number of nodes to use for TS
[K_TS,M_TS,node_TS] = beamkm_ur(N_TS,(N_TS-1)*Beam_L/Beam_Elems,Beam_b*15,Beam_h*15,200e9,7800,.3,'2D');
node_TS = node_TS+node_A(end-N_TS+1);
DOF_TS = reshape((1:size(K_TS,1)/2)+[.1;.5],[],1)-1+floor(DOF_A(end-size(K_TS,1)+1));

K_TS = K_TS+diag(diag(K_TS))*BCr*1e-3; % Soft Distributed BC 

[phi_TS,fn_TS] = eig(K_TS,M_TS,'vector'); % Generalized Eigen Solution
phi_TS = phi_TS./sqrt(diag(phi_TS'*M_TS*phi_TS))'; % Mass normalize mode shapes
fn_TS = sqrt(abs(fn_TS))/2/pi; % Natural Frequencies

% Proportional Damping
C_TS = beta.*K_TS + gamma.*M_TS; 
zt_TS = diag(phi_TS'*C_TS*phi_TS)./(4*pi*fn_TS);

% % FRF synthesis
Hm_TS = H_DVA./(-ws.^2 + 2i.*ws.*zt_TS'.*fn_TS'*2*pi + (fn_TS*2*pi)'.^2); % Modal FRFs
H_TS = pagemtimes(phi_TS,phi_TS'.*permute(Hm_TS,[2 3 1]));  % Physical FRF matrix


%% Create a modal model of A+TS
% Superimpose TS over the end of A

M_A_TS = M_A;
M_A_TS((end-N_TS*2+1):end,(end-N_TS*2+1):end) = ...
    M_A_TS((end-N_TS*2+1):end,(end-N_TS*2+1):end) + M_TS;

C_A_TS = C_A;
C_A_TS((end-N_TS*2+1):end,(end-N_TS*2+1):end) = ...
    C_A_TS((end-N_TS*2+1):end,(end-N_TS*2+1):end) + C_TS;

K_A_TS = K_A;
K_A_TS((end-N_TS*2+1):end,(end-N_TS*2+1):end) = ...
    K_A_TS((end-N_TS*2+1):end,(end-N_TS*2+1):end) + K_TS;

[phi_A_TS,fn_A_TS] = eig(K_A_TS,M_A_TS,'vector');
[fn_A_TS, idx] = sort(fn_A_TS); phi_A_TS = phi_A_TS(:, idx);
phi_A_TS = phi_A_TS./sqrt(diag(phi_A_TS'*M_A_TS*phi_A_TS))';
fn_A_TS = sqrt(abs(fn_A_TS))/2/pi;
zt_A_TS = diag(phi_A_TS'*C_A_TS*phi_A_TS)./(4*pi*fn_A_TS);
DOF_A_TS = DOF_A;
node_A_TS = node_A;

% % FRF synthesis
Hm_A_TS = H_DVA./(-ws.^2 + 2i.*ws.*zt_A_TS'.*fn_A_TS'*2*pi + (fn_A_TS*2*pi)'.^2); % Modal FRFs
H_A_TS = pagemtimes(phi_A_TS,phi_A_TS'.*permute(Hm_A_TS,[2 3 1]));  % Physical FRF matrix


%% Create a modal model of TS+B
% Append TS on to the end of B

M_B_TS = zeros(N_TS*2+N_B*2-2);
M_B_TS(1:N_TS*2,1:N_TS*2) = M_TS;
M_B_TS((N_TS*2-1):end,(N_TS*2-1):end) = ...
    M_B_TS((N_TS*2-1):end,(N_TS*2-1):end) + M_B;

C_B_TS = zeros(N_TS*2+N_B*2-2);
C_B_TS(1:N_TS*2,1:N_TS*2) = C_TS;
C_B_TS((N_TS*2-1):end,(N_TS*2-1):end) = ...
    C_B_TS((N_TS*2-1):end,(N_TS*2-1):end) + C_B;

K_B_TS = zeros(N_TS*2+N_B*2-2);
K_B_TS(1:N_TS*2,1:N_TS*2) = K_TS;
K_B_TS((N_TS*2-1):end,(N_TS*2-1):end) = ...
    K_B_TS((N_TS*2-1):end,(N_TS*2-1):end) + K_B;

[phi_B_TS,fn_B_TS] = eig(K_B_TS,M_B_TS,'vector');
[fn_B_TS, idx] = sort(fn_B_TS); phi_B_TS = phi_B_TS(:, idx);
phi_B_TS=phi_B_TS./sqrt(diag(phi_B_TS'*M_B_TS*phi_B_TS))';
fn_B_TS = sqrt(abs(fn_B_TS))/2/pi;
zt_B_TS = diag(phi_B_TS'*C_B_TS*phi_B_TS)./(4*pi*fn_B_TS);
DOF_B_TS = [DOF_TS(1:end-2) ; DOF_B];
node_B_TS = [node_TS(1:end-1) ; node_B];

% % FRF synthesis
Hm_B_TS = H_DVA./(-ws.^2 + 2i.*ws.*zt_B_TS'.*fn_B_TS'*2*pi + (fn_B_TS*2*pi)'.^2); % Modal FRFs
H_B_TS = pagemtimes(phi_B_TS,phi_B_TS'.*permute(Hm_B_TS,[2 3 1]));  % Physical FRF matrix


%% Transmission Simulator Method for CMS Substructuring

% Use modes under 1000Hz
num_A_TS = abs(fn_A_TS)<1000;
num_B_TS = abs(fn_B_TS)<1000;
num_TS   = abs(fn_TS  )<1000;

% Manually Select Modes
% num_A_TS=1:10; 
% num_B_TS=1:10; 
% num_TS=1:3;

M_Global_TSM = blkdiag(eye(nnz(num_A_TS)), eye(nnz(num_B_TS)), -2*eye(nnz(num_TS)));
C_Global_TSM = blkdiag(diag(2*zt_A_TS(num_A_TS).*fn_A_TS(num_A_TS)*2*pi), ...
    diag(2*zt_B_TS(num_B_TS).*fn_B_TS(num_B_TS)*2*pi), -2*diag(2*zt_TS(num_TS).*fn_TS(num_TS)*2*pi));
K_Global_TSM = blkdiag(diag((fn_A_TS(num_A_TS)*2*pi).^2), ...
    diag((fn_B_TS(num_B_TS)*2*pi).^2), -2*diag((fn_TS(num_TS)*2*pi).^2));
Phi_Global_TSM = blkdiag(phi_A_TS(:,num_A_TS),phi_B_TS(:,num_B_TS),phi_TS(:,num_TS));
Phi_soften = phi_TS(:,num_TS);
DOF_TSM = [DOF_A_TS;DOF_B_TS;DOF_TS];

% Full Constraints at all overlapping DOF between A, B & TS
B2 = zeros(2*N_TS*2,size(K_A,1)-2+size(K_B,1)+size(K_TS,1)+size(K_TS,1));
for ii=1:N_TS*2
    B2(ii,ii+N_A*2-N_TS*2)=1;
    B2(ii,ii+size(B2,2)-N_TS*2)=-1;
    B2(ii+N_TS*2,ii+N_A*2)=1;
    B2(ii+N_TS*2,ii+size(B2,2)-N_TS*2)=-1;
end
% Reduce to only Translation Constraints
B2 = B2(1:2:end,:);
Phi_soften = Phi_soften(1:2:end,:);

B_TSM = B2*Phi_Global_TSM; % CMS Constraints
B_TSM = blkdiag(pinv(Phi_soften),pinv(Phi_soften))*B_TSM; % Softened Modal Constraints
L_TSM = null(B_TSM);

M_TSM = L_TSM.'*M_Global_TSM*L_TSM;
C_TSM = L_TSM.'*C_Global_TSM*L_TSM;
K_TSM = L_TSM.'*K_Global_TSM*L_TSM;

[phi_TSM,fn_TSM] = eig(K_TSM,M_TSM,'vector');
[fn_TSM, idx] = sort(fn_TSM); phi_TSM = phi_TSM(:, idx);
phi_TSM=phi_TSM./sqrt(diag(phi_TSM'*M_TSM*phi_TSM))';
fn_TSM = sqrt(abs(fn_TSM))/2/pi;
zt_TSM = diag(phi_TSM'*C_TSM*phi_TSM)./(4*pi*fn_TSM); % Damping Ratios
phi_TSM = Phi_Global_TSM*L_TSM*phi_TSM; % Back to physical
[~,inds,~]=unique(DOF_TSM,'stable'); % Find unique nodes
phi_TSM = phi_TSM(inds,:); % Remove repeated nodes

% FRF Synthesis
Hm_TSM = H_DVA./(-ws.^2 + 2i.*ws.*zt_TSM'.*fn_TSM'*2*pi + (fn_TSM*2*pi)'.^2); % Modal FRFs
H_TSM = pagemtimes(phi_TSM,phi_TSM'.*permute(Hm_TSM,[2 3 1]));  % Physical FRF matrix

% TSM Results
fn_compare_TSM = [fn_C(1:length(fn_TSM)) fn_TSM (fn_TSM-fn_C(1:length(fn_TSM)))./fn_C(1:length(fn_TSM))*100]; 
zt_compare_TSM = [zt_C(1:length(zt_TSM)) zt_TSM (zt_TSM-zt_C(1:length(zt_TSM)))./zt_C(1:length(zt_TSM))*100]; 
MAC_compare(phi_C(:,1:min([size(phi_C,2) size(phi_TSM,2) 50])),...
    phi_TSM(:,1:min([size(phi_C,2) size(phi_TSM,2) 50])),2,4001,fn_C,fn_TSM);
xlabel('Truth'); ylabel('TSM'); 


% % MAC Plots of interface DOF in TSM
% inds_A_TS=find(sum(abs(B2(:,1:size(phi_A_TS,1))),1)>eps);
% inds_B_TS=find(sum(abs(B2(:,size(phi_A_TS,1)+1:size(phi_A_TS,1)+size(phi_B_TS,1))),1)>eps);
% inds_TS=find(sum(abs(B2(:,size(phi_A_TS,1)+size(phi_B_TS,1)+1:end)),1)>eps);
% MAC_compare(phi_A_TS(inds_A_TS,:),phi_TS(inds_TS,:),3,123);
% MAC_compare(phi_B_TS(inds_B_TS,:),phi_TS(inds_TS,:),3,123);
% MAC_compare(phi_TS(inds_TS,:),phi_TS(inds_TS,:),3,123);


%% Plot Results 

% Drivepoint FRFs
% DOF1=5; DOF2=15; % Manually pick DOF for FRF plotting
% DOF1=N_A*2-1; DOF2=DOF1; % Constraint Drivepoint DOF
DOF1=size(H_C,2)-1; DOF2=DOF1; % Translation Beam Tip DOF
figure(4002); clf;
semilogy(fs, abs(squeeze(abs(H_C(DOF1,DOF2,:)))),'linewidth',3); hold on
semilogy(fs, abs(squeeze(abs(H_CMS(DOF1,DOF2,:)))),':','linewidth',2)
semilogy(fs, abs(squeeze(abs(H_TSM(DOF1,DOF2,:)))),'linewidth',2); hold off;
xlabel('Frequency (Hz)'); ylabel('FRF Magnitude'); 
legend('Truth','CMS','TSM','location','northeast')
title('FRFs at Tip Translation DOF'); grid on; axis tight;


% Setup Mode Shapes for plotting
phi_p_C   = phi_C(1:2:end,:).*sign(phi_C(end-1,:));
phi_p_CMS = phi_CMS(1:2:end,:).*sign(phi_CMS(end-1,:));
phi_p_TSM = phi_TSM(1:2:end,:).*sign(phi_TSM(end-1,:));
N_p = 1:4; % modes shapes to plot
figure(4003); clf;
plot(node_C,phi_p_C(:,N_p),'color',"#0072BD",'linewidth',3); hold on;
plot(node_C,phi_p_CMS(:,N_p),':','color',"#D95319",'linewidth',2);
plot(node_C,phi_p_TSM(:,N_p),'color',"#EDB120",'linewidth',2); hold off;
xlabel('Beam Length (m)'); ylabel('Modal Deflection'); axis tight; grid on; 
title('Mode Shape Results'); ylim([-1.6 2]);
xline(node_A_TS(end),'linewidth',1.5,'color','k');
xline(node_B_TS(1),'linewidth',1.5,'color','k');


% See what mode shapes were used in truncated forms of A and B 
figure(4004); clf;
plot(node_TS,phi_TS(1:2:end,num_TS).*sign(phi_TS(end-1,num_TS)),'color',"#EDB120",'linewidth',3); hold on;
plot(node_A_TS,phi_A_TS(1:2:end,num_A_TS).*sign(phi_A_TS(end-1,num_A_TS)),'color',"#0072BD",'linewidth',2); 
plot(node_B_TS,phi_B_TS(1:2:end,num_B_TS).*sign(phi_B_TS(end-1,num_B_TS)),'color',"#D95319",'linewidth',2); hold off;
xlabel('Beam Length (m)'); ylabel('Modal Deflection'); axis tight;
title('Mode Shapes of Beam A+TS, Beam B+TS and TS'); grid on; 
xline(node_A(end),'linewidth',1.5,'color','k'); ylim([-2.07 2.07]);
% legend('Beam A+TS','','','','Beam B+TS','','','','TS')


% Subsystem FRF with and without the TS
figure(4005); clf;
semilogy(fs, abs(squeeze(abs(H_A(end-1,end-1,:)))),'linewidth',3); hold on
semilogy(fs, abs(squeeze(abs(H_B(1,1,:)))),'linewidth',3); hold off;
xlabel('Frequency (Hz)'); ylabel('FRF Magnitude'); 
legend('Beam A','Beam B'); axis tight; grid on; 
title('FRFs - Standard A Beam and B Beam');

figure(4006); clf;
semilogy(fs, abs(squeeze(abs(H_A_TS(end-1,end-1,:)))),'linewidth',3); hold on
semilogy(fs, abs(squeeze(abs(H_B_TS(1,1,:)))),'linewidth',3);
semilogy(fs, abs(squeeze(abs(H_TS(end-1,end-1,:)))),'linewidth',3); hold off;
xlabel('Frequency (Hz)'); ylabel('FRF Magnitude'); grid on
legend('Beam A+TS','Beam B+TS','TS'); axis tight; 
title('FRFs - A and B Beams with TS');
