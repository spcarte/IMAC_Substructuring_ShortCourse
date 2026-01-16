% NUMERICAL EXAMPLE CASES FOR IMAC SHORT COUSR ON DYNAMIC SUBSTRUCTURING
%   This script goes through several example cases of utilziing different
%   forms of dynamic substructuring to couple two beams together to form a
%   single longer beam. The first section contains various parameters that
%   can be set to change the dynamic characteristics of the beams. The 
%   remaining sections have partially completed code to assemble Beam A
%   and Beam B using Physical Substructuring, Lagrange Multiplier Frequency
%   Based Substructuring, and Component Mode Synthesis. 
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
%% Section 1: PARAMETERS TO DEFINE BEAM SUBSYSTEMS
clc; clear; % Initialize Workspace

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

% % Proportional Damping Parameters
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

%% Section 2: CONSTRUCT THE MODELS OF BEAM A, BEAM B, AND THE TRUTH ASSEMBLY BEAM C
%  This section does not need to be changed. It gets all necessary 
%  information from Section 1 above.

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


%% Section 3: Physical Substructuring via M,C,K matrices of A & B 


% Q1: Form the Signed Boolean Matrix
%   We have two constraint equations that join the translation and rotation
%   DOFs at the right end of Beam A and the left end of Beam B
B = zeros(2,N_A*2+N_B*2); 
B(1,N_A*2-1) =  1;
B(1,N_A*2+1) = -1;
B(2,N_A*2  ) =  1;
B(2,N_A*2+2) = -1;


% Q2: Generate the Localization Matrix as the null space of B
L_PHS = null(B,'rational'); % Note: By default, null uses singular value 
% decomposition to robustly compute a null space. In this special case,
% that process actually introduces some numerical noise/rounding that
% propagates through and shows up as error in the results. Using the 
% 'rational' input computes the null space via Reduced Row Echelon Form 
% operations, which in general is less accurate, but for this its better. 
% Run this to see the difference: L_PHS = null(B);


% Forming the Uncoupled Global EOM M,C,K Matrices
M_Global_PHS = blkdiag(M_A,M_B);
C_Global_PHS = blkdiag(C_A,C_B);
K_Global_PHS = blkdiag(K_A,K_B);


% Q3: Compute the Coupled System Matrices via Primal Assembly
M_PHS = L_PHS.'*M_Global_PHS*L_PHS;
C_PHS = L_PHS.'*C_Global_PHS*L_PHS;
K_PHS = L_PHS.'*K_Global_PHS*L_PHS;


% Eigenvalues and Eigenvectors of Assembled System
[phi_PHS,fn_PHS] = eig(K_PHS,M_PHS,'vector'); % Generalized Eigen Solution
[fn_PHS, idx] = sort(fn_PHS); phi_PHS = phi_PHS(:,idx); % Sort Result
phi_PHS = phi_PHS./sqrt(diag(phi_PHS'*M_PHS*phi_PHS))'; % Mass normalize mode shapes
fn_PHS = sqrt(abs(fn_PHS))/2/pi; % Natural Frequencies
zt_PHS = diag(phi_PHS'*C_PHS*phi_PHS)./(4*pi*fn_PHS); % Damping Ratios
phi_PHS = L_PHS*phi_PHS; % Convert shapes back to physical coordinates
phi_PHS = phi_PHS([1:N_A*2 (N_A*2+3):end],:); % Remove repeated DOF

% FRF synthesis
Hm_PHS = H_DVA./(-ws.^2 + 2i.*ws.*zt_PHS'.*fn_PHS'*2*pi + (fn_PHS*2*pi)'.^2); % Modal FRFs
H_PHS = pagemtimes(phi_PHS,phi_PHS'.*permute(Hm_PHS,[2 3 1]));  % Physical FRF matrix

% Compare the Physical Substructuring result to superimposing the matrices 
% - They should be exactly the same. 
fn_compare_PHS = [fn_C(1:length(fn_PHS)) fn_PHS (fn_PHS-fn_C(1:length(fn_PHS)))./fn_C(1:length(fn_PHS))*100]; 
zt_compare_PHS = [zt_C(1:length(zt_PHS)) zt_PHS (zt_PHS-zt_C(1:length(zt_PHS)))./zt_C(1:length(zt_PHS))*100]; 
phi_compare_PHS= (abs(phi_C)-abs(phi_PHS))./abs(phi_C)*100;
% % Mass Matrix (and C & K) are also identical
% tmp = L_PHS*M_PHS*L_PHS.'; % Convert back to physical space
% tmp = tmp([1:N_A*2 (N_A*2+3):end],[1:N_A*2 (N_A*2+3):end]); % remove dof
% % figure; surf((abs(M_C-tmp)./abs(M_C)*100),'facecolor','interp','edgecolor','none')
% figure; plot((abs(M_C-tmp)./abs(M_C)*100))


% % % % Plots % % % % 
% Drivepoint FRFs compared to truth
% DOF1=5; DOF2=15; % Manually pick DOF
DOF1=N_A*2-1; DOF2=DOF1; % Translation Constraint Drivepoint DOF
figure(1000); clf;
semilogy(fs,abs(squeeze(abs(H_C(DOF1,DOF2,:)))),'color','k','linewidth',3); hold on
semilogy(fs,abs(squeeze(abs(H_PHS(DOF1,DOF2,:)))),'-.','linewidth',3); 
semilogy(fs,abs(squeeze(abs(H_C(DOF1,DOF2,:)-H_PHS(DOF1,DOF2,:)))),'.','markersize',8);
xlabel('Frequency (Hz)'); ylabel('FRF Magnitude'); grid on
legend('Truth','Physical','Error','location','northeast')
title('Beam C FRFs at Translation Constraint DOF');

% Mode Shapes compared to Truth
N_p = 1:3; % modes shapes to plot
figure(1001); clf;
plot(node_C,phi_C(1:2:end,N_p).*sign(phi_C(end-1,N_p)),'k','linewidth',3); hold on;
plot(node_C,phi_PHS(1:2:end,N_p).*sign(phi_PHS(end-1,N_p)),'-.','color','#D95319','linewidth',3); hold off;
xlabel('Beam Length (m)'); ylabel('Modal Deflection');
title('Mode Shape Results'); axis tight; grid on; 


%% Section 4: Frequency Based Substructuring using Subsystem FRFs
% LM-FBS (Lagrange Multiplier Frequency Based Substructuring)

% Generate Noise to add to the FRFs
% H_n_amp = 1e-6; % Noise Amplitude relative to median of FRF 
H_n_amp = 5e-2; % Noise Amplitude relative to median of FRF 
% -> 1e-6 is effectively no noise. It significantly effects results at 1e-1
H_An = H_n_amp*median(abs(H_A),3).*abs(randn(size(H_A))).*exp(1i*(rand(size(H_A))*2*pi));
H_Bn = H_n_amp*median(abs(H_B),3).*abs(randn(size(H_B))).*exp(1i*(rand(size(H_B))*2*pi));
% Plot of Constraint DOF FRFs and the Noise to be added
figure(2000); semilogy(fs,abs([squeeze(H_A(end-1,end-1,:)) ...
    squeeze(H_B(1,1,:)) squeeze(H_An(end-1,end-1,:)) ...
    squeeze(H_Bn(1,1,:)) ]),'linewidth',2); grid on; axis tight; 
xlabel('Frequency (Hz)'); legend('Beam A','Beam B','A Noise','B Noise');
ylabel('Magnitude'); title('Constraint Translation DOF - FRFs and Noise'); 

% figure(998);
% subplot(2,1,1); semilogy(fs,abs([squeeze(H_A(end-1,end-1,:)) ...
%     squeeze(H_B(1,1,:)) squeeze(H_An(end-1,end-1,:)) squeeze(H_Bn(1,1,:))]),'linewidth',2);
% xlabel('Frequency (Hz)'); legend('Beam A','Beam B','A Noise','B Noise'); grid on; axis tight; 
% ylabel('Magnitude'); title('Constraint Translation DOF - FRFs and Noise'); ylim([5e-5 3]); 
% subplot(2,1,2); semilogy(fs,abs([squeeze(H_A(end-1,end-1,:))+squeeze(H_An(end-1,end-1,:)) ...
%     squeeze(H_B(1,1,:))+squeeze(H_Bn(1,1,:))]),'linewidth',2); grid on; axis tight; ylim([5e-5 3]); 
% xlabel('Frequency (Hz)'); ylabel('Magnitude'); legend('Noisy Beam A','Noisy Beam B');


% Q1 - Form the Reference and Response Signed Boolean Matrices
B_ref = B; 
B_res = B; 
% In this setup, the reference and response constraint DOF are the same


% Q2 - Implement the LM-FBS Assembly Equation
H_FBS = zeros(N_A*2+N_B*2, N_A*2+N_B*2, length(ws));
for ii = 1:length(ws)
    H_Global = blkdiag(H_A(:,:,ii)+H_An(:,:,ii),H_B(:,:,ii)+H_Bn(:,:,ii)); % Subsystem FRFs at current frequency value
    H_FBS(:,:,ii) = H_Global-H_Global*B_ref.'*pinv(B_res*H_Global*B_ref.')*B_res*H_Global; % LM-FBS
end
H_FBS = H_FBS([1:N_A*2 (N_A*2+3):end],[1:N_A*2 (N_A*2+3):end],:); % Remove repeated DOF


% Drivepoint FRFs Plot
% DOF1=5; DOF2=15; % Manually pick DOF for FRF plotting
DOF1=N_A*2-1; DOF2=DOF1; % Translation Constraint Drivepoint DOF
figure(2001); clf;
semilogy(fs, abs(squeeze(abs(H_C(DOF1,DOF2,:)))),'color','k','linewidth',3); hold on
semilogy(fs, abs(squeeze(abs(H_FBS(DOF1,DOF2,:)))),'color',"#EDB120",'linewidth',3)
semilogy(fs, abs(squeeze(abs(H_C(DOF1,DOF2,:)-H_FBS(DOF1,DOF2,:)))),'r.','markersize',8); hold off;
xlabel('Frequency (Hz)'); ylabel('FRF Magnitude'); grid on; axis tight; 
legend('Truth','FBS','Error','location','northeast')
title('Beam C FRFs at Translation Constraint DOF');


%% Section 5: Component Mode Synthesis Substructuring

% Keep All Modes
num_A = true(length(fn_A),1);
num_B = true(length(fn_B),1);

% % Modal Truncation -> Try using modes under 1000Hz
% num_A = abs(fn_A)<1000;
% num_B = abs(fn_B)<1000;

% Form Global Modal M,C,K Matrices
M_Global_CMS = blkdiag(eye(nnz(num_A)), eye(nnz(num_B)));
C_Global_CMS = blkdiag(diag(2*zt_A(num_A).*fn_A(num_A)*2*pi), diag(2*zt_B(num_B).*fn_B(num_B)*2*pi));
K_Global_CMS = blkdiag(diag((fn_A(num_A)*2*pi).^2), diag((fn_B(num_B)*2*pi).^2));
Phi_Global_CMS = blkdiag(phi_A(:,num_A),phi_B(:,num_B));


% Q1: Convert Physical Constraints to Modal
B_CMS = B*Phi_Global_CMS;


% Q2: Generate Localization Matrix
L_CMS = null(B_CMS);


% Q3: Apply Primal Constraints to Modal EOM Matrices 
M_CMS = L_CMS.'*M_Global_CMS*L_CMS;
C_CMS = L_CMS.'*C_Global_CMS*L_CMS;
K_CMS = L_CMS.'*K_Global_CMS*L_CMS;


% Compute Modal Parameters of Assembled System 
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


% % % Plots % % %

% Drivepoint FRFs
% DOF1=5; DOF2=7; % Manually pick DOF for FRF plotting
DOF1=size(H_C,2)-1; DOF2=DOF1; % Translation Beam Tip DOF
figure(3000); clf;
semilogy(fs, abs(squeeze(abs(H_C(DOF1,DOF2,:)))),'color','k','linewidth',3); hold on
semilogy(fs, abs(squeeze(abs(H_CMS(DOF1,DOF2,:)))),'color',"#7E2F8E",'linewidth',3); 
xlabel('Frequency (Hz)'); ylabel('FRF Magnitude'); grid on; axis tight; hold off;
legend('Truth','CMS','location','northeast'); title('Beam C FRFs at Tip Translation DOF');


% Setup Mode Shapes for plotting
N_p = 1:3; % modes shapes to plot
phi_p_C   = phi_C(1:2:end,N_p).*sign(phi_C(end-1,N_p));
phi_p_CMS = phi_CMS(1:2:end,N_p).*sign(phi_CMS(end-1,N_p));
figure(3001); clf;
plot(node_C,phi_p_C(:,:),'k','linewidth',3); hold on;
plot(node_C,phi_p_CMS(:,:),'g','linewidth',2); hold off;
xlabel('Beam Length (m)'); ylabel('Modal Deflection'); axis tight; grid on; 
title('Mode Shape Results');


% MAC Plot of mode shapes
MAC_compare(phi_C(:,1:min([size(phi_C,2) size(phi_CMS,2) 20])),...
    phi_CMS(:,1:min([size(phi_C,2) size(phi_CMS,2) 20])),2,3002,fn_C,fn_CMS);
xlabel('Truth'); ylabel('CMS'); 


% See what mode shapes were used in truncated forms of A and B 
figure(3003); clf;
plot(node_A,phi_A(1:2:end,1:3).*sign(phi_A(end-1,1:3)),'color',"#0072BD",'linewidth',2); hold on;
plot(node_B,phi_B(1:2:end,1:3).*sign(phi_B(end-1,1:3)),'color',"#D95319",'linewidth',2); hold off;
xlabel('Beam Length (m)'); ylabel('Modal Deflection'); axis tight;
title('Mode Shapes of Beam A and Beam B'); grid on; 
xline(node_A(end),'linewidth',3,'color','k');
legend('Beam A','','','Beam B')
