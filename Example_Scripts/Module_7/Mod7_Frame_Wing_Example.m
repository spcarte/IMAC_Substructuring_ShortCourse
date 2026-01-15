clear; clc; % Initialize

% % Load System Data % % 
load('FrameWing_Data.mat');

% % Put the systems into a structure array
ss(1)=SS1_FramePlateThinWing;
ss(2)=SS2_PlateThinWing;
ss(3)=SS3_PlateThickWing;
ss(4)=SS4_FramePlateThickWing;

% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % Pick subassembly modes manually,  20 14 12 with pinv
ss(1).Modes = [1:20]; % Frame + Plate + Thin Wing
ss(2).Modes = [1:14]; % Plate + Thin Wing
ss(3).Modes = [1:12]; % Plate + Thick Wing
f_num = 100; % Figure Number to plot to

% % % Pick modes based on Frequency range
% ss(1).Modes = find(1000 > SS1_FramePlateThinWing.wn/2/pi); % Frame + Plate + ThinWing
% ss(2).Modes = find(1000 > SS2_PlateThinWing.wn/2/pi); % Plate + ThinWing
% ss(3).Modes = find(1000 > SS3_PlateThickWing.wn/2/pi); % Plate + ThickWing
% f_num = 200; % Figure Number to plot to
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Reduce systems to the specified modes
ss(1).phi=ss(1).phi(:,ss(1).Modes); ss(1).wn=ss(1).wn(ss(1).Modes); ss(1).zt=ss(1).zt(ss(1).Modes); 
ss(2).phi=ss(2).phi(:,ss(2).Modes); ss(2).wn=ss(2).wn(ss(2).Modes); ss(2).zt=ss(2).zt(ss(2).Modes); 
ss(3).phi=ss(3).phi(:,ss(3).Modes); ss(3).wn=ss(3).wn(ss(3).Modes); ss(3).zt=ss(3).zt(ss(3).Modes); 


% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% By default, constraints are only defined at the DOF on the plate in the
% overlap region of the assembly. This ends up making the problem somewhat
% challenging. However,since the thin and thick wings have DOF that are
% located at physically the same positions, we can expand the constraints
% to also include those. We can do that by just simply replacing the names.
% ss(3).names=ss(2).names; % Uncomment to apply constraints on wings
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



% % % % % % Constraint Equations % % % % % % 

% Determine Overlapping Nodes
[con21,ConnThinWing,ConnFrame] = intersect(ss(2).names,ss(1).names);
[con23,ConnThinWing2,ConnThickWing] = intersect(ss(2).names,ss(3).names);

% Generate Signed Boolean Matrix  
B = zeros(length(con21)+length(con23),length(ss(1).names)+length(ss(2).names)+length(ss(3).names));
for ii = 1:length(con21)
    B(ii, ConnFrame(ii)) = 1;
    B(ii, length(ss(1).names)+ConnThinWing(ii)) = -1;        
end
for ii = 1:length(con23)
    B(ii+length(con21), length(ss(1).names)+ConnThinWing2(ii)) = -1;
    B(ii+length(con21), length(ss(1).names)+length(ss(2).names)+ConnThickWing(ii)) = 1;    
end



% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% % % % % % Constraint Softening Matrix % % % % % % 

% % % MAC at interface - Can help with selecting what modes to include
% % Subsystem 1 to 2 constraints
% MAC_compare(ss(1).phi(ConnFrame,:),ss(1).phi(ConnFrame,:),2,1000); % Frame Thin Self MAC
% MAC_compare(ss(2).phi(ConnThinWing,:),ss(2).phi(ConnThinWing,:),2,1000); % Thin Self MAC
% MAC_compare(ss(1).phi(ConnFrame,:),ss(2).phi(ConnThinWing,:),2,1000); % Frame Thin to Thin Cross MAC
% % Subsystem 2 to 3 constraints
% MAC_compare(ss(2).phi(ConnThinWing2,:),ss(2).phi(ConnThinWing2,:),2,1000); % Thin Self MAC
% MAC_compare(ss(3).phi(ConnThickWing,:),ss(3).phi(ConnThickWing,:),2,1000); % Thick Self MAC
% MAC_compare(ss(2).phi(ConnThinWing2,:),ss(3).phi(ConnThickWing,:),2,1000); % Thin to Thick Cross MAC


% % Using mode shapes of SS2 for the constraint softening
phip = blkdiag(pinv(ss(2).phi(ConnThinWing,:)),pinv(ss(2).phi(ConnThinWing2,:)));


% % SVD of SS1 & SS2 for decoupling, and SS2 & SS3 for coupling
% [a11,a12,a13]=svd([ss(1).phi(ConnFrame,:) ss(2).phi(ConnThinWing,:)]);
% svd_val1 = (diag(a12)./max(a12(:)))*100; inds1=svd_val1>30; %10 - percentage of values
% [a21,a22,a23]=svd([ss(2).phi(ConnThinWing2,:) ss(3).phi(ConnThickWing,:)]);
% svd_val2 = (diag(a22)./max(a22(:)))*100; inds2=svd_val2>30; %18
% phip = blkdiag(a11(:,inds1).' , a21(:,inds2).');
% 
% % Looking at the reconstructed mode shapes from the truncated singular vectors
% % SS1 and SS2
% shps_recon1 = a11(:,inds1)*a12(inds1,inds1)*a13(:,inds1)';
% plot(diag(MAC_compare([ss(1).phi(ConnFrame,:) ss(2).phi(ConnThinWing,:)],shps_recon1)),'.'); xline(size(ss(1).phi,2)+.5);
% MAC_compare([ss(1).phi(ConnFrame,:) ss(2).phi(ConnThinWing,:)],shps_recon1,3,1000);
% MAC_compare(ss(1).phi(ConnFrame,:),ss(1).phi(ConnFrame,:),3,1000);
% MAC_compare(ss(1).phi(ConnFrame,:),shps_recon1(:,1:size(ss(1).phi,2)),3,1000);
% MAC_compare(ss(2).phi(ConnThinWing,:),ss(2).phi(ConnThinWing,:),3,100);
% MAC_compare(ss(2).phi(ConnThinWing,:),shps_recon1(:,size(ss(1).phi,2)+1:end),3,1000);
% 
% % SS2 and SS3
% shps_recon2 = a21(:,inds2)*a22(inds2,inds2)*a23(:,inds2)';
% MAC_compare([ss(2).phi(ConnThinWing2,:) ss(3).phi(ConnThickWing,:)],shps_recon2,3,1000);
% MAC_compare(ss(2).phi(ConnThinWing2,:),ss(2).phi(ConnThinWing2,:),3,1000);
% MAC_compare(ss(2).phi(ConnThinWing2,:),shps_recon2(:,1:size(ss(2).phi,2)),3,1000);
% MAC_compare(ss(3).phi(ConnThickWing,:),ss(3).phi(ConnThickWing,:),3,1000);
% MAC_compare(ss(3).phi(ConnThickWing,:),shps_recon2(:,size(ss(2).phi,2)+1:end),3,1000);
% 
% % ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


% Shouldnt have to change anything below here

% % % % % % Primal CMS to Assemble the Subsystems % % % % % % 

% Modal Equation of Motion Terms for all Subsystems
mM = blkdiag(eye(length(ss(1).wn)),-eye(length(ss(2).wn)),eye(length(ss(3).wn))); % Mass
mC = blkdiag(diag(2*ss(1).wn.*ss(1).zt), -diag(2*ss(2).wn.*ss(2).zt), diag(2*ss(3).wn.*ss(3).zt)); % Damping
mK = blkdiag(diag(ss(1).wn.^2), -diag(ss(2).wn.^2), diag(ss(3).wn.^2)); % Stiffness
Phi_stack = blkdiag(ss(1).phi,ss(2).phi,ss(3).phi); % Compile All Mode Shapes

B_bar = B*Phi_stack; % Modal Constraints
B_bar2 = phip*B_bar; % Softened Constraints
L = null(B_bar2); % Localization Matrix

% Transform systems to coupled coordinates
Meta = L.'*mM*L;
Ceta = L.'*mC*L;
Keta = L.'*mK*L;

% Modal Parameters of Assembled System
[Psi,lam] = eig(Keta,Meta,'vector');
[~,ind]=sort(abs(lam)); lam=lam(ind); Psi=Psi(:,ind);
Psi = Psi./sqrt(diag(Psi.'*Meta*Psi))'; 
PHI=Phi_stack*L*Psi;
fn_est = sqrt(lam)/2/pi;
zt_est = diag(Psi.'*Ceta*Psi)/2./(fn_est*2*pi)*100;
% PHI_est = PHI([1:75 286:339],:); % Frame and plate DOFs, and Thick Wing SS DOFs
% PHI_est = PHI([1:24 130:180 286:339],:); % Frame DOFs, thin plate DOFs, and Thick Wing SS DOFs
PHI_est = PHI([1:24 235:339],:); % Frame DOFs and Plate+Thick Wing SS DOFs



% % % % % % Results % % % % % % 
MAC_compare(ss(4).phi(:,ss(4).wn/2/pi<max(abs(fn_est))),PHI_est,2,f_num,ss(4).wn( ...
    ss(4).wn/2/pi<max(abs(fn_est)))/2/pi,fn_est); xlabel('CMS Result'); ylabel('Truth'); 


% % % Find what Predicted mode best matches each Truth Mode
MAC_matrix = MAC_compare(ss(4).phi,PHI_est); % Based on MAC values
[val,inds]=max(MAC_matrix,[],2);
% [Truth wn, Predicted wn, wn Error % , Truth zt, Predicted zt, zt Error % , MAC Value]
SS_Pred = [ss(4).wn/2/pi fn_est(inds) (fn_est(inds)-ss(4).wn/2/pi)./(ss(4).wn/2/pi)*100 ...
           ss(4).zt*100  zt_est(inds) (zt_est(inds)-ss(4).zt*100)./(ss(4).zt*100)*100 val];
SS_Pred(val>.5,:); % Show the ones that have a MAC Value > 0.5


% % % FRF Plot - Mean of all DOF
fs = linspace(1,1000,2000)'; ws=fs*2*pi; % Frequency Range

% Initial Data - SS1
zt_ss1=SS1_FramePlateThinWing.zt(SS1_FramePlateThinWing.wn>5); 
fn_ss1=SS1_FramePlateThinWing.wn(SS1_FramePlateThinWing.wn>5)/2/pi; 
phi_ss1=SS1_FramePlateThinWing.phi(:,SS1_FramePlateThinWing.wn>5);
Hm_ss1 = 1./(-ws.^2 + 2i.*ws.*zt_ss1'.*fn_ss1'*2*pi + (fn_ss1*2*pi)'.^2); % Modal FRFs
H_ss1 = pagemtimes(phi_ss1,phi_ss1'.*permute(Hm_ss1,[2 3 1]));  % Physical FRF matrix

% Truth Assembly Data - SS4
zt_ss4=SS4_FramePlateThickWing.zt; 
fn_ss4=SS4_FramePlateThickWing.wn/2/pi; 
phi_ss4=SS4_FramePlateThickWing.phi;
Hm_ss4 = 1./(-ws.^2 + 2i.*ws.*zt_ss4'.*fn_ss4'*2*pi + (fn_ss4*2*pi)'.^2); % Modal FRFs
H_ss4 = pagemtimes(phi_ss4,phi_ss4'.*permute(Hm_ss4,[2 3 1]));  % Physical FRF matrix

% CMS Results - SS1-SS2+SS3
zt_CMS=zt_est(abs(fn_est)>5)/100; fn_CMS=fn_est(abs(fn_est)>5); phi_CMS=PHI_est(:,abs(fn_est)>5);
Hm_CMS = 1./(-ws.^2 + 2i.*ws.*zt_CMS'.*fn_CMS'*2*pi + (fn_CMS*2*pi)'.^2); % Modal FRFs
H_CMS = pagemtimes(phi_CMS,phi_CMS'.*permute(Hm_CMS,[2 3 1]));  % Physical FRF matrix

figure(f_num+1); clf;
semilogy(fs, squeeze(mean(abs(H_ss1),[1 2])),'linewidth',1); hold on
semilogy(fs, squeeze(mean(abs(H_ss4),[1 2])),'linewidth',3); 
semilogy(fs, squeeze(mean(abs(H_CMS),[1 2])),'linewidth',3); hold off;
xlabel('Frequency (Hz)'); ylabel('FRF Magnitude'); 
legend('Initial: SS1','Truth: SS4','CMS Result','location','northeast')
title('FRFs'); grid on; axis tight;



% % Drivepoint FRFs - Pick a specific DOF
% fs = linspace(1,1000,4000)'; ws=fs*2*pi;
% 
% % DOF1=5; DOF2=15; % 101Y+: Top of Frame Y
% % DOF1=24; DOF2=DOF1; % 108Z+: Bottom of Frame Z
% DOF1=78; DOF2=DOF1; % 400Z+: Left Bottom Wing Tip Z
% 
% % FRF Synthesis
% zt_CMS=zt_est(7:end)/100; fn_CMS=fn_est(7:end); phi_CMS=PHI_est(:,7:end);
% Hm_CMS = (1i.*ws)./(-ws.^2 + 2i.*ws.*zt_CMS'.*fn_CMS'*2*pi + (fn_CMS*2*pi)'.^2); % Modal FRFs
% H_CMS = pagemtimes(phi_CMS(DOF1,:),phi_CMS(DOF2,:)'.*permute(Hm_CMS,[2 3 1]));  % Physical FRF matrix
% 
% % FRF Synthesis
% zt_ss1=ss(1).zt(7:end); fn_ss1=ss(1).wn(7:end)/2/pi; phi_ss1=ss(1).phi(:,7:end);
% Hm_ss1 = (1i.*ws)./(-ws.^2 + 2i.*ws.*zt_ss1'.*fn_ss1'*2*pi + (fn_ss1*2*pi)'.^2); % Modal FRFs
% H_ss1 = pagemtimes(phi_ss1(DOF1,:),phi_ss1(DOF2,:)'.*permute(Hm_ss1,[2 3 1]));  % Physical FRF matrix
% 
% % FRF Synthesis
% zt_ss4=ss(4).zt; fn_ss4=ss(4).wn/2/pi; phi_ss4=ss(4).phi;
% Hm_ss4 = (1i.*ws)./(-ws.^2 + 2i.*ws.*zt_ss4'.*fn_ss4'*2*pi + (fn_ss4*2*pi)'.^2); % Modal FRFs
% H_ss4 = pagemtimes(phi_ss4(DOF1,:),phi_ss4(DOF2,:)'.*permute(Hm_ss4,[2 3 1]));  % Physical FRF matrix
% 
% figure(4002); clf;
% semilogy(fs, abs(squeeze(abs(H_ss1))),'linewidth',1); hold on
% semilogy(fs, abs(squeeze(abs(H_ss4))),'linewidth',3); 
% semilogy(fs, abs(squeeze(abs(H_CMS))),'linewidth',3); hold off;
% xlabel('Frequency (Hz)'); ylabel('FRF Magnitude'); 
% legend('Initial','Truth','CMS','location','northeast')
% title('FRFs'); grid on; axis tight;
