function [K,M,nodedef] = beamkm_ur(Nelem,L,b,h,E,rho,nu,varargin);
% Simplified FE model of a uniform beam with a uniform rectangular cross section
%
% [K,M,nodedef] = beamkm_ur(Nelem,L,b,h,E,rho,nu);
% 
% Creates a 3-D FE model for a free beam with properties L,b,h,E,rho,nu.
% Note that this m-file uses a highly simplified formula for the torsional
% stiffness constant that may not be correct for non-circular bars.
%
% The global displacement vector u is ordered as follows:
%   u = [u1 (node1) u2 (node1) u3 (node1) th1 (node1) th2 (node1) th3 (node1) ...
%        u1 (node2) u2 (node2) u3 (node2) th1 (node2) th2 (node2) th3 (node2) ...
%        ...
%
% %%%%%% Alternate calling sequence for 2D FE models. %%%%%%
%
% [K,M,x0] = beamkm_ur(Nelem,L,b,h,E,rho,nu,'2D');
%
% This returns a 1-D beam (bending only) and the corresponding vector of
% node locations x0.  The DOF are sorted as: x1, theta1, x2, theta2, etc...
%
% To obtain a cantilever, use the following:
%
% DOF = DOF(3:end);
% M = M(DOF,DOF); K = K(DOF,DOF);
% 
% M.S. Allen, based on Clark Dohrman's beamkm.m program.
%

Iarea = b*h^3/12;
Iarea2 = b^3*h/12;
G = E/(2*(1-nu));
J = Iarea+Iarea2;
Ixx_perunitL = (1/12)*rho*b*h*(b^2+h^2); % Mass-moment per unit length

elemdef = [(1:Nelem-1)',(2:Nelem)'];
nodedef = [zeros(Nelem,2),(0:(Nelem-1))'*L/(Nelem-1)];
dirdef = [ones(Nelem,1),zeros(Nelem,2)];
propdef = zeros(Nelem-1,6);
    propdef(:,1) = b*h*E;        %AE
    propdef(:,2) = J*G;          %JG
    propdef(:,3) = E*Iarea;      %EI1
    propdef(:,4) = E*Iarea2;     %EI2
    propdef(:,5) = rho*b*h;      %mass/unit length
    propdef(:,6) = Ixx_perunitL; %torsional mass moment/unit length

[K,M] = beamkm(elemdef,nodedef,dirdef,propdef);

if nargin > 7
    DOF = ([(1:6:6*Nelem),(5:6:6*Nelem)]);    % pulling out just DOFs 1 & 5 for each grid point
    DOF = sort(DOF);    % now DOF are x1, theta1, x2, theta2, etc...
    K = K(DOF,DOF);     % new K and M matrices for beam, 2dof per node
    M = M(DOF,DOF);
    nodedef = nodedef(:,3);
end

