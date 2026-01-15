% This is a small example on how to derive the matrices using Substructuring
% Primal and Dual Assemblies. Substructure A is a three degree-of freedom
% system clamped at its left end and substructure B is a two
% degree-of-freedom structure clamped at its right end. The systems are
% coupled in their free ends.
%
% 2026-01-12
% Andreas Linderholt, LNU
%--------------------------------------------------------------------------
%                                                                Initialize
%--------------------------------------------------------------------------
%
m1=1; m2=1; m3=1; m4=1; m5=1;
c1=1; c2=2; c3=3; c4=4; c5=5;
k1=10000; k2=20000; k3=30000; k4=40000; k5=50000;
%
%--------------------------------------------------------------------------
%                                                            Substructure A
%--------------------------------------------------------------------------
%
MA=[m1 0 0;0 m2 0;0 0 m3];
CA=[(c1+c2) -c2 0;-c2 (c2+c3) -c3;0 -c3 c3];
KA=[(k1+k2) -k2 0;-k2 (k2+k3) -k3;0 -k3 k3];
%
FA=[1 0 0]';
GA=['g1';'g2';'g3'];
[rsubA,csubA]=size(MA);
%
%--------------------------------------------------------------------------
%                                                            Substructure B
%--------------------------------------------------------------------------
%
MB=[m4 0;0 m5];
CB=[c4 -c4;-c4 (c4+c5)];
KB=[k4 -k4;-k4 (k4+k5)];
%
FB=[0 0]';
GB=['g4';'g5'];
[rsubB,csubB]=size(MB);
%
%--------------------------------------------------------------------------
%                                       Block diagonal matrices and vectors
%--------------------------------------------------------------------------
%
M=[MA zeros(rsubA,csubB);zeros(rsubB,csubA) MB];
C=[CA zeros(rsubA,csubB);zeros(rsubB,csubA) CB];
K=[KA zeros(rsubA,csubB);zeros(rsubB,csubA) KB];
F=[FA;FB];
G=[GA;GB];
%
%--------------------------------------------------------------------------
%                                               The compatibility matrix, B
%--------------------------------------------------------------------------
%
B=[0 0 1 -1 0];
[rB,cB]=size(B);
%
%--------------------------------------------------------------------------
%                                                The localization matrix, L
%--------------------------------------------------------------------------
%
L=[1 0 0 0;0 1 0 0;0 0 1 0;0 0 1 0;0 0 0 1];
%
%--------------------------------------------------------------------------
%                                               The three field formulation
%--------------------------------------------------------------------------
%
% Mx_dot_dot+Cx_dot+K=F+G       (1)            Governing equation of motion
% Bx=0                          (2)                           Compatibility
% L'G=0                         (3)                             Equilibrium
%
%--------------------------------------------------------------------------
%                                                           Primal assembly
%--------------------------------------------------------------------------
%
% x=Lxglobal=Lxg                (4)
% (4) in (1) and pre-multiply with L' ==>
% L'MLxg_dot_dot+L'CLxg_dot+L'KLxg=L'F+L'G
% using (3) ==> L'MLxg_dot_dot+L'CLxg_dot+L'KLxg=L'F
%
Mprimal=L'*M*L , Cprimal=L'*C*L  ,  Kprimal=L'*K*L  , Fprimal=L'*F
%
%--------------------------------------------------------------------------
%                                                             Dual assembly
%--------------------------------------------------------------------------
%
% G=-B'Lambda ==> (3) ==> L'(-B'Lambda)=-L'B'Lambda=-(BL)'Lambda==0
% (3) is, by the choice of G, fulfilled and (1)==>
% Mx_dot_dot+Cx_dot+K=F-B'Lambda ==> Mx_dot_dot+Cx_dot+K+B'Lambda=F
% With the displacement vector [x;Lambda] etc.:
%
Mdual=[M zeros(cB,rB);zeros(rB,cB) zeros(rB,rB)]
Cdual=[C zeros(cB,rB);zeros(rB,cB) zeros(rB,rB)]
Kdual=[K B';B zeros(rB,rB)]
Fdual=[F;zeros(rB,1)]
