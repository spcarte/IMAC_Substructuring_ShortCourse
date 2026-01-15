function [K,M]=beamkm(elemdef,nodedef,dirdef,propdef)
%function [K,M]=beamkm(elemdef,nodedef,dirdef,propdef)
% Input:
%  elemdef(i,1)   = element i node 1
%  elemdef(i,2)   = element i node 2
%  nodedef(j,1:3) = x1, x2, x3 coordinates of node j
%  dirdef(i,1:3)  = direction for bending axis 1
%                   (bending axis 2 is obtained by taking cross
%                   product of axial direction and bending axis 1)
%  propdef(i,1)   = AE of element i
%  propdef(i,2)   = JG of element i
%  propdef(i,3)   = EI for bending axis 1 of element i
%  propdef(i,4)   = EI for bending axis 2 of element i
%  propdef(i,5)   = mass per unit length of element i
%  propdef(i,6)   = torsional mass moment of inertia per unit
%                   length of element i
% Output:
%  K              = assembled stiffness matrix
%  M              = assembled mass matrix
%
% Comments: the global displacement vector u is ordered as follows
%
% u = [u1 (node1) u2 (node1) u3 (node1) th1 (node1) th2 (node1) th3 (node1) ...
%      u1 (node2) u2 (node2) u3 (node2) th1 (node2) th2 (node2) th3 (node2) ...
%      ...
%      u1 (nodem) u2 (nodem) u3 (nodem) th1 (nodem) th2 (nodem) th3 (nodem)]
%
% where m is the number of nodes in the model and uj and thj are the
% displacement and rotation, respectively, about global axis j;
nelem=size(elemdef,1);
nnode=size(nodedef,1);
K=zeros(6*nnode,6*nnode);
M=zeros(6*nnode,6*nnode);
PA= [1 0 0 0 0 0 0 0 0 0 0 0; ...
     0 0 0 0 0 0 1 0 0 0 0 0];
PT= [0 0 0 1 0 0 0 0 0 0 0 0; ...
     0 0 0 0 0 0 0 0 0 1 0 0];
PB1=[0 1 0 0 0 0 0 0 0 0 0 0; ...
     0 0 0 0 0 1 0 0 0 0 0 0; ...
     0 0 0 0 0 0 0 1 0 0 0 0; ...
     0 0 0 0 0 0 0 0 0 0 0 1];
PB2=[0 0 1 0  0 0 0 0 0 0  0 0; ...
     0 0 0 0 -1 0 0 0 0 0  0 0; ...
     0 0 0 0  0 0 0 0 1 0  0 0; ...
     0 0 0 0  0 0 0 0 0 0 -1 0];
for i=1:nelem
    nodes=elemdef(i,1:2);
    dx=nodedef(nodes(2),:)-nodedef(nodes(1),:);
    L=norm(dx);
    AE=propdef(i,1);
    JG=propdef(i,2);
    EI1=propdef(i,3);
    EI2=propdef(i,4);
    rho=propdef(i,5);
    rhoT=propdef(i,6);
    KelemA=AE/L*[ 1 -1; ...
                 -1  1];
    KelemT=JG/L*[ 1 -1; ...
                 -1  1];
    KelemB1=EI1/L^3*[ 12  6*L  -12    6*L; ...
                     6*L 4*L^2  -6*L 2*L^2; ...
                     -12 -6*L   12   -6*L; ...
                     6*L 2*L^2  -6*L  4*L^2];
    KelemB2=EI2/L^3*[ 12  6*L  -12    6*L; ...
                     6*L 4*L^2  -6*L 2*L^2; ...
                     -12 -6*L   12   -6*L; ...
                     6*L 2*L^2  -6*L  4*L^2];
    Kelem=PA'*KelemA*PA + ...
          PT'*KelemT*PT + ...
         PB1'*KelemB1*PB1 + ...
         PB2'*KelemB2*PB2;
    %
    MelemA=rho*L*[1/3 1/6; ...
                  1/6 1/3];
    MelemT=rhoT*L*[1/3 1/6; ...
                   1/6 1/3];
    MelemB1=rho*L*[   13/35    11/210*L   9/70     -13/420*L; ...
                    11/210*L  1/105*L^2  13/420*L -1/140*L^2; ...
                       9/70    13/420*L  13/35     -11/210*L; ...
                   -13/420*L -1/140*L^2 -11/210*L  1/105*L^2];
    MelemB2=rho*L*[   13/35    11/210*L   9/70     -13/420*L; ...
                    11/210*L  1/105*L^2  13/420*L -1/140*L^2; ...
                       9/70    13/420*L  13/35     -11/210*L; ...
                   -13/420*L -1/140*L^2 -11/210*L  1/105*L^2];
    Melem=PA'*MelemA*PA + ...
          PT'*MelemT*PT + ...
         PB1'*MelemB1*PB1 + ...
         PB2'*MelemB2*PB2;
    %
    d1=dx/norm(dx);
    d2=dirdef(i,:)/norm(dirdef(i,:));
    if abs(d1*d2')>1e-6
        fprintf(1,'Warning: axial and bending 1 directions are not orthogonal\n');
        fprintf(1,'dot product of two directions is %g\n',d1*d2');
        fprintf(1,'axial direction to be projected of bending 1 direction\n');
        d2=d2-(d1*d2')*d1;
        d2=d2/norm(d2);
    end
    d3(1)=d1(2)*d2(3)-d1(3)*d2(2);
    d3(2)=d1(3)*d2(1)-d1(1)*d2(3);
    d3(3)=d1(1)*d2(2)-d1(2)*d2(1);
    C=[d1;d2;d3];
    A=zeros(12,12);
    A(1:3,1:3)=C;
    A(4:6,4:6)=C;
    A(7:9,7:9)=C;
    A(10:12,10:12)=C;
    i1=6*(nodes(1)-1)+1:6*nodes(1);
    i2=6*(nodes(2)-1)+1:6*nodes(2);
    K([i1 i2],[i1 i2])=K([i1 i2],[i1 i2])+A'*Kelem*A;
    M([i1 i2],[i1 i2])=M([i1 i2],[i1 i2])+A'*Melem*A;
end
