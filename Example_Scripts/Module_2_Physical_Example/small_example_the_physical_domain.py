"""
This is a small example on how to derive the matrices using Substructuring
Primal and Dual Assemblies. Substructure A is a three degree-of freedom
system clamped at its left end and substructure B is a two
degree-of-freedom structure clamped at its right end. The systems are
coupled in their free ends.

2026-01-12
Andreas Linderholt, LNU
"""
import numpy as np
#--------------------------------------------------------------------------
#                                                                Initialize
#--------------------------------------------------------------------------

m1 = 1; m2 = 1; m3 = 1; m4 = 1; m5 = 1;
c1 = 1; c2 = 2; c3 = 3; c4 = 4; c5 = 5;
k1 = 10000; k2 = 20000; k3 = 30000; k4 = 40000; k5 = 50000;

#--------------------------------------------------------------------------
#                                                            Substructure A
#--------------------------------------------------------------------------

MA = np.array([[m1, 0,  0],
               [0,  m2, 0],
               [0,  0,  m3]])

CA = np.array([[(c1+c2), -c2,      0],
               [-c2,     (c2+c3), -c3],
               [ 0,       -c3,      c3]])

KA = np.array([[(k1+k2), -k2,      0],
               [-k2,     (k2+k3), -k3],
               [ 0,      -k3,      k3]])

FA = np.array([1, 0, 0])[..., np.newaxis] # np.newaxis makes this a column vector
GA = np.array(['g1', 'g2', 'g3'])[..., np.newaxis] # np.newaxis makes this a column vector
rsubA, csubA = MA.shape

#--------------------------------------------------------------------------
#                                                            Substructure B
#--------------------------------------------------------------------------

MB = np.array([[m4, 0],
               [0,  m5]])

CB = np.array([[ c4, -c4],
               [-c4, (c4+c5)]])

KB = np.array([[ k4, -k4],
               [-k4, (k4+k5)]])

FB = np.array([0, 0])[...,np.newaxis] # np.newaxis makes this a column vector
GB = np.array(['g4', 'g5'])[...,np.newaxis] # np.newaxis makes this a column vector
rsubB, csubB = MB.shape

#--------------------------------------------------------------------------
#                                       Block diagonal matrices and vectors
#--------------------------------------------------------------------------

M = np.vstack((np.column_stack((MA, np.zeros((rsubA,csubB)))),
               np.column_stack((np.zeros((rsubB,csubA)), MB))))

C = np.vstack((np.column_stack((CA, np.zeros((rsubA,csubB)))),
               np.column_stack((np.zeros((rsubB,csubA)), CB))))

K = np.vstack((np.column_stack((KA, np.zeros((rsubA,csubB)))),
               np.column_stack((np.zeros((rsubB,csubA)), KB))))

F = np.vstack((FA, FB))

G = np.vstack((GA, GB))

#--------------------------------------------------------------------------
#                                               The compatibility matrix, B
#--------------------------------------------------------------------------

B = np.array([0, 0, 1, -1, 0])[np.newaxis,...] # np.newaxis makes this a row vector
rB, cB = B.shape

#--------------------------------------------------------------------------
#                                                The localization matrix, L
#--------------------------------------------------------------------------

L = np.array([[1, 0, 0, 0],
              [0, 1, 0, 0],
              [0, 0, 1, 0],
              [0, 0, 1, 0],
              [0, 0, 0, 1]])

#--------------------------------------------------------------------------
#                                               The three field formulation
#--------------------------------------------------------------------------

# Mx_dot_dot+Cx_dot+K=F+G       (1)            Governing equation of motion
# Bx=0                          (2)                           Compatibility
# L'G=0                         (3)                             Equilibrium

#--------------------------------------------------------------------------
#                                                           Primal assembly
#--------------------------------------------------------------------------

# x=Lxglobal=Lxg                (4)
# (4) in (1) and pre-multiply with L' ==>
# L'MLxg_dot_dot+L'CLxg_dot+L'KLxg=L'F+L'G
# using (3) ==> L'MLxg_dot_dot+L'CLxg_dot+L'KLxg=L'F

Mprimal = L.transpose()@M@L
print('Mprimal:')
print(Mprimal)

Cprimal = L.transpose()@C@L
print('Cprimal:')
print(Cprimal)

Kprimal = L.transpose()@K@L
print('Kprimal:')
print(Kprimal)

Fprimal = L.transpose()@F
print('Fprimal:')
print(Fprimal)

#--------------------------------------------------------------------------
#                                                             Dual assembly
#--------------------------------------------------------------------------

# G=-B'Lambda ==> (3) ==> L'(-B'Lambda)=-L'B'Lambda=-(BL)'Lambda==0
# (3) is, by the choice of G, fulfilled and (1)==>
# Mx_dot_dot+Cx_dot+K=F-B'Lambda ==> Mx_dot_dot+Cx_dot+K+B'Lambda=F
# With the displacement vector [x;Lambda] etc.:

Mdual = np.vstack((np.column_stack((M, np.zeros((cB,rB)))),
                   np.column_stack((np.zeros((rB,cB)), np.zeros((rB,rB))))))
print('Mdual:')
print(Mdual)

Cdual = np.vstack((np.column_stack((C, np.zeros((cB,rB)))),
                   np.column_stack((np.zeros((rB,cB)), np.zeros((rB,rB))))))
print('Cdual:')
print(Cdual)

Kdual = np.vstack((np.column_stack((K, B.transpose())),
                   np.column_stack((B, np.zeros((rB,rB))))))
print('Kdual:')
print(Kdual.astype(int)) # Printed as integer to avoid scientific notation

Fdual=np.vstack((F, np.zeros((rB,1))))
print('Fdual:')
print(Fdual)