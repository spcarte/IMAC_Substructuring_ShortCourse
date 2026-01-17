import numpy as np
from scipy.linalg import eigh, null_space

#%% Setting Up System A

# Defining Mass and Stiffness Parameters
m1  = 10
m2  = 10
k1 = 100000
k2 = 100000

# Assembling Mass and Stiffness Matrices
M_A = np.array([[m1, 0],
                [0, m2]])

K_A = np.array([[k1+k2, -k2],
                [-k2, k2]])

# Solving For Modal Parameters
wn_A, phi_A = eigh(K_A, M_A)
fn_A = np.sqrt(wn_A)/2/np.pi

#%% Setting Up System B

# Defining Mass and Stiffness Parameters
m3  = 8
m4  = 8
k3 = 120000
k4 = 120000

# Assembling Mass and Stiffness Matrices
M_B = np.array([[m3, 0],
                [0, m4]])

K_B = np.array([[k3, -k3],
                [-k3, k3+k4]])

wn_B, phi_B = eigh(K_B, M_B)
fn_B = np.sqrt(wn_B)/2/np.pi

#%% Setting Up System Truth

# Assembling Truth Matrices
M_Truth = np.array([[m1, 0,     0],
                    [0,  m2+m3, 0],
                    [0,  0,     m4]])

K_Truth = np.array([[ k1+k2, -k2,     0],
                    [-k2,     k2+k3, -k3],
                    [ 0,     -k3,     k3+k4]])

# Solving for Truth Modal Parameters
wn_Truth, phi_Truth = eigh(K_Truth, M_Truth)
fn_Truth = np.sqrt(wn_Truth)/2/np.pi

#%% Physical Substructuring

# Assembling uncoupled equations in block diagonal form
M_blk_Phys = np.vstack((np.column_stack((M_A, np.zeros((2,2)))),
                        np.column_stack((np.zeros((2,2)), M_B))))

K_blk_Phys = np.vstack((np.column_stack((K_A, np.zeros((2,2)))),
                        np.column_stack((np.zeros((2,2)), K_B))))

# Defining Constraints
B = np.array([0, 1, -1, 0])[np.newaxis,...]
L= null_space(B)

# Synthesizing equations with Constraints
M_Phys = L.T@M_blk_Phys@L
K_Phys = L.T@K_blk_Phys@L

# Solving for Synthesized Modal Parameters
wn_Phys, phi_Phys = eigh(K_Phys, M_Phys)
fn_Phys = np.sqrt(wn_Phys)/2/np.pi


#%% CMS Modal Substructuring

# Assembling uncoupled equations in a MODAL  block diagonal form
M_blk_Modal = np.vstack((np.column_stack((np.eye(2,2), np.zeros((2,2)))),
                         np.column_stack((np.zeros((2,2)), np.eye(2,2)))))

K_blk_Modal = np.vstack((np.column_stack((np.eye(2,2)*wn_A, np.zeros((2,2)))),
                         np.column_stack((np.zeros((2,2)), np.eye(2,2)*wn_B))))

# Defining physical constraints
B = np.array([0, 1, -1, 0])[np.newaxis,...]

# Casting constraints into the modal domain
Phi_stack = np.vstack((np.column_stack((phi_A, np.zeros((2,2)))),
                       np.column_stack((np.zeros((2,2)), phi_B))))
B_bar =  B@Phi_stack
L = null_space(B_bar)

# Synthesizing uncoupled modal equations
M_Modal = L.T@M_blk_Modal@L
K_Modal = L.T@K_blk_Modal@L

# Solving for synthesized modal parameters
wn_Modal, phi_Phys = eigh(K_Modal, M_Modal)
fn_Modal = np.sqrt(wn_Modal)/2/np.pi

#%% Comparison of Results

print('Truth Natural Frequencies:')
print(fn_Truth)

print('Physical Substructuring Natural Frequencies:')
print(fn_Phys)

print('Modal Substructuring Natural Frequencies:')
print(fn_Modal)