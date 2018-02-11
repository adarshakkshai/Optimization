import numpy as np
from scipy.sparse import csr_matrix, find
import math

#This alogorithm solves a csvm problem by formulating it as a constrained QP problem with linear inequality constraints and solving using the QP version Mehrotra Predictor-Corrector algorithm (see ch.16, Nocedal,Wright Numerical Optimization.
#The standard form is min .5*x'Gx + c'x s.t. Ax - y  = b where slack variables y>=0 may be needed to turn inequality constraints into equality constraints

# importing the parameters of vector length and number of data vectors
num_x = math.floor(np.loadtxt('num_vectors.txt'))	#there are 679 data vectors
len_x = math.floor(np.loadtxt('length_data.txt'))	#each data vector has 9 components
n= len_x + num_x + 1	#n = (size of each x + 1 + number of x vectors) = dimension of the solution vector
m = 2*num_x	#number of rows in constraint matrix A
#build the G matrix
G11 = np.identity(len_x)
G12 = np.zeros((len_x,num_x+1))
G_row1 = np.c_[G11,G12]
G_row2 = np.zeros((num_x+1,n))
G = np.r_[G_row1,G_row2]
#build the n-vector c
C = 1
c11 = np.zeros((1,len_x+1))
c12 = C * np.ones((1,num_x))
c = np.c_[c11,c12]
#build the A matrix. The final 1 in each row Ai1 corrsponds to the y_i*(greek y) in the constraint
#x1 through x9 including the final component of each being 1

data_vector_array = np.loadtxt('data_vectors.txt')	#679 x 10 matrix. columns 1-9 are the data vector, column 10 is y_i = +-1
#create the matrix of data vectors each multiplied by its y_1 = +-1
#data_vectors_adjusted = np.zeros((num_x,len_x))
for i in range(num_x):
    for j in range(len_x):
        data_vector_array[i][j] = data_vector_array[i][len_x] * data_vector_array[i][j]
    j = 0
i = 0	

A12 = np.identity(num_x)
A_top = np.c_[data_vector_array,A12]
A21 = np.zeros((num_x,len_x+1))
A22 = np.identity(num_x)
A_bottom = np.c_[A21,A22]
A_top_tran = A_top.transpose()
A_bottom_tran = A_bottom.transpose()
A_tran = np.c_[A_top_tran,A_bottom_tran] #A_tran = A	#A is actually the transpose of A
A = A_tran.transpose()

#build the RHS constraint vector b
b_top = np.ones((num_x,1))
b_bottom = np.zeros((num_x,1))
b = np.r_[b_top,b_bottom]

# setting up initial values
x_start= np.ones((n,1)) 
y_start= np.ones((m,1))
lambda_start= np.ones((m,1))

# Initializing x_k,lambda_k,s_k
x_k = x_start
y_k = y_start
lambda_k = lambda_start

# Algorithm 16.4 (Predictor-Corrector Algorithm for QP)
for k in range(0,20):  
    # Set (x,lambda,s)=(x_k, lambda_k, s_k)
    x = x_k
    y = y_k
    lam = lambda_k

    # make diagonal matrix with y_i and lambda-i as diagonal elements, i = 1,2,...,m
    Y = np.diag(y.transpose()[0])
    L = np.diag(lam.transpose()[0])

    # setting equation 16.59 for solving equation 16.58
    r_d = G@x - A_tran@lam + c.transpose() 
    r_p= A@x - y - b
    # End defining r_d and r_p

    # making matrix of equation 16.58
    m_11 = G
    m_12 = np.zeros((n, m))
    m_13 = -1 * A_tran
    r_1 = np.c_[m_11,m_12,m_13]

    m_21 = A
    m_22 = -1 * np.identity(m)
    m_23 = np.zeros((m,m))
    r_2 = np.c_[m_21,m_22,m_23]

    m_31 = np.zeros((m,n))
    m_32 = L 
    m_33 = Y
    r_3 = np.c_[m_31,m_32,m_33]

    Square1 = np.r_[r_1,r_2,r_3]

    # computing 'delta_x_affine', 'delta_y_affine', 'delta_lambda_affine' by using equation 16.58 with sigma = 0
    e = np.ones((m,1))
    r_yl = -L * Y @ e
    R1 = np.r_[-r_d, -r_p, r_yl]
    inverse = np.linalg.inv(Square1)
    delta_total = inverse@R1

    delta_x_affine = delta_total[0:n] #0 to n-1
    delta_y_affine = delta_total[n:n+m]	#n to n+m-1
    delta_lambda_affine = delta_total[n+m:n+2*m]	#n+m to n+2m-1

#calculate the complementarity measure mu
    mu = y.transpose()@lam	/ m

#calculate ahat_aff
    result_vec1 = np.ones((m,1))
    for i in range(m):
        if delta_y_affine[i] < 0:
            result_vec1[i] = -1 * y[i] / delta_y_affine[i]
    if result_vec1.min() < 1:
        ahat_affine = result_vec1.min()
    else:
        ahat_affine = 1
    i = 0
    result_vec2 = np.ones((m,1))
    for i in range(m):
        if delta_lambda_affine[i] < 0:
            result_vec2[i] = -1 * lam[i] / delta_lambda_affine[i]
    if result_vec2.min() < 1:
        ahat_affine_lambda = result_vec2.min()
    else:
        ahat_affine_lambda = 1
    i = 0
    if ahat_affine_lambda < ahat_affine:
        ahat_affine = ahat_affine_lambda
    dy_vec = y + (ahat_affine * delta_y_affine)
    dlam_vec = lam + (ahat_affine * delta_lambda_affine)
    mu_affine = dy_vec.transpose().dot(dlam_vec) / m

# setting centering parameter 'sigma' in equation 14.34
    sigma = (mu_affine/mu)*(mu_affine/mu)*(mu_affine/mu)
    # solve equation 16.67 for delta_x, delta_lambda, delta_s
    # setting equation 16.67
    delta_Y_affine = np.diag(delta_y_affine.transpose()[0])
    delta_L_affine = np.diag(delta_lambda_affine.transpose()[0])
    end_column = -1*L*Y@e - delta_L_affine * delta_Y_affine@e + sigma*mu*e	
    R2 = np.r_[-r_d, -r_p, end_column]
    delta = inverse@R2

    delta_x = delta[0:n]
    delta_y = delta[n:n+m]
    delta_lambda = delta[n+m:n+2*m]

#-----------------------------------------------------------------------------------
#calculate a_tau_pri and a_tau_dual
    #set value for tau in (0,1). Larger tau ---> larger step length
    tau = .95
    result_vec1 = np.ones((m,1))
    for i in range(m):
        if delta_y[i] < 0:
            result_vec1[i] = -tau * y[i] / delta_y[i]
    a_tau_pri = result_vec1.min()
    i = 0
 
    result_vec2 = np.ones((m,1))
    for i in range(m):
        if delta_lambda[i] < 0:
            result_vec2[i] = -tau * lam[i] / delta_lambda[i]
    a_tau_dual = result_vec2.min()
    i = 0
    #choose the minimum of these two step lengths
    if a_tau_pri < a_tau_dual:
        alpha = a_tau_pri
    else:
        alpha = a_tau_dual

    
#---------------------------------------------------------------------------------------


    # updating x_k, lambda_k, s_k
    x_k = x_k + alpha * delta_x
    y_k = y_k + alpha * delta_y
    lambda_k = lambda_k + alpha * delta_lambda
    opt_func_value= 0
    for i in range(len_x):
        opt_func_value = opt_func_value + x[i]*x[i]
    print(k)
    print(opt_func_value)

else:
    print('Optimal Solution  >>')
    print('w')   
    print(x_k[0:8])
    print('gamma')
    print(x_k[9])
    print('zeta')
    print(x_k[10:])
    print('Optimal objective function value   >>')
    print(opt_func_value)
    np.savetxt('csvm_data_output.txt',x_k, delimiter=',')	
