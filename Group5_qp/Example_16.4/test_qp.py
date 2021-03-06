import numpy as np
from scipy.sparse import csr_matrix, find

#This alogorithm solves a particular constrained QP problem with linear inequality constraints using the QP version Mehrotra Predictor-Corrector algorithm (see ch.16, Nocedal,Wright
#Numerical Optimization.
#The standard form is min .5*x'Gx + c'x s.t. Ax - y  = b where slack variables y>=0 may be needed to turn inequality constraints into equality constraints

# setting the parameters
n=2
m = 5
G = np.array([[2,0],[0,2]])
A=np.array([[1, -2],[-1,-2],[-1,2],[1,0],[0,1]])
b=np.array([[-2],[-6],[-2],[0],[0]])
e=np.ones((m,1))
c=np.array([[-2],[-5]])  

# setting up initial values
x_start= np.array([[2],[0]]) 
y_start= np.ones((m,1))
lambda_start= np.ones((m,1))

# Initializing x_k,lambda_k,s_k
x_k = x_start
y_k = y_start
lambda_k = lambda_start

# Algorithm 16.4 (Predictor-Corrector Algorithm for QP)
for k in range(0,10):  
    # Set (x,lambda,s)=(x_k, lambda_k, s_k)
    x = x_k
    y = y_k
    lam = lambda_k

    # make diagonal matrix with y_i and lambda-i as diagonal elements, i = 1,2,...,m
    Y = np.diag(y.transpose()[0])
    L = np.diag(lam.transpose()[0])

    # setting equation 16.59 for solving equation 16.58
    A_tran=A.transpose()
    r_d = G@x - A_tran@lam + c 
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
    print('Square1')
    print(Square1)

    # computing 'delta_x_affine', 'delta_y_affine', 'delta_lambda_affine' by using equation 16.58 with sigma = 0
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
    ahat_affine = result_vec1.min()
    i = 0
 
    result_vec2 = np.ones((m,1))
    for i in range(n):
        if delta_lambda_affine[i] < 0:
            result_vec2[i] = -1 * lam[i] / delta_lambda_affine[i]
    ahat_affine_lambda = result_vec2.min()
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
    end_column = -1*L*Y@e - delta_L_affine * delta_Y_affine@e		
    tt = sigma * mu
    end_column = end_column + tt*e
    R2 = np.r_[-r_d, -r_p, end_column]
    delta = inverse@R2

    delta_x = delta[0:n]
    print('delta_x')
    print(delta_x)
    delta_y = delta[n:n+m]
    print('delta_y')
    print(delta_y)
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
    for i in range(n):
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
    print('alpha')
    print(alpha)
    print('delta_x_end')
    print(delta_x)
    y_k = y_k + alpha * delta_y
    lambda_k = lambda_k + alpha * delta_lambda
    print('x_k')
    print(x_k)
    print('y_k')
    print(y_k)
    print('lambda_k')
    print(lambda_k)
    opt_func_value= (x_k[0] - 1)*(x_k[0] - 1) + (x_k[1] -2.5) * (x_k[1] - 2.5)
    print(k)
    print(opt_func_value)

else:
    print('Optimal Solution  >>')
    print(x_k)
    print('Optimal objective function value   >>')
    print(opt_func_value)
