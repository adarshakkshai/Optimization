import numpy as np
from scipy.sparse import csr_matrix, find

#This alogorithm solves a particular constrained LP problem using the Mehrotra Predictor-Corrector algorithm (see ch.14, Nocedal,Wright
#Numerical Optimization.
#The standard form is min c'x s.t. Ax = b where slack variables may be needed to turn inequality constraints into equality constraints

# adding positive slack variables x3, x4, x5.
n=5
eta_k=0.99
A=np.array([[-2, 1, 1, 0, 0], [-1, 2, 0, 1, 0], [1, 0, 0, 0, 1]])
b=np.array([[2],[7],[3]])
e=np.ones((n,1))
c=np.array([[-1],[-2],[0],[0],[0]])  # -X1-2X2+0X3+0X4+0X5

# setting up initial values
x_start= 10 * np.ones((n,1))
lambda_start=np.zeros((3,1))
s_start=10 * np.ones((n,1))

# Initializing x_k,lambda_k,s_k
x_k = x_start
lambda_k = lambda_start
s_k = s_start

# Algorithm 14.3 (Predictor-Corrector Algorithm)
for k in range(0,165):  #
    # Set (x,lambda,s)=(x_k, lambda_k, s_k)
    x = x_k
    lam = lambda_k
    s = s_k

    # make diagonal matrix
    S = np.diag(s.transpose()[0])
    X = np.diag(x.transpose()[0])

    # setting equation 14.7 for solving equation 14.30
    r_b = A@x - b
    A_tran=A.transpose()
    r_c = A_tran@lam + s - c
    # End defining r_b and r_c

    # making matrix of equation 14.30
    m_11 = np.zeros((5, 5))
    m_12 = A_tran
    m_13 = np.identity((5))
    r_1 = np.c_[m_11,m_12,m_13]

    m_21 = A
    m_22 = np.zeros((3,3))
    m_23 = np.zeros((3,5))
    r_2 = np.c_[m_21,m_22,m_23]

    m_31 = S
    m_32 = np.zeros((5,3))
    m_33 = X
    r_3 = np.c_[m_31,m_32,m_33]

    Square1 = np.r_[r_1,r_2,r_3]


    # computing 'delta_x_affine', 'delta_lambda_affine', 'delta_s_affine' by using equation 14.30
    r_d = -X * S @ e
    R1 = np.r_[-r_c, -r_b, r_d]
    inverse = np.linalg.inv(Square1)
    delta_total = inverse@R1

    delta_x_affine = delta_total[0:5]
    delta_lambda_affine = delta_total[5:8]
    delta_s_affine = delta_total[8:13]

#   # calculate alpha_affine_primal', 'alpha_affine_dual', 'mu_affine' as in equation 14.32 and 14.33
#   for i in range(0, 5):  # i=0;
#        if delta_x_affine[i] < 0:
#            result_vec1 = -x[i] / delta_x_affine[i]
#       else:
#           if delta_x_affine[i] > 0:
#               result_vec1=np.inf
#       if delta_s_affine[i] < 0:
#           result_vec2 = -s[i] / delta_s_affine[i]
#       else:
#           if delta_s_affine[i] > 0:
#               result_vec2=np.inf

#---------------------------------------------------------------------    -------------------------
#Matt's attempt at computing 'alpha_affine_primal', 'alpha_affine_dual    ', 'mu_affine'
    result_vec1 = np.ones((n,1))
    for i in range(n):
        if delta_x_affine[i] < 0:
            result_vec1[i] = -1 * x[i] / delta_x_affine[i]
    alpha_affine_primal = result_vec1.min()
    i = 0
 
    result_vec1 = np.ones((n,1))
    for i in range(n):
        if delta_s_affine[i] < 0:
            result_vec1[i] = -1 * s[i] / delta_s_affine[i]
    alpha_affine_dual = result_vec1.min()
    i = 0

    primal_vec = x + (alpha_affine_primal * delta_x_affine)
    dual_vec = s + (alpha_affine_dual * delta_s_affine)
    mu_affine = primal_vec.transpose().dot(dual_vec) / n
 #End Matt's attempt
 #-----------------------------------------------------------------------------------------------



    # equation 14.32a and 14.32b
#    alpha_affine_primal = min(1,result_vec1.min())
#    alpha_affine_primal = min(1,min(result_vec1))
#    alpha_affine_dual = min(1,result_vec2.min())
#    mu_affine = (x+(alpha_affine_primal*delta_x_affine)).transpose()@(s+(alpha_affine_dual*delta_s_affine))/n

    # setting centering parameter 'sigma' in equation 14.34
    mu = (x.transpose()@s)/n   # equation 14.6
    sigma = (mu_affine/mu)*(mu_affine/mu)*(mu_affine/mu)

    # solve equation 14.35 for delta_x, delta_lambda, delta_s
    # setting equation 14.35
    delta_X_affine = np.diag(delta_x_affine.transpose()[0])
    delta_S_affine = np.diag(delta_s_affine.transpose()[0])

    end_column=-X * S @ e - delta_X_affine * delta_S_affine @ e + sigma @ mu * e
    R2 = np.r_[-r_c, -r_b, end_column]
    delta = inverse@R2

    delta_x = delta[0:5]
    delta_lambda = delta[5:8]
    delta_s = delta[8:13]

    # calculate 'alpha_k_primal' and 'alpha_k_dual' from equation (14.38)
#    for i in range(0, 5):  # i=0
#        if delta_x[i] < 0:
#            result_vec3 = -x[i] / delta_x[i]
#        else:
#            if delta_x[i] > 0:
#                result_vec3 = np.inf
#        if delta_s[i] < 0:
#            result_vec4 = -s[i] / delta_s[i]
#        else:
#            if delta_x_affine[i] > 0:
#                result_vec4=np.inf


#   alpha_k_max_primal = result_vec3.min()
#    alpha_k_max_dual = result_vec4.min()
#    alpha_k_primal = min(1,eta_k*alpha_k_max_primal)
#    alpha_k_dual = min(1,eta_k*alpha_k_max_dual)

#---------------------------------------------------------------------    --------------------------
 #Matt's attempt at calculating 'alpha_k_primal' and 'alpha_k_dual' fro    m equation (14.38)
    result_vec1 = 10000 * np.ones((n,1))
    for i in range(n):
        if delta_x[i][0] < 0:
            result_vec1[i] = -1 * x[i] / delta_x[i][0]
    alpha_max_primal = result_vec1.min()
    i = 0
 
    result_vec1 = 10000 * np.ones((n,1))
    for i in range(n):
        if delta_s[i][0] < 0:
            result_vec1[i] = -1 * s[i] / delta_s[i][0]
    alpha_max_dual = result_vec1.min()
    i = 0
 
    if eta_k * alpha_max_primal < 1:
        alpha_k_primal = eta_k * alpha_max_primal
    else:
        alpha_k_primal = 1
 
    if eta_k * alpha_max_dual < 1:
        alpha_k_dual = eta_k * alpha_max_dual
    else:
        alpha_k_dual = 1
 #End Matt's attempt
 #---------------------------------------------------------------------    --------------------------




    # updating x_k, lambda_k, s_k
    x_k = x_k + alpha_k_primal * delta_x
    lambda_k = lambda_k + alpha_k_dual * delta_lambda
    s_k = s_k + alpha_k_dual * delta_s

    opt_func_value=-x_k[0]-2*x_k[1]

else:
    print('Optimal Solutions  >>')
    print(x)
    print('Optimla objective function value   >>')
    print(opt_func_value)
