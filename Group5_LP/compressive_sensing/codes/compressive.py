import numpy as np
from scipy.sparse import csr_matrix, find
from scipy.io import loadmat
from scipy.io import savemat

#This algorithm solves the compressive sensing problem where Phi is the invertible, linear DCT.
#We have Phi(x) = f, but then f --> b where we lose the last n-m entries of f. Thus we wish to solve the inverse problem
#A*xhat = b for xhat. We hope that xhat is a good approximation to the sparse vector x. We are trying to minimize the 1-norm of 
#xhat in the hopes that this gives us a sparse solution to the inverse problem to reflect the sparsity of x. We needed to phrase 
#this problem in the stadard form of a constrained LP minimization problem.


def compressive():
    #retrieve that parameters of the problem as created using the matlab script parameters.m that came from cs_ex.m
    data = loadmat('/home/matthew/Documents/education/fall17/apm523/projects/lp_qp/lp/compression_sensing/working_codes/compression_data_n500_m400.mat',squeeze_me=1);
    N = data['n']
    m = data['m']
    A= data['A'];
    b_row =data['b'];
    print('A')
    print(A)
    #------------------------------------------------------------------------
    
    n = 2*N	#this is because of the nature of the compression problem setup A = [A, -A]
#turn b into a column vector
    b = np.zeros((m,1))
    for i in range(m):
    	b[i][0] = b_row[i]
    i = 0
    print('b')
    print(b) 
    eta_k=0.99
    
    B = -1 * A
    A = np.c_[A,B]	
    A_tran = A.transpose()
    e=np.ones((n,1))
    c = np.ones((n,1))	
    
    # setting up initial values
    x_start= 10 * np.ones((n,1))
    lambda_start=np.zeros((m,1))
    s_start=10 * np.ones((n,1))
    
    # Initializing x_k,lambda_k,s_k
    x_k = x_start
    lambda_k = lambda_start
    s_k = s_start
    
    # Algorithm 14.3 (Predictor-Corrector Algorithm)
    for k in range(0,80):  #
        # Set (x,lambda,s)=(x_k, lambda_k, s_k)
        x = x_k
        lam = lambda_k
        s = s_k
    
        # make diagonal matrix
        S = np.diag(s.transpose()[0])
        X = np.diag(x.transpose()[0])
    
        # setting equation 14.7 for solving equation 14.30
        test = A@x
        r_b = A@x - b
        r_c = A_tran@lam + s - c
        # End defining r_b and r_c
    
        # making matrix of equation 14.30
        m_11 = np.zeros((n, n))
        m_12 = A_tran
        m_13 = np.identity((n))
        r_1 = np.c_[m_11,m_12,m_13]
    
        m_21 = A
        m_22 = np.zeros((m,m))
        m_23 = np.zeros((m,n))
        r_2 = np.c_[m_21,m_22,m_23]
    
        m_31 = S
        m_32 = np.zeros((n,m))
        m_33 = X
        r_3 = np.c_[m_31,m_32,m_33]
    
        Square1 = np.r_[r_1,r_2,r_3]
    
    
        # computing 'delta_x_affine', 'delta_lambda_affine', 'delta_s_affine' by using equation 14.30
        r_d = -X * S @ e
        R1 = np.r_[-r_c, -r_b, r_d]
        inverse = np.linalg.inv(Square1)
        delta_total = inverse@R1
    
        delta_x_affine = delta_total[0:n]
        delta_lambda_affine = delta_total[n:n+m]
        delta_s_affine = delta_total[n+m:2*n+m]
    
        #computing 'alpha_affine_primal', 'alpha_affine_dual    ', 'mu_affine'
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
    
        delta_x = delta[0:n]
        delta_lambda = delta[n:n+m]
        delta_s = delta[n+m:2*n+m]
    
        #calculating 'alpha_k_primal' and 'alpha_k_dual' fro    m equation (14.38)
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
    
    
        # updating x_k, lambda_k, s_k
        x_k = x_k + alpha_k_primal * delta_x
        lambda_k = lambda_k + alpha_k_dual * delta_lambda
        s_k = s_k + alpha_k_dual * delta_s
    
        print('k =')
        print(k)	
    
    else:
    #extract the nonzero portions of [x+,x-]
        x_opt = np.zeros((N,1))
        for i in range(N):
            if x_k[i] > 1e-10:
                x_opt[i] = x_k[i]
            elif x_k[i+N] > 1e-10:
                x_opt[i] = -1 * x_k[i+N]
            else:
                x_opt[i] = 0.0
        #set the true objective function to represent min ||x||_1
        c_obj = np.ones((N,1))
        opt_func_value = np.dot(x_opt.transpose(),c_obj)	    
        print('Optimal Solutions  >>')
        print(x_opt)
        print('Optimal objective function value   >>')
        print(opt_func_value)
    
    savemat('x_opt.mat',{'xr':x_opt})	

    return x_opt

if __name__ == '__main__':
    Hongjun_compression()
