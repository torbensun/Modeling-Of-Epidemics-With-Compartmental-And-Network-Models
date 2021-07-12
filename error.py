import numpy as np


def delta_kronecker(i, j):
    if i == j:
        return 1
    else:
        return 0

def ds(state,n,i,method,alpha):
    # method ist Vektor mit zwei Einträgen 0 ist SIGMA1_1 ist SIGMA_2
    S_i=state[i]
    I_i=state[i+n-1]
    res=0
    #P[i][j] ist P_{[i->j]}
    for j in range(n):
        res=res+P[i][j]*I_eff(j)
    return -method[0] * alpha * S_i * I_i - method[1] * alpha * (N_eff(i)* S_i *I_eff(i)  + S_i * res) / N(i)


def dI(state, n, i,method, alpha, beta):
    I_i=state[i+n-1]
    return ds(state,n,i,method,alpha) - beta * I_i

def dR(state ,n , beta, p):
    I_i = state[i + n - 1]
    return (1-p) * beta * I_i

def dD(state ,n , beta, p):
    I_i = state[i + n - 1]
    return (p) * beta * I_i

def Jacobian_SS(state, n):
    # n ist die Anzahl der Landkreise
    S=state[0:n-1]
    I=state[n:2*n-1]
    R=state[2*n:3*n-1]
    J_SS=np.identity(n)
    for i in range(n):
        for j in range(n):
            J_SS[i][j]=delta_kronecker(i, j)*ds(state,n,i)/S[j]
    return J_SS

def Jacobian_SI(state, n, method, alpha):
    J_SI = np.identity(n)
    S = state[0:n - 1]
    for i in range(n):
        for l in range(n):
            res1 = -alpha * delta_kronecker(i, l)*method[0]* S[i]
            sum1=0
            for j in range(n):
                sum1=sum1+P[i][j]
            res2 = -alpha*method[1]*N_eff(i)/N(i)*S[i] *(delta_kronecker(i, l) + P[l][i]/N(i)-delta_kronecker(i, l)*sum1/N(i))
            sum2=0
            for j in range(n):
                sum3=0
                for k in range(n):
                    sum3=sum3+(P[j][n]+P[l][j])/N(j)
                sum2=sum2+P[i][j]*(delta_kronecker(j, l) - delta_kronecker(j, l)*sum3)
            res3 = -alpha*method[1]*S[i]/N(i)*sum2
            J_SI[i][l]=res3+res2+res1
    return J_SI

def Jacobian_RI(n, beta, p):
    J_RI = np.identity(n)
    for i in range(n):
        for l in range(n):
            J_RI[i][l]=delta_kronecker(i, l) *(1-p)*beta
    return J_RI

def Jacobian_DI(n, beta, p):
    J_DI = np.identity(n)
    for i in range(n):
        for l in range(n):
            J_DI[i][l]=delta_kronecker(i, l) *(p)*beta
    return J_DI


def Jacobian_II(state, n, method, alpha, beta):
    J_II = np.identity(n)
    S = state[0:n - 1]
    for i in range(n):
        for l in range(n):
            sum1=0
            sum2=0
            sum3=0
            for j in range(n):
                sum1=sum1+P[i][j]
                sum2=sum2+P[l][j]/N(j)
                sum3=sum3+P[l][j]/N(l)
            res=method[0]*delta_kronecker(i, l)+method[1]/N(i)*(N_eff(i)*((delta_kronecker(i, l)+P[l][i])/N(i)-sum1*delta_kronecker(i, l)/N(i)) +P[i][l]+sum2-P[i][l]*sum3)
            J_II[i][l]=-beta*delta_kronecker(i, l) + alpha* S[i]*res


def Jacobian_IS(state, n, method, alpha):
    return -1*Jacobian_SI(state, n, method, alpha)

def Jacobian(state, n, method, alpha, beta,p):
    J_SS=Jacobian_SS(state,n)
    J_II = Jacobian_II(state, n, method, alpha, beta)
    J_RI = Jacobian_RI(n, beta, p)
    J_IS = Jacobian_IS(state, n, method, alpha)
    J_SI = Jacobian_SI(state, n, method, alpha)
    J_DI = Jacobian_DI(n, beta, p)
    J_SR = np.identity(n)*0
    J_SD = np.identity(n) * 0
    J_IR = np.identity(n)*0
    J_ID = np.identity(n)*0
    J_RS = np.identity(n)*0
    J_DS = np.identity(n)*0
    J_RR = np.identity(n)*0
    J_DR = np.identity(n)*0
    J_RD = np.identity(n)*0
    J_DD = np.identity(n)*0
    Z=[J_SS, J_SI, J_SR, J_SD]
    J =  = np.identity(n*4)
    for i in range(n):
        for j in range(n):
            J[i][j]=J_SS[i][j]
            J[i+n][j]= J_SI[i][j]
            J[i + 2*n][j] = J_SR[i][j]
            J[i + 3 * n][j] = J_SD[i][j]
            J[i][j+ n] = J_IS[i][j]
            J[i][j + 2* n] = J_RS[i][j]
            J[i][j + 3 * n] = J_DS[i][j]
            J[i + n][j + n] = J_II[i][j]
            J[i + n][j + 2 * n] = J_RI[i][j]
            J[i + n][j + 3 * n] = J_DI[i][j]
            J[i + 2 * n][j + n] = J_IR[i][j]
            J[i + 2 * n][j + 2*n] = J_RR[i][j]
            J[i + 2 * n][j + 3 * n] = J_DR[i][j]
            J[i + 3 * n][j +  n] = J_ID[i][j]
            J[i + 3 * n][j + 2*n] = J_RD[i][j]
            J[i + 3 * n][j + 3 * n] = J_DD[i][j]
    return J




####Fehlerfotpflanzung
def error_propagation(sol, init_error):
    sigma_F = np.diag(init_error, k=0)
    for t in range(len(sol)):
        J = Jacobian(sol[t])
        sigma_F = J*sigma_F*J.T
    return sigma_F


##Todo:
# define / find effective N_eff, I_eff
# import P
#define N(i)