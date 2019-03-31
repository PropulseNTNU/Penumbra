import numpy as np
#------------------------------------Kalmanfilter-----------------------
#Initiates matrices needed in the kalman filter
C_d=np.matrix([1,0])
#const Matrix<2,2> E_d={0.00001,0.0000002,0.000000004,0.0000003};
Q=np.matrix([[300,0],[0,1]])
R=[1];
I=np.matrix([[1,0],[0,1]])
P_k_bar=np.matrix([[1,0],[0,0.1]])
x_hat_bar=np.matrix([[0],[0]])#skal egentlig være ca.{2000,300};
#x_hat = np.array([[0],[0]])
#K_k=np.array([[0],[0]])
#P_k=np.matrix([[0,0],[0,0]])
#estimates=[0,0]
def kalman(altitude, acceleration, dt):
    global P_k_bar#=np.matrix([[1,0],[0,0.1]])
    global x_hat_bar#=np.array([[0],[280]])

  #Updating variables
    A_d=np.matrix([[1,dt],[0,1]])
    B_d=np.matrix([[0],[dt]])
  #Computing kalman gain------------------------------------------------------
    K_k = P_k_bar*(C_d.T)*((C_d*P_k_bar*(C_d.T)+R).I)
  #Serial << "K_k" << K_k << '\n';
  #Update estimate with measurement-------------------------------------------
    x_hat = x_hat_bar + K_k*(altitude -(x_hat_bar[0][0]))
  #Serial << "x_hat: " << x_hat << '\n';
  #//updatet estimates----------------------------------------------------------
    #estimates[0]=x_hat[0][0]#kanskje ikke bruke pekere
    #estimates[1]=x_hat[1][0]
  #//Compute  error  covariance  for  updated  estimate-------------------------
    P_k = (I - K_k*C_d)*P_k_bar *((I-K_k*C_d).T) + K_k*R*(K_k.T)
  #Serial << "P_k: " << P_k << '\n';
  #//project ahead--------------------------------------------------------------
    x_hat_bar = A_d * x_hat + B_d * (acceleration)
    P_k_bar = A_d * P_k * (A_d.T) + Q #//Mulig Q skal byttes ut med: E_d * Q * (~E_d); Uten pådrag blir det bare Q;
  #Serial << "P_k_bar: " << P_k_bar << '\n';
  #Serial << "x_hat_bar: " << x_hat_bar << '\n';
    return x_hat[0][0], x_hat[1][0]




#-----------
