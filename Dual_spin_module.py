
import matplotlib.pyplot as plt
import numpy as np

def plotting(t,sol):
    plt.plot(t,sol[:,0],'y',label='om1(t)')
    plt.plot(t,sol[:,1],'m',label='om2(t)')
    plt.plot(t,sol[:,2],'orange',label='om3(t)')
    plt.legend(loc='best')
    plt.xlabel('t')
    plt.grid()
    plt.show()


def rota(y,t,omw,tup,I1,I2,I3,Iw):
    #here, the wheel rotation speed omw(t) is increasing with t (omw(t) = k*t)
    #function following the variations of w
    #tup: duration of continuous increase of omw
    om1,om2,om3=y
    k=omw/tup
    dydt=[(1/I1)*((I2-I3)*om2*om3+Iw*om3*k*t),(1/I2)*((I3-I1)*om1*om3-Iw*k),(1/I3)*((I1-I2)*om1*om2-Iw*om1*k*t)]

    return dydt

def rota_up(y,t,omw,I1,I2,I3,Iw):
    #here, the wheel rotation speed has reached for its commanded value, omw 
    #function following the variations of w
    om1,om2,om3=y
    dydt=[(1/I1)*((I2-I3)*om2*om3+Iw*om3*omw),(1/I2)*((I3-I1)*om1*om3),(1/I3)*((I1-I2)*om1*om2-Iw*om1*omw)]

    return dydt

def rota_beta(y,t,omw,tup,I1,I2,I3,Iw):
    #here, the wheel rotation speed omw(t) is increasing with t (omw(t) = k*t)
    #function following the variations of w and of the quaternions (ie Euler) parameters
    #tup: duration of continuous increase of omw
    om1,om2,om3,b0,b1,b2,b3=y
    
    BETA=np.mat([[b0,-b1,-b2,-b3],[b1,b0,-b3,b2],[b2,b3,b0,-b1],[b3,-b2,b1,b0]])
    #Quaternions derivatives:
    qd=BETA*np.matrix.transpose(np.mat([0,om1,om2,om3]))
    db0,db1,db2,db3=np.double(np.array(qd)[0]),np.double(np.array(qd)[1]),np.double(np.array(qd)[2]),np.double(np.array(qd)[3])

    k=omw/tup
    dy_b_dt=[(1/I1)*((I2-I3)*om2*om3+Iw*om3*k*t),(1/I2)*((I3-I1)*om1*om3-Iw*k),(1/I3)*((I1-I2)*om1*om2-Iw*om1*k*t),db0,db1,db2,db3] 
   
    return dy_b_dt

def rota_up_beta(y,t,omw,I1,I2,I3,Iw):
    #here, the wheel rotation speed has reached for its commanded value, omw 
    #function following the variations of w
    om1,om2,om3,b0,b1,b2,b3=y

    BETA=np.mat([[b0,-b1,-b2,-b3],[b1,b0,-b3,b2],[b2,b3,b0,-b1],[b3,-b2,b1,b0]])
    #Quaternions derivatives:
    qd=BETA*np.matrix.transpose(np.mat([0,om1,om2,om3]))
    db0,db1,db2,db3=np.double(np.array(qd)[0]),np.double(np.array(qd)[1]),np.double(np.array(qd)[2]),np.double(np.array(qd)[3])


    dy_b_dt=[(1/I1)*((I2-I3)*om2*om3+Iw*om3*omw),(1/I2)*((I3-I1)*om1*om3),(1/I3)*((I1-I2)*om1*om2-Iw*om1*omw),db0,db1,db2,db3]

    return dy_b_dt
