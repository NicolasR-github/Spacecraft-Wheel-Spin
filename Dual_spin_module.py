
import matplotlib.pyplot as plt

def plotting(t,sol):
    plt.plot(t,sol[:,0],'y',label='om1(t)')
    plt.plot(t,sol[:,1],'m',label='om2(t)')
    plt.plot(t,sol[:,2],'orange',label='om3(t)')
    plt.legend(loc='best')
    plt.xlabel('t')
    plt.grid()
    plt.show()


def rota(y,t,omw,tup,I1,I2,I3,Iw):
    #tup: duration of continuous increase of omw
    om1,om2,om3=y
    k=omw/tup
    dydt=[(1/I1)*((I2-I3)*om2*om3+Iw*om3*k*t),(1/I2)*((I3-I1)*om1*om3-Iw*k),(1/I3)*((I1-I2)*om1*om2-Iw*om1*k*t)]
    return dydt

def rota_up(y,t,omw,I1,I2,I3,Iw):
    om1,om2,om3=y
    dydt=[(1/I1)*((I2-I3)*om2*om3+Iw*om3*omw),(1/I2)*((I3-I1)*om1*om3),(1/I3)*((I1-I2)*om1*om2-Iw*om1*omw)]
    return dydt
