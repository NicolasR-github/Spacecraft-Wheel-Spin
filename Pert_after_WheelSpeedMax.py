
# coding: utf-8


#Let's visualize the stabilization of a spacecraft
#by a spinning wheel, rotating at a constant speed, in case a perturbation intervenes after the spinning wheel reaches for its 
#constant target speed.

#the Ii are the inertial moments around the principal axis of the spacecraft
#Iw is the spin inertial moment of the wheel (whose, here, is chosen along the 2nd principal axis of the spacecraft)
#omw is the rotation speed of the wheel

#The "full" equations are [d(om1)/dt, d(om2)/dt, d(om3)/dt] =  
#[(1/I1)*((I2-I3)*om2*om3+Iw*om3*k*t),(1/I2)*((I3-I1)*om1*om3-Iw*k),(1/I3)*((I1-I2)*om1*om2-Iw*om1*k*t)]
#and the "simplified" ones, ie when the spinning wheel is at its 
#constant target speed, are [d(om1)/dt, d(om2)/dt, d(om3)/dt] =  
#[(1/I1)*((I2-I3)*om2*om3+Iw*om3*omw),(1/I2)*((I3-I1)*om1*om3),(1/I3)*((I1-I2)*om1*om2-Iw*om1*omw)]

import numpy as np

from scipy.integrate import odeint

import Dual_spin_module as wheel

global I1,I2,I3,Iw,omw,tup


#for example
I1,I2,I3=630000,1365000,1665000
Iw=96
#For a stabilization along the intermediate (here the second one) principal axis, with om2
#chosen = np.pi/50 rd.s-1, the value of a positive omw would be, if determined with the simplification mentionned above, 
#at least (I2/I3)*om2*(I3-I2)/Iw > 160 rd.s-1, and
#the value of a negative omw, would be less than (I2/I1)*(I1-I2)/Iw < -1043 rd.s-1

#duration to reach for omw
tup=15

print('Case of stabilization with omw > 160 rd.s-1. 200 rd.s-1 for example')

omw=200

y0=[0,np.pi/50,0]

t1=np.linspace(0,15,1000)


sol=odeint(wheel.rota,y0,t1,(omw,tup,I1,I2,I3,Iw))

wheel.plotting(t1,sol)

#Of course, om1 (=sol[999,0]) and om3 (=sol[999,2]) are left until now equal to 0 rad.s-1
#Then let's change [sol[999,0],sol[999,1],sol[999,2]] = [0,sol[999,1],0] into
#[np.pi/400, sol[999,1]+np.pi/350, np.pi/500], introducing arbitrary slight perturbations

y0=[np.pi/400,sol[999,1]+np.pi/350,np.pi/500]

t2=np.linspace(15,1500,100000)


sol=odeint(wheel.rota_up,y0,t2,(omw,I1,I2,I3,Iw))

wheel.plotting(t2,sol)

print('Case of unrealized stabilization because of an incorrect value of wheel rotation speed\
omw < 160 rd.s-1 and > -1042 rd.s-1; omw = 150 rad.s-1 for example')

omw=150


y0=[0,np.pi/50,0]

t1=np.linspace(0,15,1000)

sol=odeint(wheel.rota,y0,t1,(omw,tup,I1,I2,I3,Iw))

wheel.plotting(t1,sol)

#Of course, om1 (=sol[999,0]) and om3 (=sol[999,2]) are left until now equal to 0 rad.s-1
#Then let's change [sol[999,0],sol[999,1],sol[999,2]] = [0,sol[999,1],0] into
#[np.pi/400,sol[999,1]+np.pi/350,np.pi/500], introducing arbitrary slight perturbations

y0=[np.pi/400,sol[999,1]+np.pi/350,np.pi/500]

t2=np.linspace(15,1500,100000)

sol=odeint(wheel.rota_up,y0,t2,(omw,I1,I2,I3,Iw))

wheel.plotting(t2,sol)


print('Case of stabilization with omw < -1043 rd.s-1. -1250 rd.s-1 for example')

omw=-1250

y0=[0,np.pi/50,0]

t1=np.linspace(0,15,1000)

sol=odeint(wheel.rota,y0,t1,(omw,tup,I1,I2,I3,Iw))

wheel.plotting(t1,sol)

#Of course, om1 (=sol[999,0]) and om3 (=sol[999,2]) are left until now equal to 0 rad.s-1
#Then let's change [sol[999,0],sol[999,1],sol[999,2]] = [0,sol[999,1],0] into
#[np.pi/400,sol[999,1]+np.pi/350,np.pi/500], introducing arbitrary slight perturbations

y0=[np.pi/400,sol[999,1]+np.pi/350,np.pi/500]

t2=np.linspace(15,1500,100000)

sol=odeint(wheel.rota_up,y0,t2,(omw,I1,I2,I3,Iw))

wheel.plotting(t2,sol)

