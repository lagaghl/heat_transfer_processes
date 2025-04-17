from numpy import pi, exp, log
from scipy.integrate import solve_ivp as ode
from scipy.optimize import fsolve

def odes(L,Ts,D0,w1,w2):
    t = Ts[0]
    T = Ts[1]
    [Cp1,Cp2] = calc_Cp(t,T)
    U0 = calc_U0(t,T,[Cp1,Cp2])
    

    dt_dL =  (U0*pi()*D0*(T-t))/(w1*Cp1)
    dT_dL = -(U0*pi()*D0*(T-t))/(w2*Cp2)

    return [dt_dL,dT_dL]

def calc_U0(t,T,Cps):
    D0 = 1.1
    Di = 0.8
    k  = 1200
    kf1 = 100
    kf2 = 100
    vel1 = 100
    vel2 = 100
    rho1 = 1000
    rho2 = 1000
    Cp1 = 27.5
    Cp2 = 27.5
    Re1 = rho1 * vel1 * Di/miu1
    Re2 = rho2 * vel2 * D0/miu2
    Pr1 = miu1 * Cp1 / kf1
    Pr2 = miu2 * Cp2 / kf2

    [tw,Tw] = fsolve(calc_h,[t,T],[Re1,Re2],[Pr1,Pr2])
    
    #Recalcular hi y h0:
    [miu1,miu2] = calc_Miu(t,T)
    [miu_w,miu_W] = calc_Miu(tw,Tw)
    [hi,h0] = Sieder_Tate([Re1,Re2],[Pr1,Pr2],[miu1,miu2],[miu_w,miu_W])

    #Calcular U0 a esa t y T:
    U0 = (1/h0 + D0*log(D0/Di)/(2*k) + 1/hi)
    return U0

def Sieder_Tate(Re,Pr,Miu,Miuw):
    [Re1,Re2] = Re
    [Pr1,Pr2] = Pr
    pass

def calc_h (Tws,t,T,Re,Pr):
    #Datos:
    D0 = 1.1
    Di = 0.8
    k  = 1200
    kf1 = 100
    kf2 = 100
    vel1 = 100
    vel2 = 100
    rho1 = 1000
    rho2 = 1000
    Cp1 = 27.5
    Cp2 = 27.5
    #Calculos:
    [Re1,Re2] = Re
    [Pr1,Pr2] = Pr

    [tw,Tw] = Tws
    [miu1,miu2] = calc_Miu(t,T)
    [miu_w,miu_W] = calc_Miu(tw,Tw)

    [hi,h0] = Sieder_Tate([Re1,Re2],[Pr1,Pr2],[miu1,miu2],[miu_w,miu_W])

    # Buscar un Tw y tw que haga que se cumplan las ecuaciones siguientes ¿? 
    f1 = tw - (t + (h0*D0/(hi*Di))*(T-Tw))
    f2 = Tw - (((h0*D0)/(hi*Di))*(1+hi*Di*log(D0/Di)/(2*k))*T + t)/(1 +(h0*D0/(hi*Di))+(hi*Di*log(D0/Di)/(2*k)))

    return [f2,f1]

def calc_Miu(t,T):
    t += 273.15
    T += 273.15

    a1 = -6.944
    b1 =  2036.8
    miu1 = exp(a1+b1/t)

    a2 = -6.944
    b2 =  2036.8
    miu2 = exp(a2+b2/T)

    return [miu1, miu2]

def calc_Cp(t,T):
    t += 273.15
    T += 273.15

    A1 = 92.053
    B1 = -0.039953
    C1 = -0.00021103
    D1 = 5.3469e-7
    Cp1 = A1 + B1*t + C1*t**2 + D1*t**3

    A2 = 92.053
    B2 = -0.039953
    C2 = -0.00021103
    D2 = 5.3469e-7
    Cp2 = A2 + B2*T + C2*T**2 + D2*T**3

    return [Cp1,Cp2]