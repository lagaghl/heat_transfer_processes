from numpy import pi,exp, log, linspace, vstack
from scipy.integrate import odeint as ode
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

global D0; global Di; global k; global kf1; global kf2; global vel1; global vel2; global rho1; global rho2
D0 = 0.06
Di = 0.04
k  = 15
kf1 = 0.6
kf2 = 0.12
vel1 = 0.8
vel2 = 1.5
rho1 = 997
rho2 = 850

def odes(Ts,L):
    t = Ts[0]
    T = Ts[1]
    [Cp1,Cp2] = calc_Cp(t,T)
    U0 = calc_U0(t,T,[Cp1,Cp2])
    
    w1 = vel1*rho1*pi/4*Di**2
    w2 = vel2*rho2*pi/4*D0**2
    dt_dL =  (U0*pi*D0*(T-t))/(w1*Cp1)
    dT_dL = -(U0*pi*D0*(T-t))/(w2*Cp2)
    return [dt_dL,dT_dL]

def calc_U0(t,T,Cps):
    [Cp1,Cp2] = calc_Cp(t,T)
    [miu1,miu2] = calc_Miu(t,T)
    Re1 = rho1 * vel1 * Di/miu1
    Re2 = rho2 * vel2 * D0/miu2
    Pr1 = miu1 * Cp1 / kf1
    Pr2 = miu2 * Cp2 / kf2

    [tw,Tw] = fsolve(calc_h,[t,T], args=(t,T,[Re1,Re2],[Pr1,Pr2]))

    #Recalcular hi y h0:
    [miu_w,miu_W] = calc_Miu(tw,Tw)
    [hi,h0] = Sieder_Tate([Re1,Re2],[Pr1,Pr2],[miu1,miu2],[miu_w,miu_W])
    #Calcular U0 a esa t y T:
    U0 = (1/h0 + D0*log(D0/Di)/(2*k) + 1/hi)**-1
    return U0

def Sieder_Tate(Re,Pr,Miu,Miuw):
    [Re1,Re2] = Re
    [Pr1,Pr2] = Pr
    [Miu1,Miu2] = Miu
    [miu_w,miu_W] = Miuw
    h1 = kf1/Di * (0.027* Re1**(4/5) * Pr1**(1/3) * (Miu1/miu_w)**0.14)
    h2 = kf2/D0 * (0.027* Re2**(4/5) * Pr2**(1/3) * (Miu2/miu_W)**0.14)
    return[h1,h2]

def calc_h (Tws,t,T,Re,Pr):
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
    miu1 = exp(a1+b1/t)*1e-3

    a2 = -6.944
    b2 =  2036.8
    miu2 = exp(a2+b2/T)*1e-3

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

L = 3
L_largo = linspace(0,L,500)
res = ode(odes,([0,100]),L_largo)
t_largo = res[:,0]
T_largo = res[:,1]
plt.plot(L_largo,t_largo, label='Fluido frío (t)')
plt.plot(L_largo,T_largo, label='Fluido caliente (T)')
plt.xlabel('Longitud del intercambiador (m)')
plt.ylabel('Temperatura (°C)')
plt.title('Perfil de temperaturas en el intercambiador')
plt.grid()
plt.legend()
plt.show()

temp_matrix = vstack([T_largo, t_largo, T_largo])

# Crear el mapa de calor
plt.figure(figsize=(10, 3))
plt.imshow(temp_matrix, aspect='auto', cmap='jet', extent=[L_largo[0], L_largo[-1], 0, 3])

# Etiquetas del eje Y
plt.yticks([0.5, 1.5, 2.5], ['Fluido caliente', 'Fluido frío', 'Fluido caliente'])

plt.colorbar(label='Temperatura (°C)')
plt.xlabel('Longitud del intercambiador (m)')
plt.title('Distribución de temperaturas en el intercambiador de calor')
plt.tight_layout()
plt.show()