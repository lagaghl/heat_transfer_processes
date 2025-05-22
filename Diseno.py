from numpy import log, pi, linspace, vstack,exp, log10, ones
from scipy.integrate import solve_ivp, quad
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

class Aceite:
    def __init__(self):
        self.W = 3200 # lb/h
        self.T1 = 260
        self.T2 = 130
        self.SG = 0.865

    def rho(self, T):
        return 62.4*self.SG * (1 - 0.00061 * (T - 60))
    def miu(self, T):
        return 0.00488 * exp(4230.8/(T + 459.67)) * 2.4191 # lb/(ft-h)
    def kf(self, T):
        return -3e-5*T + 0.0794
    def Cp(self, T):
        return 6e-4*T + 0.4062
    def Pr(self, T):
        miu = self.miu(T)
        kf = self.kf(T)
        Cp = self.Cp(T)
        return miu * Cp / kf

class Agua:
    def __init__(self):
        self.t1 = 80
        self.t2 = 120

    def rho(self, T):
        T = (T - 32)*5/9 + 273.15 # K
        return 0.3470 * 0.274**(-(1-T/647.13)**0.2857) * 1000 * 0.062428 # lb/ft3
    
    def miu(self, T):
        T = (T - 32)*5/9 + 273.15 # K
        logmiu = -10.2158 + 1792.5/T + 0.01773*T - 1.263e-5*T**2 
        return 10**logmiu * 2.4191 # lb/(ft-h)
    
    def kf(self, T):
        T = (T - 32)*5/9 + 273.15
        return (-0.2758+0.004612*T-5.539e-6*T**2) * 0.57782 # Btu/(h-ft-°F)
    
    def Cp(self, T):
        T = (T - 32)*5/9 + 273.15
        return (92.053-3.9953e-2*T-2.1103e-4*T**2+5.3469e-7*T**3) * 1/18 * 0.238846 # Btu/(h-°F)
    
    def Pr(self, T):
        miu = self.miu(T)
        kf = self.kf(T)
        Cp = self.Cp(T)
        return miu * Cp / kf

aceite = Aceite()
agua = Agua()
agua.W = aceite.W * quad(lambda T:aceite.Cp(T),aceite.T2,aceite.T1)[0] / quad(lambda T:agua.Cp(T),agua.t1,agua.t2)[0]

kw = 9.418  #BTU/(h-ft-°F)
kens = 0.0008 # Btu/(h-ft-°F)

'''Dimensiones del intercambiador de calor'''
'''Ds 1.050'''

# Do = 0.920/12 #ft
# Di = Do - (2*0.065)/12 #ft
# Ds = 1.050/12 #ft

# Do = 0.884/12 #ft
# Di = Do - (2*0.083)/12 #ft
# Ds = 1.050/12 #ft

# Do = 0.884/12 #ft
# Di = Do - (2*0.083)/12 #ft
# Ds = 1.050/12 #ft

# Do = 0.824/12 #ft
# Di = Do - (2*0.113)/12 #ft
# Ds = 1.050/12 #ft
##L  = 7.140099999996698 ft
##caida de Presión anulo = 7.276869190760365 ft
##caida de Presión tubo = 6.558802679197342 ft

# Do = 0.742/12 #ft
# Di = Do - (2*0.154)/12 #ft
# Ds = 1.050/12 #ft

# Do = 0.612/12 #ft
# Di = Do - (2*0.219)/12 #ft
# Ds = 1.050/12 #ft

'''Ds 1.315'''
# Do = 1.185/12 #ft
# Di = Do - (2*0.065)/12 #ft
# Ds = 1.315/12 #ft

# Do = 1.097/12 #ft
# Di = Do - (2*0.109)/12 #ft
# Ds = 1.315/12 #ft
##L  = 9.986699999990064 ft
##caida de Presión anulo = 7.297075814083712 ft
##caida de Presión tubo = 1.5955410044895757 ft

# Do = 1.049/12 #ft
# Di = Do - (2*0.133)/12 #ft
# Ds = 1.315/12 #ft
##L  = 9.000099999992363 ft
##caida de Presión anulo = 3.7418519148099563 ft
##caida de Presión tubo = 2.409801356011151 ft

# Do = 0.957/12 #ft
# Di = Do - (2*0.179)/12 #ft
# Ds = 1.315/12 #ft
##L  = 7.353099999996202 ft
##caida de Presión anulo = 1.3397504817602766 ft
##caida de Presión tubo = 6.450887898704714 ft

# Do = 0.815/12 #ft
# Di = Do - (2*0.250)/12 #ft
# Ds = 1.315/12 #ft

'''Ds 1.660'''

# Do = 1.442/12 #ft
# Di = Do - (2*0.109)/12 #ft
# Ds = 1.660/12 #ft
## L  = 14.619099999979268 ft
## caida de Presión anulo = 6.892496871151998 ft
## caida de Presión tubo = 0.5100124343901555 ft

# Do = 1.380/12 #ft
# Di = Do - (2*0.140)/12 #ft
# Ds = 1.660/12 #ft
##L  = 12.85009999998339 ft
##caida de Presión anulo = 2.955595513416675 ft
##caida de Presión tubo = 0.725276660727566 ft

# Do = 1.278/12 #ft
# Di = Do - (2*0.191)/12 #ft
# Ds = 1.660/12 #ft
##L  = 10.383699999989139 ft
##caida de Presión anulo = 0.9953255350060957 ft
##caida de Presión tubo = 1.4695849639467085 ft

# Do = 1.160/12 #ft
# Di = Do - (2*0.250)/12 #ft
# Ds = 1.660/12 #ft
##L  = 8.07929999999451 ft
##caida de Presión anulo = 0.36981500777760345 ft
##caida de Presión tubo = 4.446797264376464 ft

# Do = 0.896/12 #ft
# Di = Do - (2*0.382)/12 #ft
# Ds = 1.660/12 #ft

'''Ds 1.9'''
# Do = 1.770/12 #ft
# Di = Do - (2*0.065)/12 #ft
# Ds = 1.9/12 #ft

Do = 1.682/12 #ft
Di = Do - (2*0.109)/12 #ft
Ds = 1.9/12 #ft
##L  = 19.000699999969058 ft
##caida de Presión anulo = 6.977142987448936 ft
##caida de Presión tubo = 0.33206473702811934 ft

# Do = 1.610/12 #ft
# Di = Do - (2*0.145)/12 #ft
# Ds = 1.9/12 #ft
##L  = 16.301899999975348 ft
##caida de Presión anulo = 2.6294343758684846 ft
##caida de Presión tubo = 0.39731835863656073 ft

# Do = 1.5/12 #ft
# Di = Do - (2*0.2)/12 #ft
# Ds = 1.9/12 #ft
##L  = 12.979699999983088 ft
##caida de Presión anulo = 0.8410230456484982 ft
##caida de Presión tubo = 0.7187622628165755 ft

# Do = 1.338/12 #ft
# Di = Do - (2*0.281)/12 #ft
# Ds = 1.9/12 #ft
##L  = 9.269699999991735 ft
##caida de Presión anulo = 0.23490275431191776 ft
##caida de Presión tubo = 2.4382038552809107 ft

'''Ds 3.50'''
# Do = 2.27/12 #ft
# Di = Do - (2*0.154)/12 #ft
# Ds = 3.50/12 #ft
##L  = 37.08130000010744 ft
##caida de Presión anulo = 0.03311098254376077 ft
##caida de Presión tubo = 0.10037814815377288 ft

gc = 32.174*(3600**2) #ft-lbm/(lbf*h^2)

Deq = (Ds**2 - Do**2)/Do #Ft
Area = pi/4 * (Ds**2 - Do**2)
area = pi/4 * Di**2

dL = 0.0001 #ft
T_largo = [aceite.T1]
t_largo = [agua.t2]
L_largo = [0]
Re_largo = []
re_largo = []
dP_anulo_largo = []
dP_tubo_largo = []
P_anulo_largo = [50]
P_tubo_largo = [50]
P_turbulento_tubo = 0
P_turbulento_anulo = 0
dP_turbulento_anulo = 0
dP_turbulento_tubo = 0
def calc_Uo(Ts):
    global P_turbulento_tubo; global P_turbulento_anulo; global dP_turbulento_anulo;global dP_turbulento_tubo
    [T,t] = Ts
    #aceite stream:
    G = aceite.W/area # Gasto masico lb/h-ft2
    Miu = aceite.miu(T)
    Re = G*Di/Miu
    #agua stream:
    g = agua.W/Area
    miu = agua.miu(t)
    re = g*Deq/miu
    #Start with the assumption Tw = T and tw = t
    [Tw,tw] = calc_Tw(T,t)

    Miu_W = aceite.miu(Tw)
    miu_w = agua.miu(tw)

    """Calculo de caida de presión"""
    re_fluido = g*(Ds-Do)/miu
    Re_largo.append(Re)
    re_largo.append(re_fluido)

    if Re < 2300:
        if P_turbulento_tubo == 0 and dP_turbulento_tubo == 0 and dP_tubo_largo != []:
            P_turbulento_tubo = P_tubo_largo[-1]
            dP_turbulento_tubo = dP_tubo_largo[-1]
        fd_tubo = 64/Re
        n_tubo = 0.25
    else:
        fd_tubo  = fsolve(lambda f: colebrook(f,Re,1.64042e-6,Di) - f, 0.02)[0]
        n_tubo = 1.4
    if re < 2300:
        if P_turbulento_anulo == 0 and dP_turbulento_anulo == 0 and dP_anulo_largo != []:
            P_turbulento_anulo = P_tubo_largo[-1]
            dP_turbulento_anulo = dP_anulo_largo[-1]
        fd_anulo = 64/re
        n_anulo = 0.25
    else:
        fd_anulo = fsolve(lambda f: colebrook(f,re_fluido,1.64042e-6,Ds-Do) - f, 0.02)[0]
        n_anulo = 0.14

    dP_anulo = (fd_anulo*g**2*(L_largo[-1]+dL))/(2*gc*agua.rho(t)*(Ds-Do)*((miu/miu_w)**n_anulo)) * 1/144
    dP_tubo  = (fd_tubo *G**2*(L_largo[-1]+dL))/(2*gc*aceite.rho(T)*Di*((Miu/Miu_W)**n_tubo)) * 1/144


    if P_turbulento_anulo != 0:
        P_anulo_largo.append(P_turbulento_anulo-dP_anulo)
        dP_anulo_largo.append(dP_turbulento_anulo+dP_anulo)
    else:
        P_anulo_largo.append(P_anulo_largo[0]-dP_anulo)
        dP_anulo_largo.append(dP_anulo)
    
    if P_turbulento_tubo != 0:
        P_tubo_largo.append(P_turbulento_tubo-dP_tubo)
        dP_tubo_largo.append(dP_turbulento_tubo+dP_tubo)
    else:
        P_tubo_largo.append(P_tubo_largo[0]-dP_tubo)
        dP_tubo_largo.append(dP_tubo)

    if Re < 2300:
        Nu = 1.86*(Re*aceite.Pr(T)*(Di/(L_largo[-1]+dL)))**(1/3)*(Miu/Miu_W)**0.14
    elif Re > 100000:
        Nu = 0.027*(Re**0.8)*(aceite.Pr(T)**(1/3))*(Miu/Miu_W)**0.14
    else:
        Nu = 0.116*(1+ Di/(L_largo[-1]+dL))**(2/3) * (Re**(2/3) - 125) * (aceite.Pr(T)**(1/3)) * (Miu/Miu_W)**0.14
    if re < 2300:
        nu = 1.86*(re*agua.Pr(t)*(Deq/(L_largo[-1]+dL)))**(1/3)*(miu/miu_w)**0.14
    elif Re > 100000:
        nu = 0.027*(re**0.8)*(agua.Pr(T)**(1/3))*(miu/miu_w)**0.14
    else:
        nu = 0.116*(1+ Deq/(L_largo[-1]+dL))**(2/3) * (Re**(2/3) - 125) * (agua.Pr(T)**(1/3)) * (miu/miu_w)**0.14
    
    ho = agua.kf(T)*Nu/Deq
    hio = (aceite.kf(T)*nu/Di)*(Di/Do)

    Uo = 1/(1/hio + Do * log(Do/Di)/(2*kw) + 1/ho + kens*(Do/Di))
    return Uo

def calc_Tw(T,t):
    #aceite stream:
    G = aceite.W/area # Gasto masico lb/h-ft2
    Miu = aceite.miu(T)
    Re = G*Di/Miu
    #agua stream:
    g = agua.W/Area
    miu = agua.miu(t)
    re = g*Deq/miu

    #Start with the assumption Tw = T and tw = t
    Tw = T
    tw = t
    tol = 1e-6
    Error = 1
    maxitr = 500
    i = 0
    while Error>tol:
        i += 1
        Miu_W = aceite.miu(Tw)
        miu_w = agua.miu(tw)

        if Re < 2300:
            Nu = 1.86*(Re*aceite.Pr(T)*(Di/(L_largo[-1]+dL)))**(1/3)*(Miu/Miu_W)**0.14
        elif Re > 100000:
            Nu = 0.027*(Re**0.8)*(aceite.Pr(T)**(1/3))*(Miu/Miu_W)**0.14
        else:
            Nu = 0.116*(1+ Di/(L_largo[-1]+dL))**(2/3) * (Re**(2/3) - 125) * (aceite.Pr(T)**(1/3)) * (Miu/Miu_W)**0.14
        if re < 2300:
            nu = 1.86*(re*agua.Pr(t)*(Deq/(L_largo[-1]+dL)))**(1/3)*(miu/miu_w)**0.14
        elif Re > 100000:
            nu = 0.027*(re**0.8)*(agua.Pr(T)**(1/3))*(miu/miu_w)**0.14
        else:
            nu = 0.116*(1+ Deq/(L_largo[-1]+dL))**(2/3) * (Re**(2/3) - 125) * (agua.Pr(T)**(1/3)) * (miu/miu_w)**0.14
        
        ho = agua.kf(T)*Nu/Deq
        hio = (aceite.kf(T)*nu/Di)*(Di/Do)

        Uo = 1/(1/hio + Do * log(Do/Di)/(2*kw) + 1/ho + kens*(Do/Di))

        Tw_n = T - 1/hio  * Uo * (T-t)
        tw_n = t + 1/ho * Uo * (T-t)
        Error = max([abs(Tw - Tw_n), abs(tw - tw_n)])*100
        Tw = Tw_n
        tw = tw_n
        if i == maxitr:
            print("Warning: Iterative process did not converge. Returning fallback for Tw and tw.")
            return [Tw,tw]
    return [Tw,tw]

def colebrook(f, Re, e, D):
    # Colebrook-White equation for turbulent flow
    return (-2*log10( e/(3.71*D) + 2.51/(Re*(f)**0.5)))**(-2)

while T_largo[-1] >= aceite.T2 or t_largo[-1] >= agua.t1:

    T = T_largo[-1]     
    t = t_largo[-1]
    Uo = calc_Uo([T,t])

    Tn = T - dL * Uo*pi*Do*(T - t) / (aceite.rho(T)*aceite.Cp(T))
    tn = t - dL * Uo*pi*Do*(T - t) / (agua.rho(t)*agua.Cp(t))
    
    T_largo.append(Tn)
    t_largo.append(tn)
    L_largo.append((L_largo[-1]+dL) + dL)

# Plotting
plt.plot(L_largo, T_largo, label='Aceite')
plt.plot(L_largo, t_largo, label='Agua')
plt.xlabel('Longitud (ft)')
plt.ylabel('Temperatura (°F)')
plt.title('Perfil de temperatura en el intercambiador de calor')
plt.legend()
plt.grid()
plt.show()

plt.plot(L_largo[1::], Re_largo, label='reynolds aceite')
plt.plot(L_largo[1::], re_largo, label='reynolds agua')
plt.xlabel('Longitud (ft)')
plt.ylabel('Reynolds')
plt.title('Perfil de Reynolds en el intercambiador de calor')
plt.legend()
plt.grid()
plt.show()

plt.plot(L_largo[1::], dP_anulo_largo, label='dP anulo')
plt.plot(L_largo[1::], dP_tubo_largo, label='dP tubo')
plt.xlabel('Longitud (ft)')
plt.ylabel('dP (psi)')
plt.title('Perfil de perdida de carga diferencial en el intercambiador de calor')
plt.legend()
plt.grid()
plt.show()

plt.plot(L_largo, P_anulo_largo, label='P anulo')
plt.plot(L_largo, P_tubo_largo, label='P tubo')
plt.xlabel('Longitud (ft)')
plt.ylabel('P (psi)')
plt.title('Perfil de perdida de carga acumulada en el intercambiador de calor')
plt.legend()
plt.grid()
plt.show()

temp_matrix = vstack([t_largo,T_largo,t_largo])

plt.figure(figsize=(10, 3))
plt.imshow(temp_matrix, aspect='auto', cmap='jet', extent=[L_largo[0], (L_largo[-1]+dL), 0, 3])

plt.yticks([0.5, 1.5, 2.5], ['Fluido caliente', 'Fluido frío', 'Fluido caliente'])

plt.colorbar(label='Temperatura (°F)')
plt.xlabel('largo del intercambiador (ft)')
plt.title('Perfil de temperatura en el intercambiador de calor')
plt.tight_layout()
plt.show()

print(f'T2 = {T_largo[-1]} °F') 
print(f't1 = {t_largo[-1]} °F')
print(f'L  = {(L_largo[-1]+dL)} ft')
print(f'caida de Presión anulo = {P_anulo_largo[0]-P_anulo_largo[-1]} ft')
print(f'caida de Presión tubo = {P_tubo_largo[0]-P_tubo_largo[-1]} ft')


