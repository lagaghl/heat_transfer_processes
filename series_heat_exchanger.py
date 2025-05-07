"""
This code use the laws of thermodynamics and heat transfer to simulate the temperature profile in a double tube series heat exchanger:
assumptions:
-Density is constant throughout the exchanger for all streams.
-The current of the inner tube is the cold one.
-The current of the outer tube is the hot one.
-The service fluid is the same for both currents (The hot ones). 
"""
from numpy import log, pi, linspace, zeros, vstack, column_stack, repeat, newaxis
from scipy.integrate import solve_ivp, quad
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

class Stream:
    def __init__(self, Q, rho, M, T0, Viscosity, Cp, Kf):
        self.W = Q*rho/60 #Lb/s
        self.Q = Q*0.1336806/60 #ft3/s
        self.T0 = T0
        self.Viscosity = Viscosity
        self.Cp = Cp
        self.Kf = Kf
        self.rho = rho*7.480519350799 #Lb/ft3
        self.M = M #molecular weight (international units)
    def miu(self,T):
        """calculates the dynamic viscosity in lb/(ft-s) given a temperature in °F"""
        T = ((T-32)*5/9) + 273.15 #This pass the T from °F to K
        [A,B,C,D] = self.Viscosity
        visc = (10**(A + B/T + C*T + D*T**2))/(1488.164) #Lb/(ft-s)
        return visc
    def cp(self,T):
        """calculates the Cp in BTU/(Lb-°F) given a temperature in °F"""
        T = ((T-32)*5/9) + 273.15 #This pass the T from °F to K
        [A,B,C,D] = self.Cp
        Cp = ((A + B*T + C*T**2 + D*T**3)*self.M)/4.1868 #BTU/(Lb-°F)
        return Cp
    def kf(self,T):
        """calculates the thermal conductivity (Kf) in BTU/(Lb-°F) given a temperature in °F"""
        T = ((T-32)*5/9) + 273.15 #This pass the T from °F to K
        [A,B,C] = self.Kf
        Kf = (A + B*T + C*T**2)*0.577789316545299 #BTU/(ft-°F)
        return Kf
    def Pr(self,T):
        """calculates the number of Prandtl given a temperature in °F"""
        return self.cp(T)*self.miu(T)/self.kf(T)
    def calc_deltah(self,T):
        Tref = 120
        h = quad(lambda T:self.cp(T),Tref,T)[0] #BTU/Lb
        return h

def T_mix(Ts,streams):
    T = fsolve(fobj_T,(Ts[0]),args=(Ts,streams))[0]
    return T

def fobj_T (T,Ts,streams):
    [T1,T2] = Ts
    [hot1,hot2] = streams
    w1 = hot1.W
    w2 = hot2.W
    x1 = w1/(w1+w2)
    x2 = w2/(w1+w2)
    h_componente = x1 * hot1.calc_deltah(T1) + x2 * hot2.calc_deltah(T2)
    h_mezcla = x1 * hot1.calc_deltah(T) + x2 * hot2.calc_deltah(T)
    return h_componente-h_mezcla

viscosity_parameters =  [-10.2158,1792.5,0.01773,-0.000012631]
Cp_parameters        =  [92.053,-0.039953,-0.00021103,5.3469E-07]
Kf_parameters        =  [-0.2758,0.004612,-5.5391E-06]
parameters_pipe      =  [228.2103,0.057999,-0.000086806]
#              Q,   rho,        M,  T0, Viscosity,            Cp,            Kf
#Nom = class ( Q,   Rho,       MW,  T0, viscosity_parameters, Cp_parameters, Kf_parameters)
cold = Stream(40, 8.337, 18.01528,  80, viscosity_parameters, Cp_parameters, Kf_parameters)
hot1 = Stream(30, 8.337, 18.01528, 200, viscosity_parameters, Cp_parameters, Kf_parameters)
hot2 = Stream(30, 8.337, 18.01528, 120, viscosity_parameters, Cp_parameters, Kf_parameters)

Do = 2.37/12 #ft
Di = (2.37 - 2*0.154)/12 #ft
Ds = 3.50/12 #ft
Deq = (Ds**2 - Do**2)/Do #Ft
Area = pi/4 * (Ds**2 - Do**2)
area = pi/4 * Di**2

def calc_kw (T,parameters_pipe = parameters_pipe):
    """calculates the thermal conductivity (Kw) in BTU/(Lb-°F) given a temperature in °F"""
    T = ((T-32)*5/9) + 273.15 #This pass the T from °F to K
    [A,B,C] = parameters_pipe
    Kw = (A + B*T + C*T**2)*0.577789316545299 #BTU/(ft-°F)
    return Kw

def calc_Uo(Ts, streams):
    '''Calculates the overall heat transfer coefficient (Uo) in BTU/(h-ft2-°F) given a temperature in °F'''
    [T,t,Q] = Ts
    [hot,cold] = streams
    
    #hot stream:
    Vel = hot.Q/Area
    Miu = hot.miu(T)
    Re = Vel*Deq*hot.rho/Miu
    #cold stream:
    vel = cold.Q/area
    miu = cold.miu(t)
    re = vel*Di*cold.rho/miu

    #Start with the assumption Tw = T and tw = t
    [Tw,tw] = calc_Tw(T,t,hot,cold)
    Tw_prom = (Tw-tw)/2

    Miu_W = hot.miu(Tw)
    miu_w = cold.miu(tw)

    Nu = 0.027*(Re**0.8)*(hot.Pr(T)**(1/3))*(Miu/Miu_W)**0.14
    nu = 0.027*(re**0.8)*(cold.Pr(T)**(1/3))*(miu/miu_w)**0.14
    ho = hot.kf(T)*Nu/Deq
    hio = (cold.kf(T)*nu/Di)*(Di/Do)
    Uo = 1/(1/hio + Do * log(Do/Di)/(2*calc_kw(Tw_prom)) + 1/ho)
    return Uo

def calc_Tw(T,t,hot,cold):
    """Calculates the wall temperature (Tw) and the tube temperature (tw) in °F given a temperature in °F"""
    #hot stream:
    Vel = hot.Q/Area
    Miu = hot.miu(T)
    Re = Vel*Deq*hot.rho/Miu
    #cold stream:
    vel = cold.Q/area
    miu = cold.miu(t)
    re = vel*Di*cold.rho/miu

    #Start with the assumption Tw = T and tw = t
    Tw = T
    tw = t
    tol = 1e-6
    Error = 1
    maxitr = 500
    i = 0
    while Error>tol:
        i += 1
        Tw_prom = (Tw-tw)/2
        Miu_W = hot.miu(Tw)
        miu_w = cold.miu(tw)

        Nu = 0.027*(Re**0.8)*(hot.Pr(T)**(1/3))*(Miu/Miu_W)**0.14
        nu = 0.027*(re**0.8)*(cold.Pr(T)**(1/3))*(miu/miu_w)**0.14
        ho = hot.kf(T)*Nu/Deq
        hio = (cold.kf(T)*nu/Di)*(Di/Do)

        Uo = 1/(1/hio + Do * log(Do/Di)/(2*calc_kw(Tw_prom)) + 1/ho)

        Tw_n = T - 1/ho * Uo * (T-t)
        tw_n = t + 1/hio * Uo * (T-t)
        Error = max([abs(Tw - Tw_n), abs(tw - tw_n)])*100
        Tw = Tw_n
        tw = tw_n
        if i == maxitr:
            print("Warning: Iterative process did not converge. Returning fallback for Tw and tw.")
            return [Tw,tw]
    return [Tw,tw]

def edos(L,Ts,hot,cold):
    """Calculates the differential equations for the heat exchanger."""
    [T,t,Q] = Ts
    streams=[hot,cold]

    Uo = calc_Uo(Ts, streams)
    dt_dL = Uo*pi*Do*(T-t)/(cold.W * cold.cp(t))
    dT_dL = Uo*pi*Do*(T-t)/( hot.W *  hot.cp(T))
    dQ_dL = Uo*pi*Do*(T-t)/3600 #BTU/Hr

    return[dT_dL, dt_dL, dQ_dL]

L = 10
n = 100
L_long = linspace(0,L,n)
T1_long = zeros(n)
T2_long = zeros(n)
t1_long = zeros(n)
t2_long = zeros(n)

def fobj (Tf, hot, cold):
    Tf = Tf[0]
    if hasattr(cold, 'Tf1'):
        T0 = cold.Tf1
    else:
        T0 = cold.T0

    solution = solve_ivp(edos, [0, L], [Tf, T0, 0], t_eval=L_long, args=(hot, cold))
    T_long = solution.y[0]
    return T_long[-1] - hot.T0

#search for Tf_1:
Tf_1 = 150 #°F
hot1.Tf = fsolve(fobj, (Tf_1), args=(hot1,cold), xtol=1e-6)[0]
 
solution = solve_ivp(edos, [0, L], [hot1.Tf, cold.T0,0], t_eval=L_long, args=(hot1, cold))
T1_long = solution.y[0]
t1_long = solution.y[1]
Q1_long = solution.y[2]
cold.Tf1 = t1_long[-1]

Tw_1 = zeros(n)
tw_1 = zeros(n)

for i in range(n):
    T = T1_long[i]
    t = t1_long[i]
    Tw, tw = calc_Tw(T, t, hot1, cold)
    Tw_1[i] = Tw
    tw_1[i] = tw


#search for Tf_2:
Tf_2 = 150 #°F
hot2.Tf = fsolve(fobj, (Tf_2), args=(hot2,cold), xtol=1e-8)[0]
 
solution = solve_ivp(edos, [0, L], [hot2.Tf, cold.Tf1,Q1_long[-1]], t_eval=L_long, args=(hot2, cold))
T2_long = solution.y[0]
t2_long = solution.y[1]
Q2_long = solution.y[2]
cold.Tf2 = t2_long[-1]

Tw_2 = zeros(n)
tw_2 = zeros(n)

for i in range(n):
    T = T2_long[i]
    t = t2_long[i]
    Tw, tw = calc_Tw(T, t, hot2, cold)
    Tw_2[i] = Tw
    tw_2[i] = tw

plt.plot(L_long,t1_long, label='cold current 1')
plt.plot(L_long,T1_long, label='hot current 1')
plt.plot(L_long+10,t2_long, label='cold current 2')
plt.plot(L_long+10,T2_long, label='hot current 2')
plt.xlabel('exchanger length (ft)')
plt.ylabel('temperature (°F)')
plt.title('temperature profile in the exchanger')
plt.grid()
plt.legend()
plt.show()

T_long = column_stack([T1_long, T2_long])
t_long = column_stack([t1_long, t2_long])

plt.plot(L_long,Q1_long, label='heat 1')
plt.plot(L_long+L,Q2_long, label='heat 2')
plt.xlabel('exchanger length (ft)')
plt.ylabel('Q (BTU/Hr)')
plt.title('heat transfer in the exchanger')
plt.grid()
plt.legend()
plt.show()


temp_matrix = vstack([T_long,t_long,T_long])

plt.figure(figsize=(10, 3))
plt.imshow(temp_matrix, aspect='auto', cmap='jet', extent=[L_long[0], L_long[-1]+L, 0, 3])

plt.yticks([0.5, 1.5, 2.5], ['Fluido caliente', 'Fluido frío', 'Fluido caliente'])

plt.colorbar(label='Temperature (°F)')
plt.xlabel('Exchanger length (ft)')
plt.title('Temperature profile in the exchanger')
plt.tight_layout()
plt.show()

print(f'Q  = {Q2_long[-1] } BTU/hr')
print(f'tf = {t2_long[-1]} °F')
    


