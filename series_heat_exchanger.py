"""
This code use the laws of thermodynamics and heat transfer to simulate the temperature profile in a double tube series heat exchanger:
assumptions:
-Density is constant throughout the exchanger for all streams.
-The current of the inner tube is the cold one
"""
from numpy import linspace, log

class Stream:
    def __init__(self, Q, rho, M, T0, Viscosity, Cp, Kf):
        self.W = Q*rho/60 #Lb/s
        self.T0 = T0
        self.Viscosity = Viscosity
        self.Cp = Cp
        self.Kf = Kf
        self.rho = rho*7.480519350799
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

viscosity_parameters =  [-10.2158,1792.5,0.01773,-0.000012631]
Cp_parameters        =  [92.053,-0.039953,-0.00021103,5,3469E-07]
Kf_parameters        =  [-0.2758,0.004612,-5.5391E-06]
parameters_pipe      =  [228.2103,0.057999,-0.000086806]

#Nom = class ( W,   Rho,       MW,  T0, viscosity_parameters, Cp_parameters, Kf_parameters)
cold = Stream(40, 8.337, 18.01528,  80, viscosity_parameters, Cp_parameters, Kf_parameters)
hot1 = Stream(30, 8.337, 18.01528, 120, viscosity_parameters, Cp_parameters, Kf_parameters)
hot2 = Stream(30, 8.337, 18.01528, 200, viscosity_parameters, Cp_parameters, Kf_parameters)

Do = 2.37/12 #ft
Di = (2.37 - 2*0.154)/12 #ft
Ds = 3.50/12 #ft
Deq = (Ds**2 - Do**2)/Do #Ft

def calc_kw (T,parameters_pipe = [228.2103,0.057999,-0.000086806]):
    """calculates the thermal conductivity (Kw) in BTU/(Lb-°F) given a temperature in °F"""
    T = ((T-32)*5/9) + 273.15 #This pass the T from °F to K
    [A,B,C] = parameters_pipe
    Kw = (A + B*T + C*T**2)*0.577789316545299 #BTU/(ft-°F)
    return Kw

def calc_Uo(Ts, streams, Diameters, Areas):
    [T,t] = Ts
    [Area,area] = Areas
    [Deq,Di] = Diameters
    [cold,hot] = streams
    
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
    #creo que acá va el while:
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
        Error = max([abs((Tw-Tw_n)/Tw_n), abs((tw-tw_n)/tw_n)])*100
        Tw = Tw_n
        tw = tw_n
        if i == maxitr:
            print(f'HPT CORRA HICIMOS {maxitr} ITERACIONES Y OBTUVIMOS\n    Tw:{Tw}\n   tw:{tw}')
            return Uo
    return Uo


    


    


