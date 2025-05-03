#This code use the laws of thermodynamics and heat transfer to simulate the temperature profile in a double tube series heat exchanger
from numpy import linspace

class Liquid:
    def __init__(self, Q, rho, M, T0, Viscosity, Cp, Kf):
        self.W = Q*rho
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
        """calculates the Cp in BTU/(Lb-°F) given a temperature in °F"""
        T = ((T-32)*5/9) + 273.15 #This pass the T from °F to K
        [A,B,C] = self.Kf
        Kf = (A + B*T + C*T**2)*0.577789316545299 #BTU/(ft-°F)
        return Kf

viscosity_parameters =  [-10.2158,1792.5,0.01773,-0.000012631]
Cp_parameters = [92.053,-0.039953,-0.00021103,5,3469E-07]
Kf_parameters = [-0.2758,0.004612,-5.5391E-06]
#Nom = class ( W,   Rho,       MW,  T0, viscosity_parameters, Cp_parameters, Kf_parameters)
Cold = Liquid(40, 8.337, 18.01528,  80, viscosity_parameters, Cp_parameters, Kf_parameters)
Hot1 = Liquid(30, 8.337, 18.01528, 120, viscosity_parameters, Cp_parameters, Kf_parameters)
Hot2 = Liquid(30, 8.337, 18.01528, 200, viscosity_parameters, Cp_parameters, Kf_parameters)
