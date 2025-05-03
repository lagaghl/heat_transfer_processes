#This code use the laws of thermodynamics and heat transfer to simulate the temperature profile in a double tube series heat exchanger
from numpy import linspace

class liquid:
    def __init__(self, W, rho, M, Viscosity, Cp, Kf):
        self.W = W
        self.Viscosity = Viscosity
        self.Cp = Cp
        self.Kf = Kf
        self.rho = rho
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
    
