#Class for calculation the thermopysical Properties of moist air
from numpy import array, arange, empty

class AirMix():
    def __init__(self, pressure=1.01325e5):
        #These parameters are valid around atmospheric pressure
        self.para = {}
        self.para[('air','cp')] = array([1034.09, 2.84887e-1, 4.97078e-7, 1.07702e-10])
        self.para[('air','mu')] = array([1.43387e-6, 6.56244e-8, 2.9905e-11])
        self.para[('air','lambda')] = array([6.69881e-4, 9.42482e-5, 3.2745e-8])
        self.para[('vapor','cp')] = array([3387.781, -5.55984, 1.94105e-2])
        self.para[('vapor','mu')] = array([9.7494e-7, 3.59061e-8, 2.41612e-13])
        self.para[('vapor','lambda')] = array([-3.5376e-3, 6.54755e-5, 1.74460e-8])
        self.molMassWater = 18.01525
        self.molMassAir = 28.949
        #Gas constant in J/(K*kmol)
        self.R = 8314.462175
        
        self.pressure = pressure

    def getSatPressure(self, T):
        #Calculate the saturation tressure of water in Pascal, T given in Kelvin
        A, B, C = 9.6543, 1435.264, -64.848
        return 10**(A-B/(T+C))

    def getMolByPhi(self, T, phi, pressure=None):
        #Calculates mol_water per mol_airmix based on T in K and RH
        #relative humidity should be smaller than 1, if given in percent
        if pressure is None:
            pressure=self.pressure
        if phi > 1: phi/=100
        return phi*self.getSatPressure(T)/pressure

    def getPhiByMolAm(self, T, x, pressure=None):
        #Compute relative humidity based on mol amount and T
        if pressure is None:
            pressure=self.pressure
        return x*pressure/self.getSatPressure(T)

    def getVaporMassByPhi(self, T, phi, pressure=None):
        #Compute kg_vapor/kg_airdry
        if pressure is None:
            pressure=self.pressure
        if phi > 1: phi/=100
        if phi > 1: phi/=100
        massPerTotalMass = phi*self.getSatPressure(T)/pressure*self.molMassWater/self.molMassAir
        return  massPerTotalMass/(1-massPerTotalMass)

    def getPureProperty(self, T, substance, prop):
        #T in K
        T = T-273.15
        poly = 0
        n = 0
        parameters = self.para[(substance, prop)]
        for parameter in parameters:
            poly+=parameter*T**n
            n+=1
        return poly

    def getPropertyByPhi(self, T, phi, prop, pressure=None):
        # i = dryair, j=watervapor
        # cp is returned in kJ/kg/K
        if pressure is None:
            pressure=self.pressure
        if phi > 1: phi/=100
        mWater = self.molMassWater
        mAirDry = self.molMassAir
        xVapor = self.getMolByPhi(T, phi, pressure)
        xAirDry = 1-xVapor
        propVapor = self.getPureProperty(T, 'vapor', prop)
        propAirDry = self.getPureProperty(T, 'air', prop)
        muVapor = self.getPureProperty(T, 'vapor', 'mu')
        muAir = self.getPureProperty(T, 'vapor', 'mu')
        theta = {}
        theta['air'] = 2**(.5)/4*(mWater/(mWater+mAirDry))**(.5) \
                *(1+(muAir/muVapor)**(.5)*(mWater/mWater)**(.25))**2
        theta['vapor'] = 2**(.5)/4*(mAirDry/(mWater+mAirDry))**(.5) \
                *(1+(muVapor/muAir)**(.5)*(mWater/mAirDry)**(.25))**2
        
        return propAirDry/(1+theta['air']*xVapor/xAirDry) \
                + propVapor/(1+theta['vapor']*xAirDry/xVapor)

    
    def getPropertyByMol(self, T, xVapor, prop):
        # i = dryair, j=watervapor
        # cp is returned in kJ/kg/K
        mWater = self.molMassWater
        mAirDry = self.molMassAir
        xAirDry = 1-xVapor
        propVapor = self.getPureProperty(T, 'vapor', prop)
        propAirDry = self.getPureProperty(T, 'air', prop)
        muVapor = self.getPureProperty(T, 'vapor', 'mu')
        muAir = self.getPureProperty(T, 'vapor', 'mu')
        theta = {}
        theta['air'] = 2**(.5)/4*(mWater/(mWater+mAirDry))**(.5) \
                *(1+(muAir/muVapor)**(.5)*(mWater/mWater)**(.25))**2
        theta['vapor'] = 2**(.5)/4*(mAirDry/(mWater+mAirDry))**(.5) \
                *(1+(muVapor/muAir)**(.5)*(mWater/mAirDry)**(.25))**2
        
        return propAirDry/(1+theta['air']*xVapor/xAirDry) \
                + propVapor/(1+theta['vapor']*xAirDry/xVapor)
    
    def getDensityByPhi(self, T, phi, pressure=None):
        #Computes the Density of moist air based on relative humidity
        #Density in kg/m^3
        if pressure is None:
            pressure=self.pressure
        if phi > 1: phi/=100
        xVapor = self.getMolByPhi(T, phi, pressure)
        xAirDry = 1-xVapor
        molMassMix = xVapor*self.molMassWater + xAirDry*self.molMassAir
        return pressure*molMassMix/self.R/T

    def getDensityByMol(self, T, xVapor, pressure=None):
        #Computes the Density of moist air based on mol_Vapor/mol_total 
        #Density in kg/m^3
        if pressure is None:
            pressure=self.pressure
        if phi > 1: phi/=100
        xAirDry = 1-xVapor
        molMassMix = xVapor*self.molMassWater + xAirDry*self.molMassAir
        return pressure*molMassMix/self.R/T

    def getEnthalpyByMol(self, T, xVapor):
        #Calulates the enthalpy of a mixture, h(T=273.15K)=0
        T = T-273.15
        return T*self.getPropertyByMol(T, xVapor, 'cp')
    
    def getReynolds(self, T, xVapor, velocity, diameter, pressure=None):
        #Computes Reynolds number based on temperatur, mol amount, velocity and diameter
        if pressure is None:
            pressure=self.pressure
        if phi > 1: phi/=100
        return self.getDensityByMol(T, xVapor, pressure)*velocity*diameter \
                / self.getPropertyByMol(T, xVapor, 'mu')

    def getPrandtl(self, T, xVapor):
        #Computes Prandtl number based on temperatur and mol amount
        return self.getPropertyByMol(T,xVapor,'cp')*self.getPropertyByMol(T,xVapor,'mu') \
                / self.getPropertyByMol(T,xVapor,'lambda')

    def getAlpha(self, T, xVapor, velocity, diameter, pressure=None):
        if pressure is None:
            pressure=self.pressure
        if phi > 1: phi/=100
        Nusselt = self.getReynolds(T, xVapor, velocity, diameter, pressure=self.pressure)**.5\
                * self.getPrandtl(T, xVapor)**(1./3.) * 0.331
        return Nusselt*self.getPropertyByMol(T, xVapor,'lambda')/diameter


    

        
