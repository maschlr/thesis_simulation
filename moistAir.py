#Class for calculation the thermopysical Properties of moist air
from numpy import array, arange, empty

class AirMix(object):
    def __init__(self, pressure=1.01325e5):
        #These parameters are valid around atmospheric pressure
        self.para = {}
        self.para[('air','cp')] = array([1.03409, 2.84887e-4, 4.97078e-10, 1.07702e-13])
        self.para[('air','mu')] = array([1.43387e-6, 6.56244e-8, 2.9905e-11])
        self.para[('air','lambda')] = array([6.69881e-4, 9.42482e-5, 3.2745e-8])
        self.para[('vapor','cp')] = array([3.387781, -5.55984e-3, 1.94105e-5])
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

    def getMolByPhi(self, T, phi, pressure=self.pressure):
        #Calculates mol_water per mol_airmix based on T in K and RH
        #relative humidity should be smaller than 1, if given in percent
        if phi > 1: phi/=100
        return phi*getSatPressure(T)/pressure

    def getPhiByMolAm(self, T, x, pressure=self.pressure):
        #Compute relative humidity based on mol amount and T
        return x*pressure/getSatPressure(T)

    def getVaporMassByPhi(self, T, phi, pressure=self.pressure):
        #Compute kg_vapor/kg_airdry
        if phi > 1: phi/=100
        massPerTotalMass = phi*getSatPressure(T)/pressure*self.molMassWater/self.molMassAir
        return  massPerTotalMass/(1-massPerTotalMass)

    def getPureProperty(self, T, substance, prop):
        #Computes heat capcity of dry air, T in K
        T = T-273.15
        poly = 0
        n = 0
        parameters = self.para[(substance, prop)]
        for parameter in parameters:
            poly+=parameter*T**n
        return poly

    def getMixPropertyByPhi(self, T, phi, prop, pressure=self.pressure):
        # i = dryair, j=watervapor
        mWater = self.molMassWater
        mAirDry = self.molMassAir
        xVapor = getMolByPhi(T, phi, pressure)
        xAirDry = 1-xVapor
        propVapor = getPureProperty(T, 'vapor', prop)
        propAirDry = getPurePropery(T, 'air', prop)
        muVapor = getPureProperty(T, 'vapor', 'mu')
        muAir = getPureProperty(T, 'vapor', 'mu')
        theta = {}
        theta['air'] = 2**(.5)/4*(mWater/(mWater+mAirDry))**(.5) \
                *(1+(muAir/muVapor)**(.5)*(mWater/mVapor)**(.25))**2
        theta['vapor'] = 2**(.5)/4*(mAirDry/(mWater+mAirDry))**(.5) \
                *(1+(muVapor/muAir)**(.5)*(mVapor/mAirDry)**(.25))**2
        
        return probAirDry/(1+theta['air']*xVapor/xAirDry) \
                + probVapor/(1+theta['vapor']*xAirDry/xVapor)

    
    def getMixPropertyByMol(self, T, xVapor, prop):
        # i = dryair, j=watervapor
        mWater = self.molMassWater
        mAirDry = self.molMassAir
        xAirDry = 1-xVapor
        propVapor = getPureProperty(T, 'vapor', prop)
        propAirDry = getPurePropery(T, 'air', prop)
        muVapor = getPureProperty(T, 'vapor', 'mu')
        muAir = getPureProperty(T, 'vapor', 'mu')
        theta = {}
        theta['air'] = 2**(.5)/4*(mWater/(mWater+mAirDry))**(.5) \
                *(1+(muAir/muVapor)**(.5)*(mWater/mVapor)**(.25))**2
        theta['vapor'] = 2**(.5)/4*(mAirDry/(mWater+mAirDry))**(.5) \
                *(1+(muVapor/muAir)**(.5)*(mVapor/mAirDry)**(.25))**2
        
        return probAirDry/(1+theta['air']*xVapor/xAirDry) \
                + probVapor/(1+theta['vapor']*xAirDry/xVapor)
    
    def getMixDensityByPhi(self, T, phi, pressure=self.pressure):
        #Computes the Density of moist air based on relative humidity
        #Density in kg/m^3
        xVapor = getMolByPhi(T, phi, pressure)
        xAirDry = 1-xVapor
        molMassMix = xVapor*self.molMassWater + xAirDry*self.molMassAir
        return pressure*molMassMix/self.R/T

    def getMixDensityByMol(self, T, xVapor, pressure=self.pressure):
        #Computes the Density of moist air based on mol_Vapor/mol_total 
        #Density in kg/m^3
        xAirDry = 1-xVapor
        molMassMix = xVapor*self.molMassWater + xAirDry*self.molMassAir
        return pressure*molMassMix/self.R/T

    def getEnthalpyByMol(self, T, xVapor):
        #Calulates the enthalpy of a mixture, h(T=273.15K)=0
        T = T-273.15
        return T*getMixPropertyByMol(T, xVapor, 'cp')
    
    def getReynolds(self, T, xVapor, velocity, diameter, pressure=self.pressure):
        #Computes Reynolds number based on temperatur, mol amount, velocity and diameter
        return getMixDensityByMol(T, xVapor, pressure)*velocity*diameter \
                / getMixPropertyByMol(T, xVapor, 'mu')

    def getPrandtl(self, T, xVapor):
        #Computes Prandtl number based on temperatur and mol amount
        return getMixPropertyByMol(T,xVapor,'cp')*getMixPropertyByMol(T,xVapor,'mu') \
                / getMixProperyByMol(T,xVapor,'lambda')

    def getAlpha(self, T, xVapor, velocity, diameter, pressure=self.pressure):
        Nusselt = self.getReynolds(T, xVapor, velocity, diameter, pressure=self.pressure)**.5\
                * self.getPrandtl(T, xVapor)**(1./3.) * 0.331

        return Nusselt*self.getMixPropertyByMol(T, xVapor,'lambda')/diameter


    

        
