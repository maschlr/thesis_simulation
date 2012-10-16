# Class to calculate sun position
from numpy import pi, sin, cos, tan, array, dot, vdot, sqrt, acosd

class sunPos(object):
    def __init__(self, latitude, longitude, julianStartDay, surfaceSlope=0, year=2012):
        self.leapyears = range(1900, 2200, 4)
        self.leapyears.remove(1900)
        self.leapyears.remove(2100)
        self.year = year
        self.lastLeapyear = findInList(year, self.leapyears)
        self.lontitudeTimezone = findInList(longitude, range(0,360,15))
        self.julianStartDay = julianStartDay
        self.hoursToDecLastLeapyear = (year-self.lastLeapyear)*8760+7

        self.latitude = latitude
        self.longitude = longitude
        self.surfaceSlope = 0
        self.omegaSun = 2*pi/8766.15265
        self.omegaEarth = 2*pi/23.9345
        self.cSin = sin(23.45*pi/180)
        self.cCos = cos(23.45*pi/180)

        #Initialize the vectors for t0
        self.nHorizontal = array([1, 0, tan(latitude*pi/180)])
        if latitude < 0:
            self.nSloped = array([1, 0, tan((latitude+surfaceSlope)*pi/180)])
        else:
            self.nSloped = array([1, 0, tan((latitude-surfaceSlope)*pi/180)])
        self.vecSun = array([1, 0 , tan(-23.45*pi/180)])

        #Vectors for exact sun position, not only angle of inclination
        self.vecNorth = array([-tan(latitude*pi/180), 0, 1])
        self.vecEast = array([0,1,0])
    
    def getJulianDay(self, locationTime):
        return self.startJulianDay+time/24

    def findInList(self, item, itemList):
        for lItem in itemList:
            if item>=lItem:
                return item
        return False

    def rotMatEart(self, time):
        #Returns the earth rotation matrix at a given time difference
        omega = self.omegaEarth
        return array([[cos(omega*time), -sin(omega*time), 0], \
                [sin(omega*time), cos(omega*time), 0], \
                [0, 0, 1]])

    def rotMatSun(self, time):
        #Returns the sun rotation matrix at a given time difference
        cSin = self.cSin
        cCos = self.cCos
        vSin = sin(self.omegaSun*time)
        vCos = cos(self.omegaSun*time)
        
        return array([[cSin**2*(1-vCos)+vCos, -cCos*vSin, cSin*cCos*(1-vCos)], \
                [cCos*vSin, vCos, -cSin*vSin], \
                [cCos*cSin*(1-vCos), cSin*vSin, cCos**2*(1-vCos)+vCos]])
    
    def getRealTime(self, locationTime):
        B = (self.getJulianDay(localTime)-1)*360/365
        E = 229.2*(7.5e-5 + 1.868e-3*cos(B) - 3.2077e-2*sin(B) \
                - 1.4615e-2*cos(2*B)-4.089e-2*sin(2*B))
        return locationTime + 4*(self.lontitudeTimezone-self.longitude)+E

    def getCos(self, vector1, vector2):
        #Computes the Cosinus of the angle beween to vectors
        return abs(vdot(vector1, vector2))/sqrt(vdot(vector2,vector2))/sqrt(vdot(vector1, vector1))

    def getSlopeFactor(self, time):
        vecST = dot(self.rotMatSun(time), self.vecSun)
        nHT = dot(self.rotMatEarth(time), self.nHorizontal)
        nSlopedT = dot(self.rotMatEarth(time), self.nSloped)
        return getCos(vecST, nSlopedT)/getCos(vecST, nHT)

    def getSunPosition(self, time):
        vecST = dot(self.rotMatSun(time), self.vecSun)
        nHT = dot(self.rotMatEarth(time), self.nHorizontal)
        vecNT = dot(self.rotMatEarth(time), self.vecNorth)
        vecrET dot(self.rotMatEarth(time), self.vecEast)
        incl = getCos(vecST, nHT)
        projection = vecST - vdot(vecST,nHT)/vdot(nHT, nHT) * nHT
        cosNorth = getCos(projection, vecNT)
        cosEast = getCos(projection, vecET)

        if cosEast>0: zeta = acosd(cosNorth)
        else: zeta = 360-acosd(cosNorth)
        declination = 90 - acosd(incl)

        return [declination zeta]




        

