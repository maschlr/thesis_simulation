# Class to calculate sun position
from numpy import pi, sin, cos, tan, array, dot, vdot, sqrt, arccos, degrees, sum

class sunPos(object):
    def __init__(self, latitude, longitude, month, startDay, startTime, surfaceSlope=0, year=2012):
        self.leapyears = range(1900, 2200, 4)
        self.leapyears.remove(1900)
        self.leapyears.remove(2100)
        self.year = year
        self.month = month
        self.startDay = startDay
        self.startTime = startTime
        self.getJulianStartDay()
        self.lastLeapyear = self.findInList(year-1, self.leapyears)
        self.lontitudeTimezone = self.findInList(longitude, range(0,360,15))
        self.hoursToDecLastLeapyear = ((year-self.lastLeapyear-1)*365+6.5)*24
        self.latitude = latitude
        self.longitude = longitude
        self.surfaceSlope = 0
        self.omegaSun = 2*pi/8766.15265
        self.omegaEarth = 2*pi/23.9345
        self.cSin = sin(23.45*pi/180)
        self.cCos = cos(23.45*pi/180)
        #Initialize the vectors for t0
        self.nHorizontal = array([cos(latitude*pi/180), 0, sin(latitude*pi/180)])
        if latitude < 0:
            self.nSloped = array([cos((latitude+surfaceSlope)*pi/180), 0, sin((latitude+surfaceSlope)*pi/180)])
            self.vecNorth = array([-sin(latitude*pi/180), 0, -cos(latitude*pi/180)])
        else:
            self.nSloped = array([cos((latitude+surfaceSlope)*pi/180), 0, sin((latitude-surfaceSlope)*pi/180)])
            self.vecNorth = array([-sin(latitude*pi/180), 0, cos(latitude*pi/180)])
        self.vecSun = array([cos(23.45*pi/180), 0 , -sin(23.45*pi/180)])
        #Vectors for exact sun position, not only angle of inclination
        self.vecEast = array([0,1,0])
    
    def getJulianDay(self, timePassed):
        return self.julianStartDay+(self.startTime+timePassed)/24.

    def getJulianStartDay(self):
        #Method to find out the day, based on Julian calendar, Jan 1st = 1
        daysInMonth = array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]) 
        if self.year in self.leapyears:
            daysInMonth[1] = 29
        self.julianStartDay = 0
        if self.month > 1:
            for i in range(self.month-1):
                self.julianStartDay += daysInMonth[i]
        self.julianStartDay += self.startDay

    def findInList(self, item, itemList):
        for i in range(len(itemList)):
            if item==itemList[i]:
                return itemList[i]
            elif item < itemList[i]:
                return itemList[i-1]
        return False

    def rotMatEarth(self, time):
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
    
    def getRealTime(self, timePassed):
        B = (self.getJulianDay(timePassed)-1)*360/365
        E = 229.2*(7.5e-5 + 1.868e-3*cos(B) - 3.2077e-2*sin(B) \
                - 1.4615e-2*cos(2*B)-4.089e-2*sin(2*B))
        return self.hoursToDecLastLeapyear + (self.julianStartDay-1)*24 + self.startTime \
                + timePassed + (4*(self.lontitudeTimezone-self.longitude)+E)/60

    def getCos(self, vector1, vector2):
        #Computes the Cosinus of the angle beween to vectors
        return vdot(vector1,vector2)/sqrt(sum(vector1**2))/sqrt(sum(vector2**2))

    def getSlopeFactor(self, timePassed):
        time = self.getRealTime(timePassed)
        vecST = dot(self.rotMatSun(time), self.vecSun)
        nHT = dot(self.rotMatEarth(time), self.nHorizontal)
        nSlopedT = dot(self.rotMatEarth(time), self.nSloped)
        return self.getCos(vecST,nSlopedT)/self.getCos(vecST,nHT)

    def getSunPosition(self, timePassed):
        time = self.getRealTime(timePassed)
        vecST = dot(self.rotMatSun(time), self.vecSun)
        nHT = dot(self.rotMatEarth(time), self.nHorizontal)
        vecNT = dot(self.rotMatEarth(time), self.vecNorth)
        vecET = dot(self.rotMatEarth(time), self.vecEast)
        incl = self.getCos(vecST, nHT)
        projection = vecST - vdot(vecST,nHT)/vdot(nHT, nHT) * nHT
        cosNorth = self.getCos(projection, vecNT)
        cosEast = self.getCos(projection, vecET)

        if cosEast>0: zeta = degrees(arccos(cosNorth))
        else: zeta = 360-degrees(arccos(cosNorth))
        declination = 90 - degrees(arccos(incl))

        return [declination, zeta]

