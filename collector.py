#Class to calculate the beloved collector with 3 paths and meander cooling fins like in the publication for Lima
from numpy import array, floor, pi, empty, arange, average, exp, sqrt, tanh
from weatherdata import weatherData
from moistAir import AirMix
from sunPosition import sunPos
from scipy import interpolate as ip
from scikits.odes.sundials import ida
from scikits.odes import dae

class Collector():
    def __init__(self, month, startDay, startTime, timeSpan, data=None):
        if data is None:
            self.v_a = .3
            self.L = 5.
            self.b = 2.
            self.d = 0.05
            self.A = self.b*self.L
            self.d_h = 2*self.b*self.d/(self.b+self.d)
    
            #Thicknesses
            self.l_s = 0.001
            self.l_i = 0.1
            self.l_p = 0.007
            self.l_g = 0.005	

            #Design parameters for the fins
            self.b_2 = 0.5
            self.d_2 = self.d
            self.d_h2 = 2*(self.b-self.b_2)*self.d_2/((self.b-self.b_2)+self.d_2)
            self.n_f = floor(self.L/(self.b-self.b_2)) - 1
            self.L_f = self.L + (self.n_f-1) * self.b_2
            self.v_a2 = self.v_a*self.b*self.d/(self.b-self.b_2)/self.d_2
            self.l_f = self.l_p
    
            #These are material constants
            self.lambda_i = 0.043
            self.lambda_p = 137.0
            self.alpha_g1 = 0.1
            self.alpha_g2 = 0.1
            self.alpha_p = 0.9
            self.sigma = 5.67e-8 
            self.epsilon_g1 = 0.8
            self.epsilon_g2 = 0.8
            self.epsilon_p = 0.95
            self.tau_g1 = 0.9
            self.tau_g2 = 0.8
            self.rho_s = 2700.0
            self.c_s = 790.0
            self.c_p = 452.0
            self.c_g = 800.0	
            self.rho_g = 2480.
            self.rho_p = 7860.
 
            #Slope angle of the collector
            self.beta = 5
            self.latitude = -5.2
            self.longitude = 80.633333
        
        self.month = month 
        self.startDay = startDay
        # startTime and timeSpan in hours from 0:00h
        self.startTime = startTime
        self.timeSpan = timeSpan
        if timeSpan>24:
            self.endDay = self.startDay + floor(timeSpan/24)
        else:
            self.endDay = None
        #Access to sun position functions
        self.sp=SunPos(self.latitude, self.longitude, self.month, self.startDay, self.startTime, surfaceSlope=self.beta)

        #Access to the database with weatherdata
        self.recData = weatherData()
        self.recData.readDB(self.month, self.startDay, self.startTime, self.timeSpan)

        #columns: 0:day, 1:month, 2:year, 3:hours, 4:minutes, 5:seconds from beginning of month
        # 6:temperature, 7:humidity, 8:windspeed, 9:insolation
        self.relevantData = array(self.recData.relevantData)
        self.timeline = arange(0, len(self.relevantData[:,1])*1800, 1800)
        
        #Calculate the spline interpolations for the relevant data
        self.rTemps = ip.splrep(self.timeline, self.relevantData[:,6]+273.15, s=0)
        self.rInsolation = ip.splrep(self.timeline, self.relevantData[:,9], s=0)
        self.rWindspeed = ip.splrep(self.timeline, self.relevantData[:,8], s=0)

        self.T0 = self.relevantData[0,6]
        
        #calculate some constants that are considered constant during simulation time
        self.am = AirMix()
        self.h_c, self.h_r, self.facDiff  = empty([3,3])
        averageTemp = average(self.relevantData[:,6])+273.15
        averageMol = self.am.getMolByPhi(averageTemp,.5)
        self.h_c[0] = self.am.getAlpha(averageTemp, averageMol, self.v_a, self.d_h)
        self.h_c[1] = self.am.getAlpha(averageTemp, averageMol, self.v_a2, self.d_h2)
        self.h_c[2] = self.h_c[0]
        self.h_r[0] = self.sigma*(2+1/self.epsilon_g1)**(-1)
        self.h_r[1] = self.sigma*(1+1/self.epsilon_g1+1/epsilon_g2)**(-1)
        self.h_r[2] = self.sigma*(1+1/self.epsolin_g2+1/epsilon_p)**(-1)
        self.facDiff[0] = self.l_g*self.c_g
        self.facDiff[1] = self.facDiff[0]
        self.facDiff[2] = self.A*self.l_s*self.c_s+(self.A+self.n_f*(self.b_2-self.b))*self.l_p*self.c_p
        self.m_a = self.v_a*self.b*self.d*self.am.getDensityByMol(averageTemp, averageMol)
        self.c_pa = self.am.getPropertyByMol(averageTemp, averageMol)
        self.mcp = self.m_a*self.c_pa
        m_f = sqrt(self.h_c[1]*(self.l_f+2*self.b_2)/self.lambda_p/self.l_f/self.b_2)
        self.eta_f = tanh(m_f*d_2)/m_f/d_2
        
class resindex(ida.IDA_RhsFunction):
    """ Residual function class as needed by the IDA DAE solver """

    def set_drysim(self, collectorData):
        """ Set the drying simulation to solve to have access to the data """
        self.cd = collectorData

    def evaluate(t, x, xdot, res, userdata):
        cd = self.cd
        ambTemp = ip.splev(t, cd.rTemps)
        h_cw = 2.8+3*ip.splev(t, cd.rWindspeed)
        insol = 1.2*ip.splev(t, cd.rInsolation)
        c_eq4 = 2*(cd.b-cd.b_2+cd.eta_f*cd.d_2)
        u_b = (1/h_cw + cd.l_i/cd.lambda_i + 1/cd.h_c[2])

        res[0] = h_cw*(ambTemp-x[0]) + cd.h_c[0]*(x[6]-x[0]) + cd.h_r[1]*(x[1]**4-x[0]**4) \
                + cd.h_r[0]*((ambTemp-6)**4-x[0]**4) - xdot[0]*cd.facDiff[0] \
                + insol*cd.alpha_g1
        res[1] = cd.h_c[1]*(x[7]-x[1]) + cd.h_c[0]*(x[6]-x[2]) + cd.h_r[1]*(x[0]**4-x[1]**4) \
                + cd.h_r[2]*(x[2]**4-x[1]**4) - xdot[1]*cd.facDiff[1] \
                + insol*cd.tau_g1*cd.alpha_g1
        res[2] = cd.h_c[1]*cd.A*(x[7]-x[2]) + cd.h_c[2]*cd.A*(x[8]-x[2]) \
                + cd.h_r[2]*cd.A*(x[1]**4-x[2]**4) - xdot[2]*cd.facDiff[2] \
                + insol*cd.tau_g1*cd.tau_g2*cd.alpha_p
        res[3] = (x[0]+x[1])/2 + (ambTemp - (x[0]+x[1])/2)*exp(-2*cd.h_c[0]*cd.b*cd.L/cd.mcp)
        res[4] = (cd.h_c[1]*(c_eq4-(cd.b-cd.b_2))*x[2] + (cd.b-cd.b_2)*x[1])/c_eq4 \
                + (x[3] - (cd.h_c[1]*(c_eq4-(cd.b-cd.b_2))*x[2] + (cd.b-cd.b_2)*x[1])/c_eq4) \
                * exp(-cd.h_c[1]*c_eq4/cd.mcp*cd.L_f)
        res[5] = (cd.h_c[2]*x[2] + u_b*ambTemp)/(cd.h_c[2]+u_b) \
                + (x[4] - (cd.h_c[2]*x[2] + u_b*ambTemp)/(cd.h_c[2]+u_b)) \
                * exp(-(cd.h_c[2]+u_b)*cd.b*cd.L/cd.mcp)
        res[6] = (x[0]+x[1])/2 + (ambTemp - (x[0]+x[1])/2) * (0.5*cd.mcp/cd.h_c[0]/cd.b/cd.L) \
                * (1-exp(-2*cd.h_c[0]*cd.b*cd.L/cd.mcp))
        res[7] = (cd.h_c[1]*(c_eq4-(cd.b-cd.b_2))*x[2] + (cd.b-cd.b_2)*x[1])/c_eq4 \
                + (x[3] - (cd.h_c[1]*(c_eq4-(cd.b-cd.b_2))*x[2] + (cd.b-cd.b_2)*x[1])/c_eq4) \
                * (cd.mcp/cd.h_c[1]/c_eq4/cd.L_f) * (1-exp(-cd.h_c[1]*c_eq4/cd.mcp*cd.L_f))
        res[8] = (cd.h_c[2]*x[2] + u_b*ambTemp)/(cd.h_c[2]+u_b) \
                + (x[4] - (cd.h_c[2]*x[2] + u_b*ambTemp)/(cd.h_c[2]+u_b)) \
                * (cd.mcp/(cd.h_c[2]+u_b)/cd.b/cd.L) * (1-exp(-(cd.h_c[2]+u_b)*cd.b*cd.L/cd.mcp))



        

if __name__=='__main__':
    collector = Collector(4, 5, 8, 32)
    print collector.recData.relevantData
