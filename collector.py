#Class to calculate the beloved collector with 3 paths and meander cooling fins like in the publication for Lima
from numpy import array, floor, pi, empty
from weatherdata import weatherData
from moistAir import AirMix as amix

class Collector():
    def __init__(self, month, startDay, startTime, timeSpan, data=None):
        if data is None:
            self.v_a = 1e-10;
            self.L = 5;
            self.b = 2;
            self.d = 0.05;
    
            #Thicknesses
            self.l_s = 0.001;
            self.l_i = 0.1;
            self.l_p = 0.007;
            self.l_g = 0.005;	

            #Design parameters for the fins
            self.b_2 = 0.5;
            self.d_2 = self.d;
            self.n_f = floor(self.L/(self.b-self.b_2)) - 1;
            self.L_f = self.L + (self.n_f-1) *self.b_2;
            self.l_f = self.l_p; 
    
    
            self.lambda_i = 0.043;
            self.lambda_p = 137.0;
            self.alpha_g1 = 0.1;
            self.alpha_g2 = 0.1;
            self.alpha_p = 0.9;
            self.epsilon_g1 = 0.8;
            self.epsilon_g2 = 0.8;
            self.epsilon_p = 0.95;
            self.tau_g1 = 0.9; 
            self.tau_g2 = 0.8;
            self.rho_s = 2700.0;
            self.c_s = 790.0;
            self.c_p = 452.0;
            self.c_g = 800.0;	
            self.rho_g = 2480;
            self.rho_p = 7860;
 
            self.beta = -18*pi/180;
    
            #positive for northern hemisphere, negative for southern
            self.phi = -5.2*pi/180;
        
        self.month = month 
        self.startDay = startDay
        # startTime and timeSpan in hours from 0:00h
        self.startTime = startTime
        self.timeSpan = timeSpan
        if timeSpan>24:
            self.endDay = self.startDay + floor(timeSpan/24)
        else:
            self.endDay = None

        #Access to the database with weatherdata
        self.recData = weatherData()
        self.recData.readDB(self.month, self.startDay, self.startTime, self.timeSpan)

        #columns: 0:day, 1:month, 2:year, 3:hours, 4:minutes, 5:seconds from beginning of month
        # 6:temperature, 7:humidity, 8:windspeed, 9:insolation
        self.relevantData = array(self.recData.relevantData)

        #calculate some constants that are considered constant during siulation time
        self.am = amix()
        self.h_c, self.h_r  = empty([2,3])
        self.h_c[0] = self.am.getAlpha(
        
class resindex(ida.IDA_RhsFunction):
    """ Residual function class as needed by the IDA DAE solver """

    def set_drysim(self, drysim):
        """ Set the drying simulation to solve to have access to the data """
        self.dsim = drysim

    def evaluate(t, x, xdot, res, userdata):

        

if __name__=='__main__':
    collector = Collector(4, 5, 8, 32)
    print collector.recData.relevantData
