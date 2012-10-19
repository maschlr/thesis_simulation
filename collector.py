#Class to calculate the beloved collector with 3 paths and meander cooling fins like in the publication for Lima
from numpy import array, floor, pi, empty, arange, average, exp, sqrt, tanh, ones, alen, zeros
from weatherdata import weatherData
from moistAir import AirMix
from sunPosition import sunPos
from scipy import interpolate as ip
from scikits.odes.sundials import ida
from scikits.odes import dae
import matplotlib.pyplot as plt

#Set following False to not compute solution with ddaspk
alsoddaspk = True
ATOL = 1e-5
RTOL = 1e-3
JACFAC = 1e-1

class Collector():
    def __init__(self, month, startDay, startTime, timeSpan, data=None):
        if data is None:
            self.v_a = .5
            self.L = 4.
            self.b = 1.
            self.d = 0.05
            self.A = self.b*self.L
            self.d_h = 2.*self.b*self.d/(self.b+self.d)
    
            #Thicknesses
            self.l_s = 0.1
            self.l_i = 0.1
            self.l_p = 0.007
            self.l_g = 0.005	

            #Design parameters for the fins
            self.b_2 = 0.75
            self.d_2 = self.d
            self.d_h2 = 2*(self.b-self.b_2)*self.d_2/((self.b-self.b_2)+self.d_2)
            self.n_f = floor(self.L/(self.b-self.b_2)) - 1
            self.L_f = self.b*(self.n_f+1)
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
            self.endDay = self.startDay + floor(timeSpan/24.)
        else:
            self.endDay = None
        #Access to sun position functions
        self.sp=sunPos(self.latitude, self.longitude, self.month, self.startDay, self.startTime, surfaceSlope=self.beta)

        #Access to the database with weatherdata
        self.recData = weatherData()
        self.recData.readDB(self.month, self.startDay, self.startTime, self.timeSpan)

        #columns: 0:day, 1:month, 2:year, 3:hours, 4:minutes, 5:seconds from beginning of month
        # 6:temperature, 7:humidity, 8:windspeed, 9:insolation
        self.relevantData = array(self.recData.relevantData)
        self.timeline = arange(0., len(self.relevantData[:,1])*1800., 1800.)
        
        #Calculate the spline interpolations for the relevant data
        self.rTemps = ip.splrep(self.timeline, self.relevantData[:,6]+273.15, s=0)
        self.rInsolation = ip.splrep(self.timeline, self.relevantData[:,9], s=0)
        self.rWindspeed = ip.splrep(self.timeline, self.relevantData[:,8], s=0)

        self.T0 = self.relevantData[0,6]+273.15
        
        #calculate some constants that are considered constant during simulation time
        self.am = AirMix()
        self.hc, self.hr  = empty([2,3])
        averageTemp = average(self.relevantData[:,6])+273.15
        averageMol = self.am.getMolByPhi(averageTemp,.5)
        self.hc[0] = self.am.getAlpha(averageTemp, averageMol, self.v_a, self.d_h)
        self.hc[1] = self.am.getAlpha(averageTemp, averageMol, self.v_a2, self.d_h2)
        self.hc[2] = self.hc[0]
        self.hr[0] = self.sigma*(2.+1./self.epsilon_g1)**(-1)
        self.hr[1] = self.sigma*(1.+1./self.epsilon_g1+1./self.epsilon_g2)**(-1)
        self.hr[2] = self.sigma*(1.+1./self.epsilon_g2+1./self.epsilon_p)**(-1)
        self.m_a = self.v_a*self.b*self.d*self.am.getDensityByMol(averageTemp, averageMol)
        self.c_pa = self.am.getPropertyByMol(averageTemp, averageMol,'cp')
        self.mcp = self.m_a * self.c_pa
        m_f = sqrt(self.hc[1]*(self.l_f+2.*self.b_2)/self.lambda_p/self.l_f/self.b_2)
        self.eta_f = tanh(m_f*self.d_2)/m_f/self.d_2

        #Variables for solver
        self.stop_t = arange(.0,self.timeline[-1],300.)
        self.algvar_idx=[3,4,5,6,7,8]
        self.exclalg_err = False
        self.neq = 9
        self.y0 = ones(9)*self.T0
        self.yp0 = ones(9)*-1

        self.solution = {}

    def set_res(self, resfunction):
        """Function to set the residual function as required by IDA
           needed for the ddaspk simulation"""
        self.res = resfunction

    def ddaspk_res(self, tres, yy, yp, res):
        """the residual function as required by ddaspk"""
        self.res.evaluate(tres, yy, yp, res, None)
    
    def setSolution(self, y, solver):
        self.solution[solver] = y
        self.plt_t = ones((len(self.stop_t),6))
        for i in range(6):
            self.plt_t[:,i]*=self.stop_t

class resindex(ida.IDA_RhsFunction):
    """ Residual function class as needed by the IDA DAE solver """

    def set_drysim(self, collectorData):
        """ Set the drying simulation to solve to have access to the data """
        self.cd = collectorData

    def evaluate(self, t, x, xdot, res, userdata):
        cd = self.cd
        ambTemp = ip.splev(t, cd.rTemps)
        h_cw = 2.8+3.*abs(ip.splev(t, cd.rWindspeed))
        insol = 1.2*ip.splev(t, cd.rInsolation)
        u_b = (1./h_cw + cd.l_i/cd.lambda_i + 1./cd.hc[2])**(-1)

        res[0] = h_cw*(ambTemp-x[0]) + cd.hc[0]*(x[6]-x[0]) + cd.hr[1]*(x[1]**4-x[0]**4) \
                + cd.hr[0]*((ambTemp-6)**4-x[0]**4) + insol*cd.alpha_g1 \
                - xdot[0]*cd.l_g*cd.rho_g*cd.c_g
        res[1] = cd.hc[1]*(x[7]-x[1]) + cd.hc[0]*(x[6]-x[1]) + cd.hr[1]*(x[0]**4-x[1]**4) \
                + cd.hr[2]*(x[2]**4-x[1]**4) + insol*cd.tau_g1*cd.alpha_g2 \
                - xdot[1]*cd.l_g*cd.rho_g*cd.c_g
        res[2] = cd.hc[2]*cd.A*(x[8]-x[2]) \
                + cd.hc[1]*(cd.A+cd.eta_f*cd.b_2*cd.d_2*cd.n_f)*(x[7]-x[2]) \
                + cd.hr[2]*cd.A*(x[1]**4-x[2]**4) + insol*cd.tau_g1*cd.tau_g2*cd.alpha_p*cd.A\
                - xdot[2]*(cd.A*cd.l_s*cd.rho_s*cd.c_s+(cd.A+cd.n_f*cd.d_2*cd.b_2)*cd.c_p*cd.rho_p)
        res[3] = (x[0]+x[1])/2. + (ambTemp-(x[0]+x[1])/2.)*exp(-2.*cd.hc[0]*cd.b*cd.L/cd.mcp)\
                -x[3]
        res[4] = ((cd.b-cd.b_2+cd.eta_f*cd.d_2*2)*x[2] + (cd.b-cd.b_2)*x[1]) \
                / (2*cd.b-2*cd.b_2+2*cd.eta_f*cd.d_2) \
                + (x[3] - ((cd.b-cd.b_2+2*cd.eta_f*cd.d_2)*x[2] + (cd.b-cd.b_2)*x[1]) \
                / (2*cd.b-2*cd.b_2+2*cd.eta_f*cd.d_2)) \
                * exp(-cd.hc[1]*(2*cd.b-2*cd.b_2+2*cd.eta_f*cd.d_2)/cd.mcp*cd.L_f) \
                -x[4]
        res[5] = (cd.hc[2]*x[2] + u_b*ambTemp)/(cd.hc[2]+u_b) \
                + (x[4] - (cd.hc[2]*x[2] + u_b*ambTemp)/(cd.hc[2]+u_b)) \
                * exp(-(cd.hc[2]+u_b)*cd.b*cd.L/cd.mcp) \
                - x[5]
        res[6] = (x[0]+x[1])/2. + (ambTemp-(x[0]+x[1])/2.) * (0.5*cd.mcp/cd.hc[0]/cd.b/cd.L) \
                * (1.-exp(-2.*cd.hc[0]*cd.b*cd.L/cd.mcp)) \
                -x[6]
        res[7] = ((cd.b-cd.b_2+2*cd.eta_f*cd.d_2)*x[2] + (cd.b-cd.b_2)*x[1]) \
                / (2*cd.b - 2*cd.b_2 + 2*cd.eta_f*cd.d_2) \
                + (x[3] - ((cd.b-cd.b_2 + 2*cd.eta_f*cd.d_2)*x[2] + (cd.b-cd.b_2)*x[1]) \
                / (2*cd.b-2*cd.b_2 + 2*cd.eta_f*cd.d_2)) \
                * (cd.mcp/cd.hc[1]/(2*cd.b-2*cd.b_2+2*cd.eta_f*cd.d_2)/cd.L_f) \
                * (1.-exp(-cd.hc[1]*(2*cd.b-2*cd.b_2+2*cd.eta_f*cd.d_2)/cd.mcp*cd.L_f)) \
                -x[7]
        res[8] = (cd.hc[2]*x[2] + u_b*ambTemp)/(cd.hc[2]+u_b) \
                + (x[4] - (cd.hc[2]*x[2] + u_b*ambTemp)/(cd.hc[2]+u_b)) \
                * (cd.mcp/(cd.hc[2]+u_b)/cd.b/cd.L) * (1.-exp(-(cd.hc[2]+u_b)*cd.b*cd.L/cd.mcp)) \
                -x[8]
        return 0

def calculateSolution():
    jac = None
    prob = Collector(month=2, startDay=4, startTime=5, timeSpan=24*10)
    res = resindex()
    
    res.set_drysim(prob)
    jfac = 1.

    solver = ida.IDA(res,
                compute_initcond='yp0',
                first_step_size=1e-18,
                atol=ATOL*jfac,
                rtol=RTOL*jfac,
                max_steps=1500,
                jacfn=jac,
                algebraic_vars_idx=prob.algvar_idx,
                exclude_algvar_from_error=prob.exclalg_err,
                )

    
    
    # Solution vectors
    y = [0]*(1+len(prob.stop_t)); yp = [0]*(1+len(prob.stop_t))
    y[0] = empty(prob.neq, float); yp[0] = empty(prob.neq, float)
    flag, t, y = solver.solve(prob.stop_t, prob.y0, prob.yp0)[:3]
    prob.setSolution(y, 'ida')
#    (flag, t0_init) = solver.init_step(0., prob.y0, prob.yp0, y[0], yp[0])
#    realtime = [t0_init]
#
#    i=1
#    error = False
#    for time in prob.stop_t[1:]:
#            #print 'at time', time
#            y[i] = empty(prob.neq, float)
#            yp[i] = empty(prob.neq, float)
#            flag, rt = solver.step(time, y[i], yp[i])
#            realtime += [rt]
#            print(time, y[i])
#
#            i += 1
#            if flag != 0:
#                error = True
#                print('Error in solver, breaking solution at time %g' % time)
#                break
    
    if alsoddaspk:
        ddaspkz = empty((alen(prob.stop_t), prob.neq), float)
        ddaspkzprime = empty((alen(prob.stop_t), prob.neq), float)

        prob.set_res(res)
        ig = dae('ddaspk', prob.ddaspk_res)
        if jac:
            prob.set_jac(jac)
            ig.set_options(jacfn=prob.ddaspk_jac)
        #first compute the correct initial condition from the values of z0
        ig.set_options(
                        algebraic_vars_idx=prob.algvar_idx,
                        compute_initcond='yp0',
                        first_step=1e-18,
                        exclude_algvar_from_error=prob.exclalg_err,
                        atol=ATOL, rtol=RTOL, max_steps=1500)
        tinit = ig.init_step(0., prob.y0, prob.yp0, ddaspkz[0], ddaspkzprime[0])

#        print('ddaspk started from y0 = ', prob.y0)
#        print('ddaspk initial condition calculated, [z,zprime] = [', ddaspkz[0],
#                    ddaspkzprime[0], ']')

        i=1
        error = False
        for time in prob.stop_t[1:]:
            flag, tout = ig.step(time, ddaspkz[i],  ddaspkzprime[i])
            i += 1

            if flag < 1:
                error = True
                break
        
        prob.setSolution(ddaspkz, 'ddaspk')
    return prob

