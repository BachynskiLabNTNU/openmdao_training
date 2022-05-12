## Modified from Code by Erin Bachynski 

import openmdao.api as om
import numpy as np

class computeMball(om.ExplicitComponent):
    def setup(self):
         # global design variables
        self.add_input('D')
        self.add_input('T')
        
        # problem constants
        self.add_input('rho',val=1025.0)
        self.add_input('rho_steel',val=7850.0)
        self.add_input('mturb',val=1244000.0)
        self.add_input('hfb',val=10.0)
        
        self.add_output('t')
        self.add_output('mball')
        self.add_output('msteel')
        
        
    def setup_partials(self):
        # Finite difference all partials.
        self.declare_partials('*', '*', method='fd')
        
    def compute(self,inputs,outputs): 
        rho = inputs['rho']
        rho_steel = inputs['rho_steel']
        mturb = inputs['mturb']
        hfb = inputs['hfb']
        
        D = inputs['D']
        T = inputs['T']
        
        # based on T, find thickness and ballast in order to meet vertical
        # equilibrium
        t = 0.03+T/3000;
        buoy = np.pi*np.power(D,2)/4*T*rho # displacement water mass
        msteel = rho_steel*(np.pi*D*(T+hfb)*t+2*np.pi*np.power(D,2)/4*t) # steel mass
        mball = buoy-msteel-mturb # required ballast

        outputs['mball'] = mball
        outputs['msteel'] = msteel


class computeTheta(om.ExplicitComponent): 
    #T,rho,rho_steel,rho_conc, mturb,zturb,FT,hhub,hfb
    def setup(self): 
        # global design variables
        self.add_input('D')
        self.add_input('T')
        
        # result from another component
        self.add_input('mball',val=0.0)
        self.add_input('msteel',val=0.0)
        
        # problem constants
        self.add_input('rho',val=1025.0)
        self.add_input('rho_conc', val=2650.0)
        self.add_input('mturb',val=1244000.0)
        self.add_input('zturb',val=92.5)
        self.add_input('FT',val= 1500000.0)
        self.add_input('hhub',val=120)
        self.add_input('hfb',val=10.0)
        

        # outputs
        self.add_output('theta',val = 0.0)

    
    def setup_partials(self):
        # Finite difference all partials.
        self.declare_partials('*', '*', method='fd')
        
    def compute(self, inputs, outputs):
        """
        
        """
        rho = inputs['rho']
        rho_conc = inputs['rho_conc']
        mturb = inputs['mturb']
        zturb = inputs['zturb']
        FT = inputs['FT']
        hhub = inputs['hhub']
        hfb = inputs['hfb']
        mball = inputs['mball']
        msteel = inputs['msteel']
        
        
        D = inputs['D']
        T = inputs['T']
        
        
        g = 9.81 # acceleration due to gravity

        # we can't have negative ballast, the function returns zero and 90 degree
        # theta
        if mball<0: 
            mball = 0
            theta = np.pi/2
            hbal = 0
        else: 
            #given the ballast mass, we can find the ballast height
            hbal = mball/(rho_conc*np.pi*np.power(D,2)/4)
            
            #total mass of the FWT
            mtot = msteel + mball  + mturb

            #hull steel contribution to center of gravity: m*z 
            mZGsteel = msteel*((-T+hfb)/2) 

            # overall CG
            zG = ( mZGsteel + mturb*zturb + mball*(-T+hbal/2))/mtot

            # waterplane moment of area
            Iwp = np.pi*np.power((D/2),4)/4

            # hydrostatic restoring force
            C55 = rho*g*Iwp - (msteel+mball+mturb)*zG*g - rho*g*np.pi*np.power(D,2)/4*T*T/2

            # 1DOF estimate of pitch angle
            theta = FT*hhub/C55 

        outputs['theta'] = theta
        
        
class computeLCOE(om.ExplicitComponent):
    def setup(self): 
        # global design variables
        self.add_input('D')
        self.add_input('T')
        
        # result from another component
        self.add_input('mball',val=0.0)
        self.add_input('msteel',val=0.0)
        self.add_input('theta',val=0.0)
        
        # problem constants
        self.add_input('mturb',val=1302000.0)
        self.add_input('hhub',val=120)

        self.add_output('LCOE',val = 0.0)

    
    def setup_partials(self):
        # Finite difference all partials.
        self.declare_partials('*', '*', method='fd')
        
    def compute(self, inputs, outputs):

        
        D = inputs['D'] #diameter in m
        T = inputs['T'] # draft in m
        msteel = inputs['msteel']
        mball = inputs['mball']
        hhub = inputs['hhub']
        theta = inputs['theta']
        mturb = inputs['mturb']
        
        
        f1 = 20/1E7
        f2 = 1/1E7
        f3 = 1/1E7 
        f4 = 1/1600
        
        
        LCOE = (msteel*f1 + mball*f2 + f3*theta*mturb*hhub+ f4*np.power(D,2))/(np.power(np.cos(theta),3))
        
        outputs['LCOE'] = LCOE


if __name__ == "__main__":

    model = om.Group()
    
    model.add_subsystem('mball_comp',computeMball(),
                        promotes_inputs=['D','T','rho','rho_steel',
                                        'mturb','hfb'], 
                        promotes_outputs=['mball','msteel'])
    
    model.add_subsystem('theta_comp', computeTheta(),
                        promotes_inputs=['D','T','mball','msteel','rho',
                                         'mturb','zturb','FT','hhub','hfb'], 
                        promotes_outputs=['theta'])
    
    model.add_subsystem('LCOE_comp', computeLCOE(),
                        promotes_inputs=['D','T','mball','msteel','theta',
                                         'mturb','hhub'],
                        promotes_outputs=['LCOE'])
    
    model.set_input_defaults('mturb',val=1244000.0)
    model.set_input_defaults('zturb',val=92.5)
    model.set_input_defaults('rho_steel',val=8500)
          
    prob = om.Problem(model)
    prob.model = model
    prob.driver = om.ScipyOptimizeDriver()
    prob.driver.options['optimizer'] = 'SLSQP'
    prob.setup()

    prob.set_val('D', 20)
    prob.set_val('T', 120.0)

    prob.run_model()
    
    
    
    
    
    
    

    prob.model.approx_totals()
    prob.model.add_design_var('D',lower=8,upper=30)
    prob.model.add_design_var('T',lower=50,upper=300)
    prob.model.add_constraint('mball', lower=0)
    prob.model.add_constraint('theta', lower=0,upper=45*np.pi/180)
    prob.model.set_input_defaults('D',val=10.0)
    prob.model.set_input_defaults('T',val=120.0)

    prob.model.add_objective('LCOE')
    
    prob.setup()
    prob.set_solver_print(level=0)
    
    prob.run_driver()
    
    
    print('minimum found at')
    print(prob.get_val('D'))
    print(prob.get_val('T'))

    print('minumum objective')
    print(prob.get_val('LCOE'))
    
    # om.n2(prob)

    # prob.set_val('parab_comp.x', 5.0)
    # prob.set_val('parab_comp.y', -2.0)

    # prob.run_model()
    # print(prob.get_val('parab_comp.f_xy'))

