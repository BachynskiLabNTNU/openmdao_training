## Adapted from Spar FWT OpenMDAO model by Erin Bachynski 

import openmdao.api as om
import numpy as np

class computeMball(om.ExplicitComponent):
    def setup(self):
         # global design variables
        self.add_input('D', val=1.0, units='m')
        self.add_input('T', val=1.0, units='m')
        
        # problem constants
        self.add_input('rho', val=1025.0, units='kg/m**3')
        self.add_input('rho_steel', val=7850.0, units='kg/m**3')
        self.add_input('mturb', val=1244000.0, units='kg')
        self.add_input('hfb', val=10.0, units='m')
        
        self.add_output('t', val=1.0, units='m')
        self.add_output('mball', val=1.0, units='kg')
        self.add_output('msteel', val=1.0, units='kg')
        
    def compute(self,inputs,outputs): 
        rho = inputs['rho']
        rho_steel = inputs['rho_steel']
        mturb = inputs['mturb']
        hfb = inputs['hfb']
        
        D = inputs['D']
        T = inputs['T']
        
        # based on T, find thickness and ballast in order to meet vertical equilibrium
        t = 0.03+T/3000
        buoy = np.pi*np.power(D,2)/4*T*rho # displacement water mass
        msteel = rho_steel*(np.pi*D*(T+hfb)*t+2*np.pi*np.power(D,2)/4*t) # steel mass
        mball = buoy-msteel-mturb # required ballast

        outputs['mball'] = mball
        outputs['msteel'] = msteel


class computeTheta(om.ExplicitComponent): 
    #T,rho,rho_steel,rho_conc, mturb,zturb,FT,hhub,hfb
    def setup(self): 
        # global design variables
        self.add_input('D', val=1.0, units='m')
        self.add_input('T', val=1.0, units='m')
        
        # result from another component
        self.add_input('mball',val=0.0, units='kg')
        self.add_input('msteel',val=0.0, units='kg')
        
        # problem constants
        self.add_input('rho',val=1025.0, units='kg/m**3')
        self.add_input('rho_conc', val=2650.0, units='kg/m**3')
        self.add_input('mturb',val=1244000.0, units='kg')
        self.add_input('zturb',val=92.5, units='m')
        self.add_input('FT',val= 1500000.0, units='N')
        self.add_input('hhub', val=120., units='m')
        self.add_input('hfb', val=10.0, units='m')
        
        # outputs
        self.add_output('theta', val = 0.0, units='rad')

        
    def compute(self, inputs, outputs):
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

        # we can't have negative ballast, the function returns zero and 90 degree theta
        if mball < 0: 
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

if __name__ == "__main__":

    # Define the model as a single group, and add subsystems.
    model = om.Group()
    model.add_subsystem('mball_comp',
        computeMball(),
        promotes_inputs=['D','T','rho','rho_steel','mturb','hfb'], 
        promotes_outputs=['mball','msteel'])
    model.add_subsystem('theta_comp', 
        computeTheta(),
        promotes_inputs=['D','T','mball','msteel','rho','mturb','zturb','FT','hhub','hfb'], 
        promotes_outputs=['theta'])
    
    # Change these inputs to see the effect on results
    model.set_input_defaults('D',val=15.)
    model.set_input_defaults('T',val=90.)
    model.set_input_defaults('mturb',val=1244000.0)
    model.set_input_defaults('zturb',val=92.5)
    model.set_input_defaults('rho_steel',val=8500)

    # Run problem (just the model in this case)      
    prob = om.Problem(model)
    prob.model = model
    prob.setup()
    prob.run_model()
    
    # Debugging printouts
    print('Diameter: %2.2f m' %prob.get_val('D'))
    print('Draft: %2.2f m' %prob.get_val('T'))
    print('Static Pitch Angle: %3.3f deg' %(prob.get_val('theta')*180/np.pi))
    
    # Create N2 diagram
    om.n2(prob)


