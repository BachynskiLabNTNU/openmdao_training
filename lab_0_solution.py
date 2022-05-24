## Based on OpenMDAO Training Lab 0, found here: https://github.com/OpenMDAO/openmdao_training
## Adapted from Spar FWT OpenMDAO model by Erin Bachynski 

import openmdao.api as om
import numpy as np

class computeMball(om.ExplicitComponent):
    '''
    A component to compute the required ballast mass for static equilibrium of spar-type FWT.  Also exports thickness of steel shell and total steel mass (assuming no stiffeners).  Values of water density, steel density, required freeboard, and wind turbine+tower mass are considered constants in the component and are defined through component options.
    '''
    def initialize(self):
        # Pass through dictionary containing constants - parameters that would be fixed during an optimization/solution
        self.options.declare('params', types=dict)

    def setup(self):
        # Design variables
        self.add_input('D', val=1.0, units='m')
        self.add_input('T', val=1.0, units='m')
        
        # Outputs (ballast mass, and steel mass and thickness)
        self.add_output('t', val=1.0, units='m')
        self.add_output('mball', val=1.0, units='kg')
        self.add_output('msteel', val=1.0, units='kg')
        
    def compute(self,inputs,outputs): 
        # Bring in problem constants
        params = self.options['params']
        rho = params['rho']
        rho_steel = params['rho_steel']
        mturb = params['mturb']
        hfb = params['hfb']
        
        D = inputs['D']
        T = inputs['T']
        
        # based on T, find thickness and ballast in order to meet vertical equilibrium
        t = 0.03+T/3000
        buoy = np.pi*np.power(D,2)/4*T*rho # displacement water mass
        msteel = rho_steel*(np.pi*D*(T+hfb)*t + 2*np.pi*np.power(D,2)/4*t) # steel mass
        mball = buoy-msteel-mturb # required ballast

        # Raise an analysis error if the dimensions provided are too small to support the turbine
        if mball < 0.0:
            raise om.AnalysisError('Negative ballast! Increase dimensions')

        outputs['mball'] = mball
        outputs['msteel'] = msteel

class computeTheta(om.ExplicitComponent): 
    '''
    A component to compute the mean static pitch angle of the spar-type FWT given a constant thrust force on the wind turbine.  Requires design variables, masses calculated in mBall component, and several parameters pass through as options.
    '''
    def initialize(self):
        # Pass through dictionary containing constants - parameters that would be fixed during an optimization/solution
        self.options.declare('params', types=dict)

    def setup(self): 
        # Design variables
        self.add_input('D', val=1.0, units='m')
        self.add_input('T', val=1.0, units='m')
        
        # Results from mBall component
        self.add_input('mball',val=0.0, units='kg')
        self.add_input('msteel',val=0.0, units='kg')
        
        # Output, static pitch angle (IN RADIANS)
        self.add_output('theta', val = 0.0, units='rad')

    def compute(self, inputs, outputs):
        # Bring in problem constants
        params = self.options['params']
        rho = params['rho']
        rho_conc = params['rho_conc']
        mturb = params['mturb']
        zturb = params['zturb']
        FT = params['FT']
        hhub = params['hhub']
        hfb = params['hfb']

        # Inputs from mBall component
        mball = inputs['mball']
        msteel = inputs['msteel']
        
        # Design variables
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
    # Create a dictionary with the options (constants)
    # Note this could be done for each constant individually, and that would give you a bit more control 
    fwt_params = {
        'rho': 1025.0, # kg/m**3
        'rho_steel': 8500.0, # kg/m**3
        'rho_conc': 2650.0, # kg/m**3
        'mturb': 1244000.0, # kg
        'zturb': 92.5, # m
        'FT':  1500000.0, # N
        'hhub': 120., # m
        'hfb': 10.0, # m
    }
    # Add each component as a subsystem, defining the inputs, outputs, and options
    model.add_subsystem('mball_comp',
        computeMball(params=fwt_params),
        promotes_inputs=['D','T'], 
        promotes_outputs=['mball','msteel'])
    model.add_subsystem('theta_comp', 
        computeTheta(params=fwt_params),
        promotes_inputs=['D','T','mball','msteel'], 
        promotes_outputs=['theta'])
    
    # Set value of design variables
    # Change these inputs to see the effect on results
    model.set_input_defaults('D',val=15.)
    model.set_input_defaults('T',val=90.)

    # Connect model to a problem and run the problem (just solving the model)    
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


