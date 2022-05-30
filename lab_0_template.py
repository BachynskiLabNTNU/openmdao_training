## Based on OpenMDAO Training Lab 0, found here: https://github.com/OpenMDAO/openmdao_training
## Adapted from Spar FWT OpenMDAO model by Erin Bachynski-Polić

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
        # TODO: define inputs for diameter and draft:
        # use the following format
        self.add_input(<name>, <default value>, units=<unit string>)
        
        # Outputs (ballast mass, and steel mass and thickness)
        # TODO: define outputs for ballast mass, steel mass and steel thickness:
        # use the following format
        self.add_output(<name>, <default value>, units=<unit string>)
        
    def compute(self,inputs,outputs): 
        # Bring in problem constants
        params = self.options['params']
        rho = params['rho']
        # TODO: connect additional parameters needed
        
        D = inputs['D']
        T = inputs['T']
        
        # based on T, find thickness and ballast in order to meet vertical equilibrium
        t = 0.03+T/3000
        # TODO: calculate buoyancy (displaced volume times water density)
        buoy = 
        # Steel mass calculation is already setup properly
        msteel = rho_steel*(np.pi*D*(T+hfb)*t + 2*np.pi*np.power(D,2)/4*t) # steel mass
        # TODO: calculate required ballast mass (buoyancy minus steel mass and turbine mass)
        mball = 

        # Raise an analysis error if the dimensions provided are too small to support the turbine
        if mball < 0.0:
            # TODO: Consider why this is raised as an error instead of added as a constraint (no code changes needed)
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
        # TODO: define inputs for diameter and draft:
        # use the following format
        self.add_input(<name>, <default value>, units=<unit string>)

        # TODO: define inputs for mball and msteel
        
        # Outputs (ballast mass, and steel mass and thickness)
        # TODO: define output for pitch angle:
        # use the following format
        self.add_output(<name>, <default value>, units=<unit string>)

    def compute(self, inputs, outputs):
        # TODO: nothing to do here! the computation is already setup for you
        # Bring in problem constants
        params = self.options['params']
        rho = params['rho']
        rho_conc = params['rho_conc']
        mturb = params['mturb']
        zturb = params['zturb']
        FT = params['FT']
        hhub = params['hhub']
        hfb = params['hfb']
        g = 9.81 # acceleration due to gravity

        # Inputs from mBall component
        mball = inputs['mball']
        msteel = inputs['msteel']
        
        # Design variables
        D = inputs['D']
        T = inputs['T']        

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
        # TODO: Change turbine properties if you're curious!
        'mturb': 1244000.0, # kg
        'zturb': 92.5, # m
        'FT':  1500000.0, # N
        'hhub': 120., # m
        'hfb': 10.0, # m
    }
    # Add each component as a subsystem, defining the inputs, outputs, and options
    model.add_subsystem('mball_comp',
        computeMball(params=fwt_params),
        # TODO: Promote the appropriate inputs and outputs OR add subsystems and make variable connections as needed
        promotes_inputs=[], 
        promotes_outputs=[])
    # TODO: Add a subsystem for the pitch angle calculation, based on on the ballast component

    # Connect model to a problem and run the problem (just solving the model)    
    prob = om.Problem(model, name='lab_0', reports='n2') # There are many other ways to generate N2 diagrams in OpenMDAO!
    prob.model = model
    prob.setup()
    
    # Set value of design variables
    # TODO: Change these inputs to see the effect on results
    prob.set_val('D',val=20.)
    prob.set_val('T',val=120.)
    
    # Run model
    prob.run_model()
    
    # Debugging printouts
    print('Diameter: %2.2f m' %prob.get_val('D'))
    print('Draft: %2.2f m' %prob.get_val('T'))
    print('Static Pitch Angle: %3.3f deg' %(prob.get_val('theta')*180/np.pi))