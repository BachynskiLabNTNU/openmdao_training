## Based on OpenMDAO Training Lab 1, found here: https://github.com/OpenMDAO/openmdao_training
## Adapted from Spar FWT OpenMDAO model by Erin Bachynski-Polić

import openmdao.api as om
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

class computeMball(om.ExplicitComponent):
    '''
    A component to compute the required ballast mass for static equilibrium of spar-type FWT.  Also exports thickness of steel shell and total steel mass (assuming no stiffeners).  Values of water density, steel density, required freeboard, and wind turbine+tower mass are considered constants in the component and are defined through component options.
    '''
    def initialize(self):
        # Pass through dictionary containing constants - parameters that would be fixed during an optimization/solution
        self.options.declare('params', types=dict)

    def setup(self):
        # Design variables
        self.add_input('D', val=20.0, units='m')
        self.add_input('T', val=120.0, units='m')
        
        # Outputs (ballast mass, and steel mass and thickness)
        self.add_output('t', val=1.0, units='m')
        self.add_output('mball', val=1.0, units='kg')
        self.add_output('msteel', val=1.0, units='kg')
    
    def setup_partials(self):
        # Declare partials, must analytically define the first set.  Mixing finite difference and complex step approximations
        self.declare_partials(of='t',wrt=['T'])
        self.declare_partials(of='mball', wrt=['D','T'], method='fd')
        self.declare_partials(of='msteel',wrt=['D','T'], method='cs')
        
    def compute(self,inputs,outputs): 
        # Bring in problem constants
        params = self.options['params']
        rho = params['rho']
        rho_steel = params['rho_steel']
        mturb = params['mturb']
        hfb = params['hfb']
        # Design variables
        D = inputs['D']
        T = inputs['T']
        
        # based on T, find thickness and ballast in order to meet vertical equilibrium
        t = 0.03+T/3000
        buoy = np.pi*np.power(D,2)/4*T*rho # displacement water mass
        msteel = rho_steel*(np.pi*D*(T+hfb)*t + 2*np.pi*np.power(D,2)/4*t) # steel mass
        mball = buoy-msteel-mturb # required ballast
        
        # Raise an analysis error if the dimensions provided are too small to support the turbine
        ## NOTE - throwing analysis errors isn't best practice in building models, ideally they should be robust enough to handle 'strange' situations the optimizer can think of.  In this case our model is very simplified and we need to include this to avoid non-physical (but mathematically possible) solutions.
        if mball < 0.0:
            raise om.AnalysisError('Negative ballast! Increase dimensions')

        outputs['mball'] = mball
        outputs['msteel'] = msteel
    
    def compute_partials(self, inputs, partials):
        # Analytical definition of this partial
        partials['t','T'] = 1/3000.
        # Note other partials are estimated, and nothing is needed here

class computeTheta(om.ExplicitComponent): 
    '''
    A component to compute the mean static pitch angle of the spar-type FWT given a constant thrust force on the wind turbine.  Requires design variables, masses calculated in mBall component, and several parameters pass through as options.
    '''
    def initialize(self):
        # Pass through dictionary containing constants - parameters that would be fixed during an optimization/solution
        self.options.declare('params', types=dict)

    def setup(self): 
        # Design variables
        self.add_input('D', val=20.0, units='m')
        self.add_input('T', val=120.0, units='m')
        
        # Results from mBall component
        self.add_input('mball',val=0.0, units='kg')
        self.add_input('msteel',val=0.0, units='kg')
        self.add_input('FT', val= 0.0, units='N')
        
        # Output, static pitch angle (IN RADIANS)
        self.add_output('theta', val = 0.0, units='rad')
    
    def setup_partials(self):
        self.declare_partials('*', '*', method='fd')

    def compute(self, inputs, outputs):
        # Bring in problem constants
        params = self.options['params']
        rho = params['rho']
        rho_conc = params['rho_conc']
        mturb = params['mturb']
        zturb = params['zturb']
        hhub = params['hhub']
        hfb = params['hfb']

        # Inputs from mBall component
        mball = inputs['mball']
        msteel = inputs['msteel']
        FT = inputs['FT']
        
        # Design variables
        D = inputs['D']
        T = inputs['T']
        
        g = 9.81 # acceleration due to gravity

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

class computeThrust(om.ExplicitComponent): 
    '''
    Compute thrust force based on static pitch angle, assuming a cosine loss from a fixed maximum thrust.  Introduces a cycle in the model!
    '''
    def initialize(self):
        # Pass through dictionary containing constants - parameters that would be fixed during an optimization/solution
        self.options.declare('params', types=dict)

    def setup(self): 
        # Result from theta component
        self.add_input('theta', val = 0.0, units='rad')
                
        # Output, thrust force
        self.add_output('FT', val = 0.0, units='N')
    
    def setup_partials(self):
        # Finite difference all partials.
        self.declare_partials(of='FT', wrt='theta')
        
    def compute(self, inputs, outputs):
        # Bring in problem constants
        params = self.options['params']
        thrust_0 = params['thrust_0']
        # Pitch angle input
        theta = inputs['theta']

        # if-statement needed to avoid a negative pitch angle
        if theta < 0.0:
            outputs['FT'] = 1.e-6 * thrust_0
        else :        
            # Simple cosine loss of thrust force based on pitch angle
            outputs['FT'] = thrust_0 * np.cos(theta)
    
    def compute_partials(self, inputs, partials):
        params = self.options['params']
        thrust_0 = params['thrust_0']
        theta = inputs['theta']
        
        if theta < 0.0:
            partials['FT', 'theta'] = 1.e-6
        else :
            partials['FT', 'theta'] = -1. * thrust_0 * np.sin(theta)

class computeLCOE(om.ImplicitComponent):
    '''
    Compute a "LCOE" of the spar-type floating wind turbine.  The output is unitless and the formula is not based on specific data or research, and is meant as a demonstration only.  

    Note this is a 'fake' implicit component, we could easily define this explicitly and in fact have done just that below!
    '''
    def initialize(self):
        # Pass through dictionary containing constants - parameters that would be fixed during an optimization/solution
        self.options.declare('params', types=dict)

    def setup(self): 
        # Design variables
        self.add_input('D', val=20.0, units='m')
        self.add_input('T', val=120.0, units='m')
        
        # Results from mBall and theta components
        self.add_input('mball', val=0.0, units='kg')
        self.add_input('msteel', val=0.0, units='kg')
        self.add_input('theta', val=0.0, units='rad')
        
        # Output, a unitless LCOE
        self.add_output('LCOE', val=0.0)

    def setup_partials(self):
        # Declare FD partials (to be used if partials aren't defined below...)
        self.declare_partials(of='*', wrt='*')

    def apply_nonlinear(self, inputs, outputs, residuals):
        # Bring in problem constants
        params = self.options['params']
        hhub = params['hhub']
        mturb = params['mturb']
        # Design Variables
        D = inputs['D'] #diameter in m
        T = inputs['T'] # draft in m
        # Inputs from other components
        msteel = inputs['msteel']
        mball = inputs['mball']
        theta = inputs['theta']

        # Cost factors - can experiment with changing these
        f1 = 20/1E7
        f2 = 1/1E7
        f3 = 1/1E7 
        f4 = 1/1600

        # implement the calculation LCOE = (msteel*f1 + mball*f2 + f3*theta*mturb*hhub+ f4*np.power(D,2))
        residuals['LCOE'] = (msteel*f1 + mball*f2 + f3*theta*mturb*hhub + f4*np.power(D,2)) - outputs['LCOE']
    
    def linearize(self, inputs, outputs, partials):
        # Bring in problem constants
        params = self.options['params']
        hhub = params['hhub']
        mturb = params['mturb']
        # Design Variables
        D = inputs['D'] #diameter in m
        T = inputs['T'] # draft in m

        # Cost factors - can experiment with changing these
        f1 = 20/1E7
        f2 = 1/1E7
        f3 = 1/1E7 
        f4 = 1/1600

        partials['LCOE', 'D'] = 2. * D * f4
        partials['LCOE', 'T'] = 0.
        partials['LCOE', 'msteel'] = f1
        partials['LCOE', 'mball'] = f2
        partials['LCOE', 'theta'] = f3*mturb*hhub
        partials['LCOE', 'LCOE'] = -1.

class expComputeLCOE(om.ExplicitComponent):
    '''
    An explicit version of the LCOE component written above. 
    '''
    def initialize(self):
        # Pass through dictionary containing constants - parameters that would be fixed during an optimization/solution
        self.options.declare('params', types=dict)

    def setup(self): 
        # Design variables
        self.add_input('D', val=20.0, units='m')
        self.add_input('T', val=120.0, units='m')
        
        # Results from mBall and theta components
        self.add_input('mball', val=0.0, units='kg')
        self.add_input('msteel', val=0.0, units='kg')
        self.add_input('theta', val=0.0, units='rad')
        
        # Output, a unitless LCOE
        self.add_output('LCOE', val=0.0)

    def setup_partials(self):
        # Declare FD partials (to be used if partials aren't defined below...)
        self.declare_partials('*', '*', method='fd')
        
    def compute(self, inputs, outputs):
         # Bring in problem constants
        params = self.options['params']
        hhub = params['hhub']
        mturb = params['mturb']
        # Design Variables
        D = inputs['D'] #diameter in m
        T = inputs['T'] # draft in m
        # Inputs from other components
        msteel = inputs['msteel']
        mball = inputs['mball']
        theta = inputs['theta']
        
        # Cost factors - can experiment with changing these
        f1 = 20/1E7
        f2 = 1/1E7
        f3 = 1/1E7 
        f4 = 1/1600

        LCOE = (msteel*f1 + mball*f2 + f3*theta*mturb*hhub + f4*np.power(D,2))       
        outputs['LCOE'] = LCOE

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
        # 'FT':  1500000.0, # N --- Thrust force is no longer constant!
        'thrust_0': 1500000.0, # N - Using a new 'baseline' thrust variable
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
        promotes_inputs=['D','T','mball','msteel','FT'], 
        promotes_outputs=['theta'])
    model.add_subsystem('thrust_comp', 
        computeThrust(params=fwt_params),
        promotes_inputs=['theta'], 
        promotes_outputs=['FT'])
    model.add_subsystem('LCOE_comp', 
        # expComputeLCOE(params=fwt_params),
        computeLCOE(params=fwt_params),
        promotes_inputs=['D','T','mball','msteel','theta'],
        promotes_outputs=['LCOE'])

    # Connect model to problem     
    prob = om.Problem(model=model, name='lab_1', reports='n2')
    prob.model = model

    # pick a solver: 'newton', 'broyden', 'nlbgs', or 'nlbjac'
    # must define a nonlinear solver since this system has a cycle
    solver_flag = 'newton'

    if solver_flag == 'newton':
        prob.model.nonlinear_solver=om.NewtonSolver(iprint=2)
        # solve_subsystems should almost always be turned on
        # it improves solver robustness
        prob.model.nonlinear_solver.options['solve_subsystems'] = True
        prob.model.nonlinear_solver.options['maxiter'] = 100
        # these options control how tightly the solver converges the system
        prob.model.nonlinear_solver.options['atol'] = 1e-8
        prob.model.nonlinear_solver.options['rtol'] = 1e-8
        # the Newton solver requires a linear solver
        prob.model.linear_solver = om.DirectSolver()

    elif solver_flag == 'broyden':
        prob.model.nonlinear_solver=om.BroydenSolver(iprint=2)
        prob.model.nonlinear_solver.options['compute_jacobian'] = True
        prob.model.nonlinear_solver.options['maxiter'] = 100
        # these options control how tightly the solver converges the system
        prob.model.nonlinear_solver.options['atol'] = 1e-8
        prob.model.nonlinear_solver.options['rtol'] = 1e-8
        # the Broyden solver requires a linear solver *if* options['compute_jacobian'] = True
        prob.model.linear_solver = om.DirectSolver()

    elif solver_flag == 'nlbgs':
        # The nonlinear block Gauss-Seidel solver is an iterative solvver
        # Requires no linear solver and works even without derivatives
        prob.model.nonlinear_solver=om.NonlinearBlockGS(iprint=2)
        prob.model.nonlinear_solver.options['maxiter'] = 400
        prob.model.nonlinear_solver.options['atol'] = 1e-8
        prob.model.nonlinear_solver.options['rtol'] = 1e-8
        # The Aitken relaxation method improves robustness at cost of some speed
        prob.model.nonlinear_solver.options['use_aitken'] = False
        prob.model.nonlinear_solver.options['use_apply_nonlinear'] = True

    elif solver_flag == 'nlbjac':
        # We don't usually recommend using nonlinear block Jacobi as it converges slower
        prob.model.nonlinear_solver=om.NonlinearBlockJac(iprint=2)
        prob.model.nonlinear_solver.options['maxiter'] = 400
        prob.model.nonlinear_solver.options['atol'] = 1e-8
        prob.model.nonlinear_solver.options['rtol'] = 1e-8

    else:
        raise ValueError("bad solver selection!")

    prob.setup()
    ### If using the Newton solver you should generally check your partial derivatives before you run the model. It won't converge if you made a mistake.
    prob.check_partials(compact_print=True, form='central')

    # Set value of design variables
    # Change these inputs to see the effect on results
    prob.set_val('D',val=20.)
    prob.set_val('T',val=120.)

    # Run model
    prob.run_model()

    # Debugging printouts
    print('Diameter: %2.2f m' %prob.get_val('D'))
    print('Draft: %2.2f m' %prob.get_val('T'))
    print('Static Pitch Angle: %3.3f deg' %(prob.get_val('theta')*180/np.pi))
    
    # # --- Debugging contour plots ---
    # x1s = np.linspace(8,40,50)
    # x2s = np.linspace(40,160,100)
    # fLCOE = np.zeros((len(x1s),len(x2s)))
    # ftheta = np.zeros((len(x1s),len(x2s)))

    # for ii in range(0,len(x1s)):  
    #     for jj in range(0,len(x2s)): 
    #         # print(x1s[ii],x2s[jj])
    #         prob.set_val('D',x1s[ii])
    #         prob.set_val('T',x2s[jj])
    #         prob.run_model()
    #         fLCOE[ii,jj] = prob.model.LCOE_comp.get_val('LCOE')
    #         ftheta[ii,jj] = prob.model.theta_comp.get_val('theta')*180/np.pi
    #         if (prob.model.theta_comp.get_val('theta')*180/np.pi) == 0.0 :
    #             print('Platform (D = %3.2f, T = %3.2f) is unstable!' %(x1s[ii],x2s[jj]))
    #         # print(x1s[ii],x2s[jj],fLCOE[ii,jj], ftheta[ii,jj])


    # fig, (ax1, ax2) = plt.subplots(1, 2)
    # csLCOE = ax1.contourf(x1s,x2s,np.transpose(fLCOE),np.linspace(1,200,100))
    # cstheta = ax2.contourf(x1s,x2s,np.transpose(ftheta),np.linspace(0,45,100))
    # fig.colorbar(csLCOE,ax=ax1)
    # fig.colorbar(cstheta,ax=ax2)

    # plt.show()
