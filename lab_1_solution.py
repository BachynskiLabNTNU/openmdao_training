## Adapted from Spar FWT OpenMDAO model by Erin Bachynski 

import openmdao.api as om
import numpy as np

class computeMball(om.ExplicitComponent):
    def setup(self):
        # Design variables
        self.add_input('D', val=1.0, units='m')
        self.add_input('T', val=1.0, units='m')
        
        # Problem constants, can still be modified
        self.add_input('rho', val=1025.0, units='kg/m**3')
        self.add_input('rho_steel', val=7850.0, units='kg/m**3')
        self.add_input('mturb', val=1244000.0, units='kg')
        self.add_input('hfb', val=10.0, units='m')
        
        # Outputs (ballast thickness, mass, and steel mass)
        self.add_output('t', val=1.0, units='m')
        self.add_output('mball', val=1.0, units='kg')
        self.add_output('msteel', val=1.0, units='kg')
        
        
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
        t = 0.03+T/3000
        buoy = np.pi*np.power(D,2)/4*T*rho # displacement water mass
        msteel = rho_steel*(np.pi*D*(T+hfb)*t+2*np.pi*np.power(D,2)/4*t) # steel mass
        mball = buoy-msteel-mturb # required ballast

        outputs['mball'] = mball
        outputs['msteel'] = msteel


class computeTheta(om.ExplicitComponent): 
    def setup(self): 
        # Design variables
        self.add_input('D', val=1.0, units='m')
        self.add_input('T', val=1.0, units='m')
        
        # Result from mBall and Thrust component
        self.add_input('mball', val=0.0, units='kg')
        self.add_input('msteel', val=0.0, units='kg')
        self.add_input('FT', val= 0.0, units='N')
        
        # Problem constants, can still be modified
        self.add_input('rho',val=1025.0, units='kg/m**3')
        self.add_input('rho_conc', val=2650.0, units='kg/m**3')
        self.add_input('mturb',val=1244000.0, units='kg')
        self.add_input('zturb',val=92.5, units='m')
        self.add_input('hhub', val=120., units='m')
        self.add_input('hfb', val=10.0, units='m')
        
        # Output, static pitch angle (IN RADIANS)
        self.add_output('theta', val = 0.0, units='rad')
    
    def setup_partials(self):
        # Finite difference all partials.
        self.declare_partials('*', '*', method='fd')
        
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
        
class computeThrust(om.ExplicitComponent): 
    def setup(self): 
        # Result from theta component
        self.add_input('theta', val = 0.0, units='rad')
        
        # Problem constants, can still be modified
        self.add_input('thrust_0',val= 1500000.0, units='N')
        
        # Output, thrust force
        self.add_output('FT', val = 0.0, units='N')
    
    def setup_partials(self):
        # Finite difference all partials.
        self.declare_partials('*', '*', method='fd')
        
    def compute(self, inputs, outputs):
        thrust_0 = inputs['thrust_0']
        theta = inputs['theta']
        
        outputs['FT'] = thrust_0 * np.cos(theta)
        
class computeLCOE(om.ExplicitComponent):
    def setup(self): 
        # Design variables
        self.add_input('D', val=1.0, units='m')
        self.add_input('T', val=1.0, units='m')
        
        # Results from mBall and theta components
        self.add_input('mball', val=0.0, units='kg')
        self.add_input('msteel', val=0.0, units='kg')
        self.add_input('theta', val=0.0, units='rad')
        
        # Problem constants, can still be modified
        self.add_input('mturb', val=1302000.0, units='kg')
        self.add_input('hhub', val=120, units='m')

        # Output, a unitless LCOE
        self.add_output('LCOE', val=0.0)

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
        
        # Cost factors - can experiment with changing these
        f1 = 20/1E7
        f2 = 1/1E7
        f3 = 1/1E7 
        f4 = 1/1600

        LCOE = (msteel*f1 + mball*f2 + f3*theta*mturb*hhub+ f4*np.power(D,2))/(np.power(np.cos(theta),3))        
        outputs['LCOE'] = LCOE

class implicitComputeLCOE(om.ImplicitComponent):
    ## Note this is a 'fake' implicit component, we could easily define this explicitly and in fact have done just that below!
    def setup(self): 
        # Design variables
        self.add_input('D', val=1.0, units='m')
        self.add_input('T', val=1.0, units='m')
        
        # Results from mBall and theta components
        self.add_input('mball', val=0.0, units='kg')
        self.add_input('msteel', val=0.0, units='kg')
        self.add_input('theta', val=0.0, units='rad')
        
        # Problem constants, can still be modified
        self.add_input('mturb', val=1302000.0, units='kg')
        self.add_input('hhub', val=120, units='m')

        # Output, a unitless LCOE
        self.add_output('LCOE', val=0.0)

        # Declare FD partials (to be used if partials aren't defined below...)
        self.declare_partials('*', '*', method='fd')

    def apply_nonlinear(self, inputs, outputs, residuals):
        D = inputs['D'] #diameter in m
        T = inputs['T'] # draft in m
        msteel = inputs['msteel']
        mball = inputs['mball']
        hhub = inputs['hhub']
        theta = inputs['theta']
        mturb = inputs['mturb']

        # Cost factors - can experiment with changing these
        f1 = 20/1E7
        f2 = 1/1E7
        f3 = 1/1E7 
        f4 = 1/1600

        # implement the calculation LCOE = (msteel*f1 + mball*f2 + f3*theta*mturb*hhub+ f4*np.power(D,2))
        residuals['LCOE'] = (msteel*f1 + mball*f2 + f3*theta*mturb*hhub+ f4*np.power(D,2)) - outputs['LCOE']
    
    def linearize(self, inputs, outputs, partials):
        D = inputs['D'] #diameter in m
        T = inputs['T'] # draft in m
        msteel = inputs['msteel']
        mball = inputs['mball']
        hhub = inputs['hhub']
        theta = inputs['theta']
        mturb = inputs['mturb']

        # Cost factors - can experiment with changing these
        f1 = 20/1E7
        f2 = 1/1E7
        f3 = 1/1E7 
        f4 = 1/1600

        partials['LCOE', 'D'] = 2. * D * f4
        partials['LCOE', 'T'] = 0.
        partials['LCOE', 'msteel'] = f1
        partials['LCOE', 'mball'] = f2
        partials['LCOE', 'hhub'] = f3*theta*mturb
        partials['LCOE', 'theta'] = f3*mturb*hhub
        partials['LCOE', 'mturb'] = f3*theta*hhub


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
    model.add_subsystem('thrust_comp', 
        computeThrust(),
        promotes_inputs=['theta', 'thrust_0'], 
        promotes_outputs=['FT'])
    model.add_subsystem('LCOE_comp', 
        # computeLCOE(),
        implicitComputeLCOE(),
        promotes_inputs=['D','T','mball','msteel','theta','mturb','hhub'],
        promotes_outputs=['LCOE'])
    
    # Change these inputs to see the effect on results
    model.set_input_defaults('D',val=15.)
    model.set_input_defaults('T',val=90.)
    model.set_input_defaults('mturb',val=1244000.0)
    model.set_input_defaults('zturb',val=92.5)
    model.set_input_defaults('rho_steel',val=8500)
    model.set_input_defaults('thrust_0',val=1500000.0)

    # Run problem (just the model in this case)      
    prob = om.Problem(model)
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
        # TODO: Try using broyden with and without a computed jacobian. What happens?
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
    prob.run_model()

    # Debugging printouts
    print('Diameter: %2.2f m' %prob.get_val('D'))
    print('Draft: %2.2f m' %prob.get_val('T'))
    print('Static Pitch Angle: %3.3f deg' %(prob.get_val('theta')*180/np.pi))
    
    # Create N2 diagram
    om.n2(prob)