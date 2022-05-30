import numpy as np

import openmdao.api as om
import matplotlib.pyplot as plt

# all of these components have already been created for you, but look in beam_comp.py if you're curious to see how
from beam_comps import (MomentOfInertiaComp, LocalStiffnessMatrixComp, FEM, ComplianceComp, VolumeComp)

# Note here we are creating a group, rather than a component!  It functions in a similar way, but there is no calculation happening here, only connection of existing components
class BeamGroup(om.Group):

    def initialize(self):
        # Declare options individually (before we have done this using a dictionary!)
        self.options.declare('E')
        self.options.declare('L')
        self.options.declare('b')
        self.options.declare('volume')
        self.options.declare('num_elements', int)

    def setup(self):
        E = self.options['E']
        L = self.options['L']
        b = self.options['b']
        volume = self.options['volume']
        num_elements = self.options['num_elements']
        num_nodes = num_elements + 1

        force_vector = np.zeros(2 * num_nodes)
        force_vector[-2] = -1.

        I_comp = MomentOfInertiaComp(num_elements=num_elements, b=b)
        self.add_subsystem('I_comp', 
            I_comp,
            promotes_inputs=['h'],
            promotes_outputs=['I'])

        k_comp = LocalStiffnessMatrixComp(num_elements=num_elements, E=E, L=L)
        self.add_subsystem('local_stiffness_matrix_comp', 
            k_comp,
            promotes_inputs=['I'],
            promotes_outputs=['K_local'])

        fem_comp = FEM(num_elements=num_elements, force_vector=force_vector)
        self.add_subsystem('FEM', 
            fem_comp,
            promotes_inputs=['K_local'],
            promotes_outputs=['u'])

        # connection for: FEM -> compliance
        # this one is tricky, because you just want the states from the nodes, but not the last 2 which relate to the clamped boundary condition on the left
        self.connect('u', 'compliance_comp.displacements', src_indices=np.arange(2*num_nodes))
        
        compli_comp = ComplianceComp(num_elements=num_elements, force_vector=force_vector)
        self.add_subsystem('compliance_comp', 
            compli_comp,
            promotes_inputs=[],
            promotes_outputs=['compliance'])

        v_comp = VolumeComp(num_elements=num_elements, b=b, L=L)
        self.add_subsystem('volume_comp', 
            v_comp,
            promotes_inputs=['h'],
            promotes_outputs=['volume'])

if __name__ == "__main__":

    import time

    E = 1.
    L = 1.
    b = 0.1
    volume = 0.01

    num_elements = 50

    prob = om.Problem(name='lab_2', reports='n2')
    prob.model.add_subsystem('beam', BeamGroup(E=E, L=L, b=b, volume=volume, num_elements=num_elements))

    # Add a bounded design variable (for thickness)
    prob.model.add_design_var('beam.h', lower=1e-2, upper=10.)
    # The goal is to minimize compliance (force times displacement)
    prob.model.add_objective('beam.compliance')
    # Maintain beam volume equal to set value
    prob.model.add_constraint('beam.volume', equals=volume)
    
    # Use the SciPy optimize functions as the problem's driver (this implements optimization)
    prob.driver = om.ScipyOptimizeDriver()
    prob.driver.options['optimizer'] = 'SLSQP'
    prob.driver.options['tol'] = 1e-9
    prob.driver.options['disp'] = True

    ## --- Recording Example ---
    recorder = om.SqliteRecorder('beam.sql')
    prob.driver.add_recorder(recorder)
    prob.driver.recording_options['includes'] = ['*']
    prob.add_recorder(recorder)
    prob.recording_options['includes'] = ['*']
    # This could then be read later to understand more about the optimization

    ######################################################
    # Use top level FD or CS to approximate derivatives
    ######################################################
    # prob.model.approx_totals(method='fd', step_calc="rel", step=1e-3)
    # prob.model.approx_totals(method='cs')

    # Setup, add timer and run the optimization
    prob.setup()
    start_time = time.time()
    prob.run_driver()
    print('Optimization time: %2.3f s' %(time.time()-start_time))
    print('Beam Thickness:')
    print(prob['beam.h'])

    # --- Plot thickness result---
    fig, ax = plt.subplots()
    csLCOE = ax.plot(np.linspace(0.,2.,num_elements),prob['beam.h'],'-ob')
    ax.set_xlabel('x')
    ax.set_ylabel('optimized thickness')
    plt.show()