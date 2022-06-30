import openmdao.api as om
from paraboloid import Paraboloid

# build the model
prob = om.Problem()

prob.model.add_subsystem('paraboloid', Paraboloid())

# setup the optimization
prob.driver = om.ScipyOptimizeDriver()
prob.driver.options['optimizer'] = 'SLSQP'

prob.model.add_design_var('paraboloid.x', lower=-50, upper=50)
prob.model.add_design_var('paraboloid.y', lower=-50, upper=50)
prob.model.add_objective('paraboloid.f_xy')

prob.setup()

# Set initial values.
prob.set_val('paraboloid.x', 3.0)
prob.set_val('paraboloid.y', -4.0)

# run the optimization
prob.run_driver()
print('Minimum: %2.2f at x = %2.2f, y = %2.2f' %(prob.get_val('paraboloid.f_xy'), prob.get_val('paraboloid.x'), prob.get_val('paraboloid.y')))
om.n2(prob)