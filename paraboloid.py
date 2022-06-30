"""
Demonstration of a model using the Paraboloid component.
"""
from __future__ import division, print_function
import openmdao.api as om

class Paraboloid(om.ExplicitComponent):
    """
    Evaluates the equation f(x,y) = (x-3)^2 + xy + (y+4)^2 - 3.
    """

    def setup(self):
        self.add_input('x', val=0.0)
        self.add_input('y', val=0.0)

        self.add_output('f_xy', val=0.0)

        # Finite difference all partials.
        self.declare_partials('*', '*', method='fd')

    def compute(self, inputs, outputs):
        """
        f(x,y) = (x-3)^2 + xy + (y+4)^2 - 3

        Minimum at: x = 6.6667; y = -7.3333
        """
        x = inputs['x']
        y = inputs['y']

        outputs['f_xy'] = (x-3.0)**2 + x*y + (y+4.0)**2 - 3.0


if __name__ == "__main__":

    model = om.Group()
    model.add_subsystem('paraboloid', Paraboloid())
    prob = om.Problem(model)
    prob.setup()

    prob.set_val('paraboloid.x', 3.0)
    prob.set_val('paraboloid.y', -4.0)
    prob.run_model()
    print('f = %2.2f at x = %2.2f, y = %2.2f' %(prob.get_val('paraboloid.f_xy'), prob.get_val('paraboloid.x'), prob.get_val('paraboloid.y')))

    prob.set_val('paraboloid.x', 5.0)
    prob.set_val('paraboloid.y', -2.0)
    prob.run_model()
    print('f = %2.2f at x = %2.2f, y = %2.2f' %(prob.get_val('paraboloid.f_xy'), prob.get_val('paraboloid.x'), prob.get_val('paraboloid.y')))       