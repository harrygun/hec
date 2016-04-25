import numpy as np
import sympy as sp
import sympy.mpmath as mt
import scipy.integrate as integ


#_mx_step_=500
_mx_step_=1000


def ode_solver(derv, xi, yi, x, args=(), solver_type='scipy.integrate.odeint', **kwargs):
    ''' ->> solver wrapper <<- 
        derv(p, x, y, *args, **kwargs):
    '''

    p = args[0]
    if solver_type=='sympy.mpmath.odefun':
        return map( mt.odefun(lambda x, y: derv(p, x, y), xi, yi, **kwargs), x)

    elif solver_type=='scipy.integrate.odeint':
        return integ.odeint(lambda y, x, p: derv(p, x, y), yi, x, args=(p,), mxstep=_mx_step_)




def ode_solver_general(derv, xi, yi, x, args=(), solver_type='scipy.integrate.odeint', **kwargs):
    ''' ->> solver wrapper <<- 
        derv(p, x, y, *args, **kwargs):
    '''

    p = args[0]
    other_args=args[1:]

    if solver_type=='sympy.mpmath.odefun':
        #return map( mt.odefun(lambda x, y: derv(p, x, y), xi, yi, **kwargs), x)
	raise Exception('sympy integrator NOT supported yet.')

    elif solver_type=='scipy.integrate.odeint':

        #return integ.odeint(lambda y, x, _args_: derv(p, x, y, other_args), yi, x, args=tuple([p]+list(other_args)), mxstep=_mx_step_)
        #return integ.odeint(lambda y, x, _args_: derv(p, x, y, other_args), yi, x, args=args, mxstep=_mx_step_)

        return integ.odeint(lambda y, x, p, other_args: derv(p, x, y, other_args), yi, x, \
                            args=(p, other_args, ), mxstep=_mx_step_)


