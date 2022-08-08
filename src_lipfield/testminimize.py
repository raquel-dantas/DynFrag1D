from scipy.optimize import minimize, rosen, rosen_der

# x0 = [1.3, 0.7, 0.8, 1.9, 1.2]
# res = minimize(rosen, x0, method='Nelder-Mead', tol=1e-6)
# pass


const =({'type': 'ineq', 'fun': lambda x:  x-1})

res = minimize(lambda x: x**4-3*x**2+2, -2, method='SLSQP', bounds=[(-3,3)], constraints=const, tol=1e-8)
pass