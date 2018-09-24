def l2_err_sqr(cgp, x_data, y_data, pars):
	out = 0.0
	n = len(x_data)

	for i in range(n):
		diff = y_data[i] - cgp.eval([x_data[i]], pars)
		out += diff*diff
	return out

def gradient(cgp, x_data, y_data, pars, h=1.0e-10):
	n = len(x_data)
	cgp_data = [cgp.eval([x_data[i]], pars) for i in range(n)]

	gradients = [0.0]*len(pars)
	for par_nr in range(len(pars)):
		grad = 0.0
		pars_plus_h = list(pars)
		pars_plus_h[par_nr] += h
		for i in range(n):
			derivative = (cgp.eval([x_data[i]], pars_plus_h)-cgp_data[i])/h
			grad += (cgp_data[i]-y_data[i])*derivative
		grad = 2.0*grad
		gradients[par_nr] = grad

	return gradients

def grad_descent(cgp, x_data, y_data, max_iter=100):
	#TODO: This is a stupid function by a stupid person.
	assert len(x_data) == len(y_data)

	pars = [1.0e-2]*cgp.nr_of_used_pars
	old_err = l2_err_sqr(cgp, x_data, y_data, pars)

	step_size = 1.0e-1 # TODO: Make adaptive.

	for itr in range(max_iter):
		grad = gradient(cgp, x_data, y_data, pars)

		new_pars = [pars[i]-grad[i]*step_size for i in range(len(pars))]
		new_err = l2_err_sqr(cgp, x_data, y_data, new_pars)

		if new_err < old_err:
			#print("err:", new_err)
			old_err = new_err
			pars = new_pars
		else:
			break
	return pars