from full_search import full_search
from tune_compressor import objective_function
from operation_table import op_table

def tune_compressor(max_nr_of_nodes, max_err_allowed, y_data_list, freq=1.0, resolution=64, nr_of_parameters=4, max_nr_of_funcs=4):
	best_cgp_list = []
	best_compression_rate = 1.0
	cgp_dims = 1

	for itr in range(max_nr_of_funcs):
		# Do a full search.
		obj_func_curry = lambda new_cgp: objective_function(new_cgp, best_cgp_list, max_err_allowed, y_data_list, freq, resolution, nr_of_parameters, op_table)
		new_gene, new_err = full_search(max_nr_of_nodes, y_data_list, nr_of_parameters, obj_func_curry)
		new_compression_rate = 1.0/new_err
		new_cgp = CGP(1, op_table, new_gene, nr_of_parameters=nr_of_parameters)

		if new_compression_rate	> best_compression_rate:
			best_compression_rate = new_compression_rate
			best_cgp_list.append(new_cgp)
		else:
			break

	return best_cgp_list, best_cgp_list


if __name__ == '__main__':
	from math import sin, cos, pi, fabs
	from random import gauss, seed

	seed(0)

	data1 = [sin(t*pi/10.3)+0.2*cos(t*pi*2.3546)+t*sin(t)/700.0 + gauss(0, 0.001) for t in range(400)]
	data2 = [0.0]*500

	data2[0] = 0.5
	for i in range(1, len(data2)):
		data2[i] = 2.9*fabs(data2[i-1]*(1.0-data2[i-1]))
		data2[i] += gauss(0, 0.001)
	#import matplotlib.pyplot as plt
	#plt.plot(data2)
	#plt.show()

	max_err_allowed	= 1.0e-2
	max_nr_of_nodes = 2
	tune_compressor(max_nr_of_nodes, max_err_allowed, [data1, data2] ,freq=1.0e-2, resolution=64, nr_of_parameters=4, max_nr_of_funcs=3)
