from compression import get_compression_rate
from cgp import CGP

def neg_fitness(max_err_allowed, cgps, y_data_list, freq, resolution):
	total_data_in = sum(len(y_data)*resolution for y_data in y_data_list)
	total_data_out = 0.0
	for data in y_data_list:
		compress_rate = get_compression_rate(cgps, data, max_err_allowed, freq=freq, resolution=resolution)
		current_data_in = len(data)*resolution
		total_data_out += current_data_in*compress_rate
	total_compression_rate = total_data_out/total_data_in

	fitness = 1.0/total_compression_rate
	return fitness

def objective_function(new_cgp_gene, old_cgps, max_err_allowed, y_data_list, freq, resolution, nr_of_parameters, op_table):
	new_cgp = CGP(1, op_table, new_cgp_gene, nr_of_parameters=nr_of_parameters)

	# TODO: Break if the CGP is constant.

	cgps = [new_cgp] + old_cgps
	return neg_fitness(max_err_allowed, cgps, y_data_list, freq, resolution)

