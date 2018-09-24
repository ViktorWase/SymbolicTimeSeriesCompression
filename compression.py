from math import fabs

from gradient_descent import grad_descent

def tune_parameters(cgp, buff):
	#TODO: Write.
	return [0.0]*cgp.nr_of_used_pars

def nr_of_output_bits(nr_of_funcs, resolution, nr_of_pars, power):
	# TODO: Add the number of output points to the output bitstream.

	#input_bits = resolution*nr_of_output_pnts
	output_bits = resolution*nr_of_pars + power
	return output_bits

def get_compression_rate(funcs, data, max_err_allowed, freq=1.0, resolution=64):
	buff = []
	n = len(data)
	nr_of_funcs = len(funcs)

	power = 0
	while 2**power < nr_of_funcs+1: # The +1 is because we can output a single point.
		power += 1

	compression_rate_of_1_pnt_in_1_pnt_out = resolution / float(resolution+power)
	assert compression_rate_of_1_pnt_in_1_pnt_out < 1.0

	total_output_bits = 0

	compression_rate_of_buff_per_func = [0.0]*nr_of_funcs
	input_pnts_per_func = [0]*len(funcs)
	output_bits_per_func = [0]*len(funcs)
	is_func_done = [False]*nr_of_funcs

	for i in range(n):
		val = data[i]

		buff.append(val)

		times = [t/float(freq) for t in range(len(buff))]
		for j in range(nr_of_funcs):
			if not is_func_done[j]:
				pars = grad_descent(funcs[j], times, buff)
				assert len(pars) == funcs[j].nr_of_used_pars

				max_diff = max(fabs(buff[k]-funcs[j].eval([times[k]], parameters=pars)) for k in range(len(buff)))

				if max_diff <= max_err_allowed:
					input_pnts_per_func[j] = len(buff)
					output_bits_per_func[j] = nr_of_output_bits(len(funcs), resolution, len(pars), power)
					compression_rate_of_buff_per_func[j] = resolution*len(buff) / float(output_bits_per_func[j])
				else:
					is_func_done[j] = True

		if all(is_func_done):
			best_comp_rate = max(compression_rate_of_buff_per_func)

			if best_comp_rate < compression_rate_of_1_pnt_in_1_pnt_out:
				# Output a single point.
				total_output_bits += resolution+power
				buff.pop(0)

			else:
				best_func_idx = compression_rate_of_buff_per_func.index(best_comp_rate)
				
				total_output_bits += output_bits_per_func[best_func_idx]
				prev_len = len(buff)
				buff = buff[input_pnts_per_func[best_func_idx]:]

				#print(prev_len , input_pnts_per_func[best_func_idx] , len(buff))
				assert prev_len - input_pnts_per_func[best_func_idx] == len(buff)

			compression_rate_of_buff_per_func = [0.0]*nr_of_funcs
			input_pnts_per_func = [0]*len(funcs)
			output_bits_per_func = [0]*len(funcs)
			is_func_done = [False]*nr_of_funcs

	# TODO: Flush buffer
	print("Unflushed buffer length:", len(buff))

	input_bits = resolution*n
	return input_bits / total_output_bits

if __name__ == '__main__':
	from math import sin, pi
	from cgp import CGP
	data_len = 100
	data = [sin(pi*i/10.0) for i in range(data_len)]

	from operation_table import op_table


	gene1 = [1,0,1,2]
	gene2 = [0,0,1,2]
	func1 = CGP(1, op_table, gene1, nr_of_parameters=1)
	func2 = CGP(1, op_table, gene2, nr_of_parameters=1)

	func1.print_function(['a'])
	func2.print_function(['a'])

	funcs = [func1, func2]

	print("Compression rate: ", get_compression_rate(funcs, data, 1.0e-3))

	