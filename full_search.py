from math import inf

from tune_compressor import neg_fitness

def get_all_possible_node_parts(prev_nodes, op_table):
	node_parts = []
	for op_nr in range(len(op_table)):
		if op_table[op_nr].is_binary:
			for d1 in range(prev_nodes):
				for d2 in range(prev_nodes):
					node_parts.append([op_nr, d1, d2])
		else:
			for d in range(prev_nodes):
				node_parts.append([op_nr, d, 0])
	return node_parts

def update_pos(pos, all_node_parts):
	pos[0] += 1
	for i in range(len(pos)):
		if pos[i] == len(all_node_parts[i]):
			if i == len(pos)-1:
				return False
			else:
				pos[i] = 0
				pos[i+1] += 1
		else:
			break
	return pos

def get_all_genes_sub(dims, nr_of_used_nodes, op_table):
	genes_sub = []
	#print("all node parts start")
	all_node_parts = [get_all_possible_node_parts(dims+i, op_table) for i in range(nr_of_used_nodes+1)]
	#print("all node parts done")
	pos = [0]*(nr_of_used_nodes+1)

	old_last_pos = pos[-1]
	while True:
		gene = []

		for i in range(nr_of_used_nodes+1):
			p = pos[i]
			assert len(all_node_parts[i][p]) == 3
			gene += all_node_parts[i][p]
		genes_sub.append(gene)

		# Update pos
		pos = update_pos(pos, all_node_parts)

		if pos is False:
			break

		if pos[-1] != old_last_pos:
			#print("here (in):", pos[-1], len(all_node_parts[-1]))
			old_last_pos = pos[-1]

	return genes_sub

def nodes_used_rec(nodes_used, gene, current_node, op_table, dims):
	op = op_table[gene[current_node*3]]

	if op.is_binary:
		node_1 = gene[current_node*3+1] - dims
		if node_1 >= 0:
			nodes_used[node_1] = True
			nodes_used_rec(nodes_used, gene, node_1, op_table, dims)

		node_2 = gene[current_node*3+2] - dims
		if node_2 >= 0:
			nodes_used[node_2] = True
			nodes_used_rec(nodes_used, gene, node_2, op_table, dims)
	else:
		node = gene[current_node*3+1] - dims
		if node >= 0:
			nodes_used[node] = True
			nodes_used_rec(nodes_used, gene, node, op_table, dims)

def all_nodes_are_used(gene, nr_of_nodes, dims, op_table):
	nodes_used = [False for _ in range(nr_of_nodes)]
	first_node = gene[-1]-dims

	if first_node >= 0:
		nodes_used[first_node] = True
		nodes_used_rec(nodes_used, gene, first_node, op_table, dims)
		return all(nodes_used)
	else:
		return False

def get_all_genes(dims, nr_of_nodes, error_func, nr_of_parameters_in_inp, op_table):
	from operation_table import op_table
	genes = []

	gene_len = nr_of_nodes*3+1

	best_err = inf
	best_gene = None

	# First we create genes that take the one of the inputs as output
	for d in range(dims):
		gene = [0]*gene_len
		gene[-1] = d

		if all_nodes_are_used(gene, nr_of_nodes, dims, op_table):
			err = error_func(gene)
			if err < best_err:
				best_err = err
				best_gene = list(gene)
				print("best", err)

	# And then the other ones
	for node in range(nr_of_nodes):
		#print("Start")
		subs = get_all_genes_sub(dims, node, op_table)
		#print("mid")
		len_of_sub = len(subs[0])
		remaining_len = gene_len-len_of_sub
		assert remaining_len >= 0
		if remaining_len != 0:
			for i in range(len(subs)):
				subs[i] = subs[i] + [0]*remaining_len
				subs[i][-1] = dims + node

				if all_nodes_are_used(subs[i], nr_of_nodes, dims, op_table):
					err = error_func(subs[i])
					if err < best_err:
						best_err = err
						best_gene = list(subs[i])
						print("best", err, "   progress:", float(i)/len(subs))
				subs[i] = None
				if i%50000 == 0:
					print("here", node, nr_of_nodes)
					print("itr", float(i)/len(subs))
			subs = None
	return best_gene, best_err

def full_search_part(data, nr_of_parameters, dims, nr_of_nodes, error_func):
	# TODO: Is this function really needed?
	
	# Add some default arguments to the error function
	from operation_table import op_table
	nr_of_funcs = len(op_table)

	return get_all_genes(dims, nr_of_nodes, error_func, nr_of_parameters, op_table)

def full_search(max_nr_of_nodes, data, nr_of_parameters, error_func, dims=1):
	best_err = inf
	best_gene = None
	total_dims = nr_of_parameters+dims
	for nr_of_nodes in range(1, max_nr_of_nodes+1):
		print("Nr of nodes:", nr_of_nodes)
		gene, err = full_search_part(data, nr_of_parameters, total_dims, nr_of_nodes, error_func)

		if err < best_err:
			best_err = err
			best_gene = gene
	assert best_gene != None
	return best_gene, best_err

