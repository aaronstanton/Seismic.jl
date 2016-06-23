function SeisProcess(in::ASCIIString,out::ASCIIString,operators,parameters;key=[])	
	# Run processing flows that read and write from disk
	#
	# f is a function that has the following syntax: d2,h2 = f(d1,h1,param), where 
	#
	# param is a dictionary (Dict) of parameters for the function.
	# note that f can be a vector of functions. They will be executed sequentially on
	# the same group of traces.
	#

	# get list of gather lengths
	ext = ReadTextHeader(in)
	nx = ext.n2*ext.n3*ext.n4*ext.n5
	filename_headers = ParseHeaderName(in)
	stream_h = open(filename_headers)
	curr = zeros(length(key),1)
	itrace = 1
	h1 = GrabHeader(stream_h,1)
	for ikey = 1 : length(key)
		curr[ikey] = getfield(h1,symbol(key[ikey]))
	end
	prev = 1*curr
	L = []
	#println("nx=",nx)
	for j = 1 : nx
		h1 = GrabHeader(stream_h,j)
		for ikey = 1 : length(key)
			curr[ikey] = getfield(h1,symbol(key[ikey]))
		end
		#println("prev=",prev,"curr=",curr,"j=",j,"itrace=",itrace)
		if curr != prev
			push!(L,j - itrace)
			itrace = j
		end
		prev = 1*curr
	end	
	push!(L,nx - itrace + 1)
	
	itrace_in = 1
	itrace_out = 1
	for igather = 1 : length(L)
		d1,h1,e1 = SeisRead(in,group="some",itrace=itrace_in,ntrace=L[igather])
		num_traces_in = size(d1,2)
		for j = 1 : length(operators)
			op = operators[j]
			d2,h2 = op(d1,h1;parameters[j]...)
			d1 = copy(d2)
			h1 = copy(h2)
		end
		num_traces_out = size(d1,2)	
		SeisWrite(out,d1,h1,e1,itrace=itrace_out)
		itrace_in += num_traces_in
		itrace_out += num_traces_out
	end

end
