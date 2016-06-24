function CalculateSampling(in,h;cutoff=1e-10)

	itrace = 1
	wd = zeros(Float32,size(in))
	for itrace = 1 : size(in[:,:],2)
		a = sqrt(sum(in[:,itrace].^2))
		if (a > cutoff) 
			wd[:,itrace] = 1.
		end
	end
	return wd,h;
end

function CalculateSampling(in::ASCIIString,wd::ASCIIString;cutoff=1e-10)
	# calculate sampling operator (1's for live traces, 0's for missing traces)

	@compat parameters = Dict(:cutoff=>cutoff)
	SeisProcess(in,wd,[CalculateSampling],[parameters])
		
end

