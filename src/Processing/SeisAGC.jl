#simple AGC with a window specified in samples
#written by GRAM 
#modified by aaron

function SeisAGC(in;L=250)

	nt = size(in,1)
	nx = size(in[:,:],2)
	out = zeros(eltype(in),size(in))
	for ix = 1 : nx
		for it = 1 : L
			rms = sqrt(sum(in[ 1 : it + L ,ix].^2))
			out[it,ix] = in[it,ix] / (rms + 1.e-20)
		end
		for it = L + 1 : nt - L
			rms = sqrt(sum(in[it - L : it + L,ix].^2))
			out[it,ix] = in[it,ix] / (rms + 1.e-20)
		end
		for it = nt - L + 1 : nt
			rms = sqrt(sum(in[ it - L : nt ,ix].^2))
			out[it,ix] = in[it,ix] / (rms + 1.e-20)
		end
	end
	return out


end

function SeisAGC(in,h::Array{Header,1};L=250)

	out = SeisAGC(in;L=L)
	return out,h

end

function SeisAGC(in::ASCIIString,out::ASCIIString;L=250)

	@compat parameters = Dict(:L=>L)
	SeisProcess(in,out,[SeisAGC],[parameters];group="some")

end

