function differentiate_traces(d;dt=0.001,pow=-2)

	nt = size(d,1)
	D = fft(d,1)
	dw = 2.*pi/nt/dt
	nw = Int(nt/2) + 1
	eps = pow < 0 ? 10.0 : 0.0
	for iw=1:nw
		D[iw,:] *= ((iw*dw) + eps)^pow
	end
	# symmetries
	for iw=nw+1:nt
		D[iw,:] = conj(D[nt-iw+2,:])
	end
	d = real(ifft(D,1))
	return d

end

function SeisDiff(in,h::Array{Header,1};pow=-2)

	out = differentiate_traces(in;dt=h[1].d1,pow=pow)
	return out,h

end

function SeisDiff(in::ASCIIString,out::ASCIIString;pow=-2)

	@compat parameters = Dict(:pow=>pow)
	SeisProcess(in,out,[SeisDiff],[parameters];group="some")

end

