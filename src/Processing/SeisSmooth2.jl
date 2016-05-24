function SeisSmooth2(A;L1=7,L2=7)
	a1 = isodd(L1) ? 0 : 1
	a2 = isodd(L2) ? 0 : 1
	kernel = ones(eltype(A),L1,L2)/(L1*L2);
	A = conv2(A,kernel)
	A = A[floor(Int,L1/2) + 1:end-floor(Int,L1/2) + a1,floor(Int,L2/2) + 1:end-floor(Int,L2/2) + a2]
	return A
end

function SeisSmooth2(A,h::Array{Header,1};L1=7,L2=7)
	A = SeisSmooth2(A,L1=L1,L2=L2)
	return A,h
end

function SeisSmooth2(in::ASCIIString,out::ASCIIString;key=["iaz";"iang"],L1=7,L2=7)
	@compat parameters = Dict(:L1=>L1,:L2=>L2)
	SeisProcess(in,out,[SeisSmooth2],[parameters];group="gather",key=key)
end
