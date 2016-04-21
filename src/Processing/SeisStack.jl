function stacktraces(d,h;normalize=true)
	h_out = [h[1]]
	h_out[1].trid = size(d,2)
	if (normalize == true)
		val = sum(d,2)/size(d,2)
	else
		val = sum(d,2)
	end
	return val, h_out
end

function SeisStack(in::ASCIIString,out::ASCIIString;key=["imx" "imy"],normalize=true)

	@compat parameters = Dict(:normalize=>true)
	SeisProcess(in,out,[stacktraces],[parameters];group="gather",key=key)
	DATAPATH = get(ENV,"DATAPATH",join([pwd(),"/"]))
	filename_data_out = join([DATAPATH out "@data@"])
	filename_headers_out = join([DATAPATH out "@headers@"])
	extent = ReadTextHeader(in)
	extent.n4 = 1
	extent.n5 = 1
	WriteTextHeader(out,extent,"native_float",4,filename_data_out,filename_headers_out)

end
