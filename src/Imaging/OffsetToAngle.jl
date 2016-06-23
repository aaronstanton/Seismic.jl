function OffsetToAngle(m_h;nz=1,oz=0,dz=1,nhx=1,ohx=0,dhx=1,npx=1,opx=0,dpx=1)

	nkx = 2*nhx
	dkx = 2*pi/nkx/dhx
	nf = 2*nz
	dw = 2.*pi/nf/dz
	nw = Int(floor(nf/2)) + 1
	m_h_pad = zeros(Float64,nf,nkx)
	m_h_pad[1:nz,1:nhx] = m_h
	M_h = fft(fft(m_h_pad,1),2)
	for iw = 1 : nw
		w = (iw-1)*dw
		for ikx = 1 : nkx
			kx =(ikx-1)*dkx
			L = exp(-kx*ohx*im)
			M_h[iw,ikx] *= L
		end
	end
	M_a = zeros(Complex{Float64},nf,npx)
	for iw = 1 : nw
		w = (iw-1)*dw
		for ipx = 1 : npx
			px = (ipx-1)*dpx + opx
			kx = -px*w
			if (kx >= 0)
				ikx = Int(floor(kx/dkx)) + 1
				b = (kx/dkx - floor(kx/dkx))
				a = 1.0 - b
			else
				ikx = nkx + Int(ceil(kx/dkx))
				a = abs(kx/dkx - ceil(kx/dkx))
				b = 1.0 - a
			end
			if (ikx < nkx && ikx > 0)
				M_a[iw,ipx] += a*M_h[iw,ikx] + b*M_h[iw,ikx+1]
			end
		end
	end
	# symmetries
	for iw=nw+1:nf
		M_a[iw,:] = conj(M_a[nf-iw+2,:])
	end 
	m_a = ifft(M_a,1)
	m_a = real(m_a[1:nz,:])
	return m_a
end

function AngleToOffset(m_a;nz=1,oz=0,dz=1,nhx=1,ohx=0,dhx=1,npx=1,opx=0,dpx=1)

	nkx = 4*nhx
	dkx = 2*pi/nkx/dhx
	nf = 4*nz
	dw = 2.*pi/nf/dz
	nw = Int(floor(nf/2)) + 1
	m_a_pad = zeros(Float64,nf,npx)
	m_a_pad[1:nz,1:npx] = m_a
	M_a = fft(m_a_pad,1)
	M_h = zeros(Complex{Float64},nf,nkx)
	for iw = 1 : nw
		w = (iw-1)*dw
		for ipx = 1 : npx
			px = (ipx-1)*dpx + opx
			kx = -px*w
			if (kx >= 0)
				ikx = Int(floor(kx/dkx)) + 1
				b = (kx/dkx - floor(kx/dkx))
				a = 1.0 - b
			else
				ikx = nkx + Int(ceil(kx/dkx))
				a = abs(kx/dkx - ceil(kx/dkx))
				b = 1.0 - a
			end
			if (ikx < nkx && ikx > 0)
				M_h[iw,ikx] += a*M_a[iw,ipx]
				M_h[iw,ikx+1] += b*M_a[iw,ipx]
			end
		end
	end
	for iw = 1 : nw
		w = (iw-1)*dw
		for ikx = 1 : nkx
			kx =(ikx-1)*dkx
			L = exp(-kx*ohx*im)
			M_h[iw,ikx] *= conj(L)
		end
	end
	# symmetries
	for iw=nw+1:nf
		M_h[iw,:] = conj(M_h[nf-iw+2,:])
	end 

	m_h = ifft(ifft(M_h,2),1)
	m_h = real(m_h[1:nz,1:nhx])
	return m_h
end

function OffsetToAngle(in,h::Array{Header,1};nhx=1,ohx=0,dhx=1,npx=1,opx=0,dpx=1)

	out = OffsetToAngle(in;nz=h[1].n1,oz=nz=h[1].o1,dz=h[1].d1,nhx=nhx,ohx=ohx,dhx=dhx,npx=npx,opx=opx,dpx=dpx)	
	h_out = Header[]
	for ipx = 1 : npx
		h_out = push!(h_out,h[1])  
		h_out[ipx].iang = ipx
		h_out[ipx].ang = ipx*dpx + opx
	end	
	
	return out,h_out
end

function AngleToOffset(in,h::Array{Header,1};nhx=1,ohx=0,dhx=1,npx=1,opx=0,dpx=1)

	out = AngleToOffset(in;nz=h[1].n1,oz=nz=h[1].o1,dz=h[1].d1,nhx=nhx,ohx=ohx,dhx=dhx,npx=npx,opx=opx,dpx=dpx)	
	h_out = Header[]
	for ihx = 1 : nhx
		h_out = push!(h_out,h[1])  
		h_out[ihx].ihx = ihx
		h_out[ihx].hx = ihx*dhx + ohx
	end	
	
	return out,h_out
end

function OffsetToAngle(m::ASCIIString,d::ASCIIString,adj;nhx=1,ohx=0,dhx=1,npx=1,opx=0,dpx=1)

	@compat parameters = Dict(:nhx=>nhx,:ohx=>ohx,:dhx=>dhx,:npx=>npx,:opx=>opx,:dpx=>dpx)
	if (adj==true)
		SeisProcess(d,m,[OffsetToAngle],[parameters];key=["imx"])
		f_d = ParseDataName(m)
		f_h = ParseHeaderName(m)		
		ext = ReadTextHeader(d)
		ext.n3 = npx; ext.o3 = opx; ext.d3 = dpx;
		WriteTextHeader(m,ext,"native_float",4,f_d,f_h)
	else
		SeisProcess(m,d,[AngleToOffset],[parameters];key=["imx"])
		f_d = ParseDataName(d)
		f_h = ParseHeaderName(d)		
		ext = ReadTextHeader(m)
		ext.n3 = nhx; ext.o3 = ohx; ext.d3 = dhx;
		WriteTextHeader(d,ext,"native_float",4,f_d,f_h)
	end

end
