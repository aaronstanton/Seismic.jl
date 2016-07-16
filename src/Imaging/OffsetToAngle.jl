function OffsetToAngle(m_h;nz=1,oz=0,dz=1,nhx=1,ohx=0,dhx=1,npx=1,opx=0,dpx=1)

	nkx = 4*nhx
	dkx = 2.0*pi/nkx/dhx
	nkx_half = Int(floor(nkx/2)) + 1
	nkz = 2*nz
	dkz = 2.0*pi/nkz/dz
	nkz_half = Int(floor(nkz/2)) + 1
	m_h_pad = zeros(Float64,nkz,nkx)
	m_h_pad[1:nz,1:nhx] = m_h
	M_h = fftshift(fft(m_h_pad),2)/sqrt(nkz*nkx)
	M_a = zeros(Complex{Float64},nkz,npx)
	okx = -pi/dhx
	for ikz = 1 : nkz
		if (ikz <= nkz_half)
			kz = (ikz-1)*dkz
		else
			kz = -(dkz*nkz -  (ikz-1)*dkz)
		end
		for ipx = 1 : npx
			px = (ipx-1)*dpx + opx
			kx = -px*kz
			ikx = Int(floor((kx - okx)/dkx)) + 1
			if (ikx < nkx && ikx > 0)
				b = (kx/dkx - floor(kx/dkx))
				a = 1.0 - b
				M_a[ikz,ipx] += (a*M_h[ikz,ikx] + b*M_h[ikz,ikx+1])*exp(-kx*ohx*im)
			end
		end
	end
	m_a = bfft(M_a,1)/sqrt(nkz)
	m_a = real(m_a[1:nz,:])
	return m_a
end

function AngleToOffset(m_a;nz=1,oz=0,dz=1,nhx=1,ohx=0,dhx=1,npx=1,opx=0,dpx=1)

	nkx = 4*nhx
	dkx = 2.0*pi/nkx/dhx
	nkx_half = Int(floor(nkx/2)) + 1
	nkz = 2*nz
	dkz = 2.0*pi/nkz/dz
	nkz_half = Int(floor(nkz/2)) + 1
	m_a_pad = zeros(Float64,nkz,npx)
	m_a_pad[1:nz,1:npx] = m_a
	M_a = fft(m_a_pad,1)/sqrt(nkz)
	M_h = zeros(Complex{Float64},nkz,nkx)
	okx = -pi/dhx
	for ikz = 1 : nkz
		if (ikz <= nkz_half)
			kz = (ikz-1)*dkz
		else
			kz = -(dkz*nkz -  (ikz-1)*dkz)
		end
		for ipx = 1 : npx
			px = (ipx-1)*dpx + opx
			kx = -px*kz
			ikx = Int(floor((kx - okx)/dkx)) + 1
			if (ikx < nkx && ikx > 0)
				b = (kx/dkx - floor(kx/dkx))
				a = 1.0 - b
				M_h[ikz,ikx] += a*M_a[ikz,ipx]*exp(kx*ohx*im)
				M_h[ikz,ikx+1] += b*M_a[ikz,ipx]*exp(kx*ohx*im)
			end
		end
	end
	m_h = bfft(ifftshift(M_h,2))/sqrt(nkz*nkx)
	m_h = real(m_h[1:nz,1:nhx])
	return m_h
end

function OffsetToAngle(in,h::Array{Header,1};nhx=1,ohx=0,dhx=1,npx=1,opx=0,dpx=1)

	out = OffsetToAngle(in;nz=h[1].n1,oz=h[1].o1,dz=h[1].d1,nhx=nhx,ohx=ohx,dhx=dhx,npx=npx,opx=opx,dpx=dpx)	
	h_out = Header[]
	for ipx = 1 : npx
		h_out = push!(h_out,h[1])  
		h_out[ipx].iang = ipx
		h_out[ipx].ang = ipx*dpx + opx
	end	

	return out,h_out
end

function AngleToOffset(in,h::Array{Header,1};nhx=1,ohx=0,dhx=1,npx=1,opx=0,dpx=1)

	out = AngleToOffset(in;nz=h[1].n1,oz=h[1].o1,dz=h[1].d1,nhx=nhx,ohx=ohx,dhx=dhx,npx=npx,opx=opx,dpx=dpx)	
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
		ext.n2 = npx; ext.o2 = opx; ext.d2 = dpx;
		WriteTextHeader(m,ext,"native_float",4,f_d,f_h)
	else
		SeisProcess(m,d,[AngleToOffset],[parameters];key=["imx"])
		f_d = ParseDataName(d)
		f_h = ParseHeaderName(d)		
		ext = ReadTextHeader(m)
		ext.n2 = nhx; ext.o2 = ohx; ext.d2 = dhx;
		WriteTextHeader(d,ext,"native_float",4,f_d,f_h)
	end

end

function OffsetToAngle(m::Array{ASCIIString,1},d::Array{ASCIIString,1},adj;nhx=1,ohx=0,dhx=1,npx=1,opx=0,dpx=1)

	for j = 1 : length(m)
		OffsetToAngle(m[j],d[j],adj;nhx=nhx,ohx=ohx,dhx=dhx,npx=npx,opx=opx,dpx=dpx)
	end     

end

