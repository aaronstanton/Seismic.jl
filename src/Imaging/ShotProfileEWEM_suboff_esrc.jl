function ShotProfileEWEM_suboff_esrc(m::Array{ASCIIString,1},d::Array{ASCIIString,1},adj=true;pspi=true,nref=5,vp="vp.seis",vs="vs.seis",wav_p="wav_p.seis",wav_s="wav_s.seis",sz=0.,gz=0.,nhx=1,ohx=0,dhx=1,fmin=0,fmax=80,padt=1,padx=1,verbose=false,sx=[0],kz_eps=0.0005)

	nshot = length(sx)	
	v,h,e = SeisRead(vp)
	min_imx = h[1].imx
	max_imx = h[end].imx
	nx = max_imx - min_imx + 1
	ox = h[1].mx
	dx = h[2].mx - ox
	nz = h[1].n1
	dz = h[1].d1
	oz = h[1].o1
	w_p,h,e = SeisRead(wav_p)
	nt = h[1].n1
	dt = h[1].d1
	ot = h[1].o1
	nsx = length(unique(sx))
	osx = sx[1]	
	dsx = nsx > 1 ? sx[2] - sx[1] : 1.	
	dsx = dsx == 0. ? 1. : dsx	

	shot_list = Array(Shot2C_suboff_esrc,nshot)
	for ishot = 1 : nshot
		shot_list[ishot] = Shot2C_suboff_esrc(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
		shot_list[ishot].ux = join([d[1] "_shot_" Int(floor(sx[ishot])) ])
		shot_list[ishot].uz = join([d[2] "_shot_" Int(floor(sx[ishot])) ])
		shot_list[ishot].mpp = join([m[1] "_shot_" Int(floor(sx[ishot]))])
		shot_list[ishot].mps = join([m[2] "_shot_" Int(floor(sx[ishot]))])
		shot_list[ishot].vp = vp
		shot_list[ishot].vs = vs
		shot_list[ishot].wav_p = wav_p
		shot_list[ishot].wav_s = wav_s
		shot_list[ishot].sx = sx[ishot]
		shot_list[ishot].sz = sz
		shot_list[ishot].gz = gz
		shot_list[ishot].nhx = nhx
		shot_list[ishot].dhx = dhx
		shot_list[ishot].ohx = ohx
		shot_list[ishot].fmin = fmin
		shot_list[ishot].fmax = fmax
		shot_list[ishot].padt = padt
		shot_list[ishot].padx = padx		
		shot_list[ishot].adj = (adj == true) ? "y" : "n"
		shot_list[ishot].pspi = (pspi == true) ? "y" : "n"
		shot_list[ishot].nref = nref
		shot_list[ishot].verbose = (verbose == true) ? "y" : "n"
		shot_list[ishot].kz_eps = kz_eps
	end

	if (adj == true)
		for ishot = 1 : nshot
			SeisWindow(d[1],shot_list[ishot].ux,key=["sx"],minval=[sx[ishot]],maxval=[sx[ishot]])
			SeisWindow(d[2],shot_list[ishot].uz,key=["sx"],minval=[sx[ishot]],maxval=[sx[ishot]])
		end
		@sync @parallel for ishot = 1 : nshot
			a = shotewem_suboff_esrc(shot_list[ishot])
		end
		j = 1
		gather = zeros(Float32,nz,nhx)
		extent = Seismic.Extent(nz,nhx,nx,1,1,
		oz,ohx,ox,0,0,
		dz,dhx,dx,1,1,
		"Depth","hx","mx","","",
		"m","m","m","","",
		"")
		for imx = 1 : nx
			h = Header[]
			for ihx = 1 : nhx
				push!(h,Seismic.InitSeisHeader())		
				h[ihx].tracenum = convert(typeof(h[1].tracenum),j + ihx)
				h[ihx].o1 = convert(typeof(h[1].o1),0)
				h[ihx].n1 = convert(typeof(h[1].n1),nz)
				h[ihx].d1 = convert(typeof(h[1].d1),dz)
				h[ihx].imx = convert(typeof(h[1].imx),imx-1 + min_imx)
				h[ihx].mx = convert(typeof(h[1].mx),dx*(imx-1 + min_imx))
				h[ihx].ihx = convert(typeof(h[1].ihx),ihx)		
			end
			SeisWrite(m[1],gather,h,extent,itrace=j)
			SeisWrite(m[2],gather,h,extent,itrace=j)
			j += nhx
		end

		mpp_all,h,ext  = SeisRead(m[1])
		mps_all,h,ext  = SeisRead(m[2])
		for ishot = 1 : nshot
			mpp_shot,h_shot,ext  = SeisRead(shot_list[ishot].mpp)



			offset = Seismic.ExtractHeader(h_shot,"mx") - shot_list[ishot].sx
			mpp_shot = SeisMute(mpp_shot,offset=offset,tmute=0.0,vmute=0.35,taper=1000,dt=h_shot[1].d1)

			mpp_all[:] = mpp_all[:] + mpp_shot[:]
			mps_shot,h_shot,ext  = SeisRead(shot_list[ishot].mps)

			offset = Seismic.ExtractHeader(h_shot,"mx") - shot_list[ishot].sx
			mps_shot = SeisMute(mps_shot,offset=offset,tmute=0.0,vmute=0.35,taper=1000,dt=h_shot[1].d1)

			mps_all[:] = mps_all[:] + mps_shot[:]
			SeisRemove(shot_list[ishot].ux)
			SeisRemove(shot_list[ishot].uz)
			SeisRemove(shot_list[ishot].mpp)
			SeisRemove(shot_list[ishot].mps)
		end	
		SeisWrite(m[1],mpp_all,h,ext)
		SeisWrite(m[2],mps_all,h,ext)

	else

		for ishot = 1 : nshot
			SeisCopy(m[1],shot_list[ishot].mpp)
			SeisCopy(m[2],shot_list[ishot].mps)
		end
		@sync @parallel for ishot = 1 : nshot
			a = shotewem_suboff_esrc(shot_list[ishot])
		end
		extent = Seismic.Extent(nt,nx,nsx,1,1,
		ot,ox,osx,0,0,
		dt,dx,dsx,1,1,
		"Time","gx","sx","","",
		"s","m","m","m","m",
		"")
		j = 1
		for ishot = 1 : nshot
			ux_shot,h_shot,e = SeisRead(shot_list[ishot].ux)
			SeisWrite(d[1],ux_shot,h_shot,extent,itrace=j)
			uz_shot,h_shot,e = SeisRead(shot_list[ishot].uz)
			SeisWrite(d[2],uz_shot,h_shot,extent,itrace=j)
			j += size(ux_shot,2)
			SeisRemove(shot_list[ishot].ux)
			SeisRemove(shot_list[ishot].uz)
			SeisRemove(shot_list[ishot].mpp)
			SeisRemove(shot_list[ishot].mps)
		end	
	end
end

type Shot2C_suboff_esrc
	ux
	uz
	mpp
	mps
	vp
	vs
	wav_p
	wav_s
	sx
	sz
	gz
	nhx
	dhx
	ohx
	fmin
	fmax
	padt
	padx
	adj
	pspi
	nref
	verbose
	kz_eps
end

function shotewem_suboff_esrc(shot)

	@compat if (Libdl.find_library(["shotewem_suboff_esrc"]) == "")
		error("Couldn't find shared library shotewem_suboff_esrc.so, make sure it is installed and check LD_LIBRARY_PATH environmental variable.")
	end
	argv = ["0",
	join(["adj=",shot.adj]), 
	join(["pspi=",shot.pspi]), 
	join(["nref=",shot.nref]), 
	join(["ux=",shot.ux]), join(["uz=",shot.uz]),
	join(["mpp=",shot.mpp]), join(["mps=",shot.mps]),
	join(["vp=",shot.vp]), join(["vs=",shot.vs]), join(["wav_p=",shot.wav_p]), join(["wav_s=",shot.wav_s]), 
	join(["sx=",shot.sx]), join(["sz=",shot.sz]),  join(["gz=",shot.gz]), 
	join(["nhx=",shot.nhx]),
	join(["dhx=",shot.dhx]),
	join(["ohx=",shot.ohx]),
	join(["fmin=",shot.fmin]),  join(["fmax=",shot.fmax]), 
	join(["padt=",shot.padt]),  join(["padx=",shot.padx]), 
	join(["verbose=",shot.verbose]), join(["kz_eps=",shot.kz_eps]) ]
	@compat a = ccall((:main, "shotewem_suboff_esrc"), Int32, (Int32, Ptr{Ptr{UInt8}}), length(argv), argv)                    
	return(a)

end

function signf(a)

	b = a < 0.0 ? -1.0 : 1.0
	return(b)

end

