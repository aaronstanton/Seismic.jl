function ShotProfileEWEM_suboff(m::Array{ASCIIString,1},d::Array{ASCIIString,1},adj=true;pspi=true,nref=5,vp="vp.seis",vs="vs.seis",wav="wav.seis",sz=0.,gz=0.,nhx=1,ohx=0,dhx=1,nhy=1,ohy=0,dhy=1,fmin=0,fmax=80,padt=1,padx=1,verbose=false,sx=[0],sy=[0],kz_eps=0.0005)

	nshot = length(sx)	
	v,h,e = SeisRead(vp)
	min_imx = h[1].imx
	max_imx = h[end].imx
	nx = max_imx - min_imx + 1
	ox = h[1].mx
	dx = h[2].mx - ox
	min_imy = h[1].imy
	max_imy = h[end].imy
	ny = max_imy - min_imy + 1
	oy = h[1].my
	dy = ny > 1 ? h[nx+1].my - oy : dx
	nz = h[1].n1
	dz = h[1].d1
	oz = h[1].o1
	w,h,e = SeisRead(wav)
	nt = h[1].n1
	dt = h[1].d1
	ot = h[1].o1
	nsx = length(unique(sx))
	osx = sx[1]	
	dsx = nsx > 1 ? sx[2] - sx[1] : 1.	
	nsy = length(unique(sy))	
	osy = sy[1]	
	dsy = nsy > 1 ? sy[2] - sy[1] : 1.
	dsx = dsx == 0. ? 1. : dsx	
	dsy = dsy == 0. ? 1. : dsy	

	shot_list = Array(Shot3C_suboff,nshot)
	for ishot = 1 : nshot
		shot_list[ishot] = Shot3C_suboff(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)
		shot_list[ishot].ux = join([d[1] "_shot_" Int(floor(sx[ishot])) "_" Int(floor(sy[ishot]))])
		shot_list[ishot].uy = join([d[2] "_shot_" Int(floor(sx[ishot])) "_" Int(floor(sy[ishot]))])
		shot_list[ishot].uz = join([d[3] "_shot_" Int(floor(sx[ishot])) "_" Int(floor(sy[ishot]))])
		shot_list[ishot].mpp = join([m[1] "_shot_" Int(floor(sx[ishot])) "_" Int(floor(sy[ishot]))])
		shot_list[ishot].mps1 = join([m[2] "_shot_" Int(floor(sx[ishot])) "_" Int(floor(sy[ishot]))])
		shot_list[ishot].mps2 = join([m[3] "_shot_" Int(floor(sx[ishot])) "_" Int(floor(sy[ishot]))])
		shot_list[ishot].vp = vp
		shot_list[ishot].vs = vs
		shot_list[ishot].wav = wav
		shot_list[ishot].sx = sx[ishot]
		shot_list[ishot].sy = sy[ishot]
		shot_list[ishot].sz = sz
		shot_list[ishot].gz = gz
		shot_list[ishot].nhx = nhx
		shot_list[ishot].nhy = nhy
		shot_list[ishot].dhx = dhx
		shot_list[ishot].dhy = dhy
		shot_list[ishot].ohx = ohx
		shot_list[ishot].ohy = ohy		
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
			SeisWindow(d[1],shot_list[ishot].ux,key=["sx","sy"],minval=[sx[ishot],sy[ishot]],maxval=[sx[ishot],sy[ishot]])
			SeisWindow(d[2],shot_list[ishot].uy,key=["sx","sy"],minval=[sx[ishot],sy[ishot]],maxval=[sx[ishot],sy[ishot]])
			SeisWindow(d[3],shot_list[ishot].uz,key=["sx","sy"],minval=[sx[ishot],sy[ishot]],maxval=[sx[ishot],sy[ishot]])
		end
		@sync @parallel for ishot = 1 : nshot
			a = shotewem_suboff(shot_list[ishot])
		end
		j = 1
		gather = zeros(Float32,nz,nhx*nhy)
		extent = Seismic.Extent(nz,nhy,nhx,ny,nx,
		oz,ohy,ohx,oy,ox,
		dz,dhy,dhx,dy,dx,
		"Depth","hy","hx","my","mx",
		"m","m","m","m","m",
		"")
		for imx = 1 : nx
			for imy = 1 : ny
				h = Header[]
				for ihx = 1 : nhx
					for ihy = 1 : nhy
						push!(h,Seismic.InitSeisHeader())		
						h[(ihx-1)*nhy + ihy].tracenum = convert(typeof(h[1].tracenum),j + (ihx-1)*nhy + ihy)
						h[(ihx-1)*nhy + ihy].o1 = convert(typeof(h[1].o1),0)
						h[(ihx-1)*nhy + ihy].n1 = convert(typeof(h[1].n1),nz)
						h[(ihx-1)*nhy + ihy].d1 = convert(typeof(h[1].d1),dz)
						h[(ihx-1)*nhy + ihy].imx = convert(typeof(h[1].imx),imx-1 + min_imx)
						h[(ihx-1)*nhy + ihy].imy = convert(typeof(h[1].imy),imy-1 + min_imy)			
						h[(ihx-1)*nhy + ihy].mx = convert(typeof(h[1].mx),dx*(imx-1 + min_imx))
						h[(ihx-1)*nhy + ihy].my = convert(typeof(h[1].my),dy*(imy-1 + min_imy))
						h[(ihx-1)*nhy + ihy].ihx = convert(typeof(h[1].ihx),ihx)		
						h[(ihx-1)*nhy + ihy].ihy = convert(typeof(h[1].ihy),ihy)		
					end
				end
				SeisWrite(m[1],gather,h,extent,itrace=j)
				SeisWrite(m[2],gather,h,extent,itrace=j)
				SeisWrite(m[3],gather,h,extent,itrace=j)
				j += nhx*nhy
			end
		end

		mpp_all,h,ext  = SeisRead(m[1])
		mps1_all,h,ext  = SeisRead(m[2])
		mps2_all,h,ext  = SeisRead(m[3])
		for ishot = 1 : nshot
			mpp_shot,h_shot,ext  = SeisRead(shot_list[ishot].mpp)
			mpp_all += mpp_shot
			mps1_shot,h_shot,ext  = SeisRead(shot_list[ishot].mps1)
			mps1_all += mps1_shot
			mps2_shot,h_shot,ext  = SeisRead(shot_list[ishot].mps2)
			mps2_all += mps2_shot
			SeisRemove(shot_list[ishot].ux)
			SeisRemove(shot_list[ishot].uy)
			SeisRemove(shot_list[ishot].uz)
			SeisRemove(shot_list[ishot].mpp)
			SeisRemove(shot_list[ishot].mps1)
			SeisRemove(shot_list[ishot].mps2)			
		end	
		SeisWrite(m[1],mpp_all,h,ext)
		SeisWrite(m[2],mps1_all,h,ext)
		SeisWrite(m[3],mps2_all,h,ext)

	else

		for ishot = 1 : nshot
			SeisCopy(m[1],shot_list[ishot].mpp)
			SeisCopy(m[2],shot_list[ishot].mps1)
			SeisCopy(m[3],shot_list[ishot].mps2)
		end
		@sync @parallel for ishot = 1 : nshot
			a = shotewem_suboff(shot_list[ishot])
		end
		extent = Seismic.Extent(nt,nsy,nsx,ny,nx,
		ot,osy,osx,oy,ox,
		dt,dsy,dsx,dy,dx,
		"Time","sy","sx","gy","gx",
		"s","m","m","m","m",
		"")
		j = 1
		for ishot = 1 : nshot
			ux_shot,h_shot,e = SeisRead(shot_list[ishot].ux)
			SeisWrite(d[1],ux_shot,h_shot,extent,itrace=j)
			uy_shot,h_shot,e = SeisRead(shot_list[ishot].uy)
			SeisWrite(d[2],uy_shot,h_shot,extent,itrace=j)
			uz_shot,h_shot,e = SeisRead(shot_list[ishot].uz)
			SeisWrite(d[1],uz_shot,h_shot,extent,itrace=j)
			j += size(ux_shot,2)
			SeisRemove(shot_list[ishot].ux)
			SeisRemove(shot_list[ishot].uy)
			SeisRemove(shot_list[ishot].uz)
			SeisRemove(shot_list[ishot].mpp)
			SeisRemove(shot_list[ishot].mps1)
			SeisRemove(shot_list[ishot].mps2)
		end	
	end
end

type Shot3C_suboff
	ux
	uy
	uz
	mpp
	mps1
	mps2
	vp
	vs
	wav
	sx
	sy
	sz
	gz
	nhx
	nhy
	dhx
	dhy
	ohx
	ohy
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

function shotewem_suboff(shot)

	@compat if (Libdl.find_library(["shotewem_suboff"]) == "")
		error("Couldn't find shared library shotewem_suboff.so, make sure it is installed and check LD_LIBRARY_PATH environmental variable.")
	end
	argv = ["0",
	join(["adj=",shot.adj]), 
	join(["pspi=",shot.pspi]), 
	join(["nref=",shot.nref]), 
	join(["ux=",shot.ux]), join(["uy=",shot.uy]), join(["uz=",shot.uz]),
	join(["mpp=",shot.mpp]), join(["mps1=",shot.mps1]), join(["mps2=",shot.mps2]),
	join(["vp=",shot.vp]), join(["vs=",shot.vs]), join(["wav=",shot.wav]), 
	join(["sx=",shot.sx]), join(["sy=",shot.sy]),  join(["sz=",shot.sz]),  join(["gz=",shot.gz]), 
	join(["nhx=",shot.nhx]), join(["nhy=",shot.nhy]), 
	join(["dhx=",shot.dhx]), join(["dhy=",shot.dhy]), 
	join(["ohx=",shot.ohx]), join(["ohy=",shot.ohy]), 
	join(["fmin=",shot.fmin]),  join(["fmax=",shot.fmax]), 
	join(["padt=",shot.padt]),  join(["padx=",shot.padx]), 
	join(["verbose=",shot.verbose]), join(["kz_eps=",shot.kz_eps]) ]
	@compat a = ccall((:main, "shotewem_suboff"), Int32, (Int32, Ptr{Ptr{UInt8}}), length(argv), argv)                    
	return(a)

end

function signf(a)

	b = a < 0.0 ? -1.0 : 1.0
	return(b)

end

