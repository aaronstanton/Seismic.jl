"""
**ComputeAngles**

*Compute angles for shot gathers. These angles can be used for mapping migrated shots into angle gathers during shot profile migration.*

**IN**   


* angx = "angx" : filename for incidence angles in the x direction for each shot
* angy = "angy" : filename for incidence angles in the y direction for each shot
* dip_flag = false : flag to subtract reflector dip from the computed angles to make them with reference to reflector normal
* vel = "vel" : seis file containing the velocity (should have same x and z dimensions as the desired image)
* wav = "wav" : seis file containing the source wavelet (in time domain) 
* sz = 0. : source depth (Dev: read this from source wavelet file for variable source depth)
* nhx = 101 : number of offset bins
* ohx = 1000. : min offset (surface offset in the data)
* dhx = 10. : offset increment
* nhy = 101 : number of offset bins
* ohy = 1000. : min offset (surface offset in the data)
* dhy = 10. : offset increment
* fmin = 0. : min frequency to process (Hz)
* fmax = 80. : max frequency to process (Hz)
* padt = 2 : pad factor for the time axis
* padx = 2 : pad factor for the spatial axes
* verbose = false : flag for error / debugging messages
* sx = [0.] : array of source X positions (meters)
* sy = [0.] : array of source Y positions (meters)

**OUT**  

*Credits: AS, 2015*

"""

function ComputeAngles(angx::ASCIIString,angy::ASCIIString;vel="vel",dip_flag=false,dipx="dipx",dipy="dipy",wav="wav",sz=0.,fmin=0,fmax=80,padt=2,padx=2,verbose=false,sx=[0],sy=[0])

	nshot = length(sx)	
	v,h,e = SeisRead(vel)
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
	min_gx = h[1].mx
	max_gx = h[end].mx
	min_gy = h[1].my
	max_gy = h[end].my

	shot_list = Array(ShotAngles,nshot)
	for ishot = 1 : nshot
		shot_list[ishot] = ShotAngles(0,0,0,0,0,0,0,0,0,0,0,0,0)
		shot_list[ishot].angx = join(["angx_shot_" Int(floor(sx[ishot])) "_" Int(floor(sy[ishot]))])
		shot_list[ishot].angy = join(["angy_shot_" Int(floor(sx[ishot])) "_" Int(floor(sy[ishot]))])
		shot_list[ishot].vel = vel
		shot_list[ishot].dipx = dipx
		shot_list[ishot].dipy = dipy
		shot_list[ishot].wav = wav
		shot_list[ishot].sx = sx[ishot]
		shot_list[ishot].sy = sy[ishot]
		shot_list[ishot].sz = sz
		shot_list[ishot].fmin = fmin
		shot_list[ishot].fmax = fmax
		shot_list[ishot].dip_flag = (dip_flag == true) ? "y" : "n"
		shot_list[ishot].verbose = (verbose == true) ? "y" : "n"	
	end

	@sync @parallel for ishot = 1 : nshot
		a = compute_angles(shot_list[ishot])
	end
	extent = Seismic.Extent(nz,nx,ny,nsx,nsy,
			oz,ox,oy,osx,osy,
			dz,dx,dy,dsx,dsy,
			"Depth","mx","my","sx","sy",
			"m","m","m","m","m",
			"")
	j = 1
	for ishot = 1 : nshot
		angx_shot,h_shot = SeisRead(shot_list[ishot].angx)
		angy_shot,h_shot = SeisRead(shot_list[ishot].angy)
		SeisWrite(angx,angx_shot,h_shot,extent,itrace=j)
		SeisWrite(angy,angy_shot,h_shot,extent,itrace=j)
		j += size(angx_shot,2)
		SeisRemove(shot_list[ishot].angx)
		SeisRemove(shot_list[ishot].angy)
	end

end

type ShotAngles
	angx
	angy
	vel
	dipx
	dipy
	wav
	sx
	sy
	sz
	fmin
	fmax
	dip_flag
	verbose
end

function compute_angles(shot)
	@compat if (Libdl.find_library(["compute_angles"]) == "")
		error("Couldn't find shared library compute_angles.so, make sure it is installed and check LD_LIBRARY_PATH environmental variable.")
	end
	argv = ["0",
	join(["angx=",shot.angx]), join(["angy=",shot.angy]), join(["vel=",shot.vel]), join(["dipx=",shot.dipx]), join(["dipy=",shot.dipy]), join(["wav=",shot.wav]),  
	join(["sx=",shot.sx]), join(["sy=",shot.sy]),  join(["sz=",shot.sz]),
	join(["fmin=",shot.fmin]),  join(["fmax=",shot.fmax]), 
	join(["dip_flag=",shot.dip_flag]),  join(["verbose=",shot.verbose]) ] 
	@compat a = ccall((:main, "compute_angles"), Int32, (Int32, Ptr{Ptr{UInt8}}), length(argv), argv)                    
	return(a)

end

