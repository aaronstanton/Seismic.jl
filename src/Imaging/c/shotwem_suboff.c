#include "seismic.h"
#include "wem_suboff.h"

void docs () {
	char * doc =
		"                                                                      \n"
		" Common shot wave equation migration for 3D acoustic, isotropic data. \n"
		"                                                                      \n";
	fprintf (stderr, "%s", doc);
}

int main (int argc, char *argv[])
{
	char d_name[512],m_name[512],vel_name[512],wav_name[512];
	struct SeisHeader *h_d=NULL;
	struct SeisHeader *h_m=NULL;
	struct SeisHeader *h_vel=NULL;
	struct SeisHeader *h_wav=NULL;
	int nx,ny,nz,nt,ix,iy,iz,it,nref,ihx,ihy;
	float **d=NULL,**m=NULL,**vel=NULL,**wav=NULL,sx,sy,sz,gz;
	float ox,dx,oy,dy,oz,dz,ot,dt,fmin,fmax;
	int ntraces;
	int padt,padx;
	bool adj,pspi,verbose;
	struct SeisFileHeader fh;
	float kz_eps;
	int nhx,nhy;
	float ohx,ohy,dhx,dhy;
	
	if (!par_read_bool(argc,argv,"adj",&adj)) adj = true;
	if (!par_read_bool(argc,argv,"pspi",&pspi)) pspi = true;
	if (!par_read_bool(argc,argv,"verbose",&verbose)) verbose = false;
	if (!par_read_string(argc,argv,"d", d_name)) { docs (); exit (1); }
	if (!par_read_string(argc,argv,"m", m_name)) { docs (); exit (1); }
	if (!par_read_string(argc,argv,"vel", vel_name)) { docs (); exit (1); }
	if (!par_read_string(argc,argv,"wav", wav_name)) { docs (); exit (1); }
	if (!par_read_float(argc,argv,"sz",&sz)) sz = 0;
	if (!par_read_float(argc,argv,"gz",&gz)) gz = 0;
	if (!par_read_float(argc,argv,"sx",&sx)) sx = 0;
	if (!par_read_float(argc,argv,"sy",&sy)) sy = 0;
	if (!par_read_int(argc,argv,"nhx",&nhx)) nhx = 1;
	if (!par_read_float(argc,argv,"dhx",&dhx)) dhx = 1;
	if (!par_read_float(argc,argv,"ohx",&ohx)) ohx = 0;
	if (!par_read_int(argc,argv,"nhy",&nhy)) nhy = 1;
	if (!par_read_float(argc,argv,"dhy",&dhy)) dhy = 1;
	if (!par_read_float(argc,argv,"ohy",&ohy)) ohy = 0;
	if (!par_read_float(argc,argv,"fmin",&fmin)) fmin = 0;
	if (!par_read_float(argc,argv,"fmax",&fmax)) fmax = 80;
	if (!par_read_int(argc,argv,"padt",&padt)) padt = 1;
	if (!par_read_int(argc,argv,"padx",&padx)) padx = 1;
	if (!par_read_int(argc,argv,"nref",&nref)) nref = 5;
	if (!par_read_float(argc,argv,"kz_eps",&kz_eps)) kz_eps = 0.0005;
	// get dimensions from velocity (nz,oz,dz,nx,ox,dx) and wavelet (nt,sx) files
	InitFileHeader(&fh);
	ReadFileHeader(wav_name,&fh);
	nt = fh.n1;
	ReadFileHeader(vel_name,&fh);
	nz = fh.n1; ntraces = fh.n2*fh.n3*fh.n4*fh.n5;
	
	//fprintf(stderr,"nx=%d\n",nx);
	//fprintf(stderr,"nz=%d\n",nz);
	//fprintf(stderr,"nt=%d\n",nt);
	h_wav = allocSeisHeader(1);
	wav = alloc2float(nt,1); 		
	SeisRead(wav_name,wav,h_wav,&fh);
	ot = fh.o1;
	dt = fh.d1;
	h_vel = allocSeisHeader(ntraces);
	vel = alloc2float(nz,ntraces);
	SeisRead(vel_name,vel,h_vel,&fh); 	
	nx = h_vel[ntraces-1].imx - h_vel[0].imx + 1;
	ny = h_vel[ntraces-1].imy - h_vel[0].imy + 1;
	oz = fh.o1;
	dz = fh.d1;
	ox = h_vel[0].mx;
	dx = h_vel[1].mx - h_vel[0].mx;
	oy = h_vel[0].my;
	dy = dx;

	if (adj){
		h_d = allocSeisHeader(nx*ny);
		d = alloc2float(nt,nx*ny);
		SeisRead(d_name,d,h_d,&fh);
	}        
	else{
		h_m = allocSeisHeader(nx*ny*nhx*nhy);
		m = alloc2float(nz,nx*ny*nhx*nhy);
		SeisRead(m_name,m,h_m,&fh);
	}
	if (adj){
		m = alloc2float(nz,nx*ny*nhx*nhy);
		for (iz=0;iz<nz;iz++) for (ix=0;ix<nx*ny*nhx*nhy;ix++) m[ix][iz] = 0.0;
		h_m = allocSeisHeader(nx*ny*nhx*nhy);
		for (ix=0;ix<nx;ix++){
			for (iy=0;iy<ny;iy++){
				for (ihx=0;ihx<nhx;ihx++){ 
					for (ihy=0;ihy<nhy;ihy++){ 
						h_m[ix*ny*nhx*nhy + iy*nhx*nhy + ihx*nhy + ihy] = h_vel[ix*ny + iy];
						h_m[ix*ny*nhx*nhy + iy*nhx*nhy + ihx*nhy + ihy].ihx = ihx;
						h_m[ix*ny*nhx*nhy + iy*nhx*nhy + ihx*nhy + ihy].ihy = ihy;
						h_m[ix*ny*nhx*nhy + iy*nhx*nhy + ihx*nhy + ihy].hx = ihx*dhx + ohx;
						h_m[ix*ny*nhx*nhy + iy*nhx*nhy + ihx*nhy + ihy].hy = ihy*dhy + ohy;
					}
				}
			}
		}
	}
	else{  	
		h_d = allocSeisHeader(nx*ny);
		d = alloc2float(nt,nx*ny);
		for (it=0;it<nt;it++) for (ix=0;ix<nx*ny;ix++) d[ix][it] = 0.0;
		for (ix=0;ix<nx*ny;ix++) h_d[ix] = h_m[ix];
		for (ix=0;ix<nx*ny;ix++) h_d[ix].n1 = nt;
		for (ix=0;ix<nx*ny;ix++) h_d[ix].d1 = dt;
		for (ix=0;ix<nx*ny;ix++) h_d[ix].o1 = ot;
		for (ix=0;ix<nx*ny;ix++) h_d[ix].gx = h_m[ix].mx;
		for (ix=0;ix<nx*ny;ix++) h_d[ix].gy = h_m[ix].my;
		for (ix=0;ix<nx*ny;ix++) h_d[ix].hx = h_m[ix].mx - sx;
		for (ix=0;ix<nx*ny;ix++) h_d[ix].hy = h_m[ix].my - sy;  		
	}
	wem_suboff(d,
	    m,
	    wav,
	    nt,ot,dt,
	    nx,ox,dx,
	    ny,oy,dy,
	    nhx,ohx,dhx,
	    nhy,ohy,dhy,
	    sx,sy,
	    nz,oz,dz,
	    gz,sz,
	    vel,nref,
	    fmin,fmax,
	    padt,padx,
	    adj,pspi,
	    verbose,
	    kz_eps);
	if (adj){
		for (ix=0;ix<nx*ny*nhx*nhy;ix++) h_m[ix].mx = h_d[ix].gx;
		for (ix=0;ix<nx*ny*nhx*nhy;ix++) h_m[ix].my = h_d[ix].gy;
		for (ix=0;ix<nx*ny*nhx*nhy;ix++) h_m[ix].imx = (int) (h_m[ix].mx - ox)/dx;
		for (ix=0;ix<nx*ny*nhx*nhy;ix++) h_m[ix].imy = (int) (h_m[ix].my - oy)/dy;
		InitFileHeader(&fh);
		fh.n1 = nz; fh.o1 = oz; fh.d1 = dz;
		fh.n2 = nhy; fh.o2 = ohy; fh.d2 = dhy;
		fh.n3 = nhx; fh.o3 = ohx; fh.d3 = dhx;
		fh.n4 = ny; fh.o4 = oy; fh.d4 = dy;
		fh.n5 = nx; fh.o5 = ox; fh.d5 = dx;
		SeisWrite(m_name,m,h_m,&fh);
	}
	else{
		InitFileHeader(&fh);
		fh.n1 = nt; fh.o1 = ot; fh.d1 = dt;
		fh.n2 = ny; fh.o2 = oy; fh.d2 = dy;
		fh.n3 = nx; fh.o3 = ox; fh.d3 = dx;
		SeisWrite(d_name,d,h_d,&fh);
	}
	free2float(d); 
	free2float(m); 
	free2float(vel); 
	free2float(wav); 
	freeSeisHeader(h_d);
	freeSeisHeader(h_m);
	freeSeisHeader(h_vel);
	freeSeisHeader(h_wav);
	return 0;
}
