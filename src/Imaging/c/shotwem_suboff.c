#include "seismic.h"
#include "wem_suboff.h"

void docs () {
	char * doc =
		"                                                                      \n"
		" Common shot wave equation migration for 2D acoustic, isotropic data. \n"
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
	int nx,nz,nt,ix,iz,it,nref,ihx;
	float **d=NULL,**m=NULL,**vel=NULL,**wav=NULL,sx,sz,gz;
	float ox,dx,oz,dz,ot,dt,fmin,fmax;
	int ntraces;
	int padt,padx;
	bool adj,pspi,verbose;
	struct SeisFileHeader fh;
	float kz_eps;
	int nhx;
	float ohx,dhx;

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
	if (!par_read_int(argc,argv,"nhx",&nhx)) nhx = 1;
	if (!par_read_float(argc,argv,"dhx",&dhx)) dhx = 1;
	if (!par_read_float(argc,argv,"ohx",&ohx)) ohx = 0;
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
	oz = fh.o1;
	dz = fh.d1;
	ox = h_vel[0].mx;
	dx = h_vel[1].mx - h_vel[0].mx;

	if (adj){
		h_d = allocSeisHeader(nx);
		d = alloc2float(nt,nx);
		SeisRead(d_name,d,h_d,&fh);
	}        
	else{
		h_m = allocSeisHeader(nx*nhx);
		m = alloc2float(nz,nx*nhx);
		SeisRead(m_name,m,h_m,&fh);
	}
	if (adj){
		m = alloc2float(nz,nx*nhx);
		for (iz=0;iz<nz;iz++) for (ix=0;ix<nx*nhx;ix++) m[ix][iz] = 0.0;
		h_m = allocSeisHeader(nx*nhx);
		for (ix=0;ix<nx;ix++){
			for (ihx=0;ihx<nhx;ihx++){ 
				h_m[ix*nhx + ihx] = h_vel[ix];
				h_m[ix*nhx + ihx].ihx = ihx;
				h_m[ix*nhx + ihx].hx = ihx*dhx + ohx;
			}
		}
	}
	else{  	
		h_d = allocSeisHeader(nx);
		d = alloc2float(nt,nx);
		for (it=0;it<nt;it++) for (ix=0;ix<nx;ix++) d[ix][it] = 0.0;
		for (ix=0;ix<nx;ix++) h_d[ix] = h_m[ix];
		for (ix=0;ix<nx;ix++) h_d[ix].n1 = nt;
		for (ix=0;ix<nx;ix++) h_d[ix].d1 = dt;
		for (ix=0;ix<nx;ix++) h_d[ix].o1 = ot;
		for (ix=0;ix<nx;ix++) h_d[ix].gx = h_m[ix].mx;
		for (ix=0;ix<nx;ix++) h_d[ix].hx = h_m[ix].mx - sx;
	}
	wem_suboff(d,
			m,
			wav,
			nt,ot,dt,
			nx,ox,dx,
			nhx,ohx,dhx,
			sx,
			nz,oz,dz,
			gz,sz,
			vel,nref,
			fmin,fmax,
			padt,padx,
			adj,pspi,
			verbose,
			kz_eps);
	if (adj){
		for (ix=0;ix<nx*nhx;ix++) h_m[ix].mx = h_d[ix].gx;
		for (ix=0;ix<nx*nhx;ix++) h_m[ix].imx = (int) (h_m[ix].mx - ox)/dx;
		InitFileHeader(&fh);
		fh.n1 = nz; fh.o1 = oz; fh.d1 = dz;
		fh.n2 = nhx; fh.o2 = ohx; fh.d2 = dhx;
		fh.n3 = nx; fh.o3 = ox; fh.d3 = dx;
		fh.n4 = 1; fh.o4 = 0.; fh.d4 = 1.;
		fh.n5 = 1; fh.o5 = 0.; fh.d5 = 1.;
		SeisWrite(m_name,m,h_m,&fh);
	}
	else{
		InitFileHeader(&fh);
		fh.n1 = nt; fh.o1 = ot; fh.d1 = dt;
		fh.n2 = nx; fh.o2 = ox; fh.d2 = dx;
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
