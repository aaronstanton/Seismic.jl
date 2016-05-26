#include "seismic.h"
#include "illum.h"

void docs () {
	char * doc =
		"                                                                      \n"
		" Common shot wave equation source illumination calculation for 3D acoustic, isotropic data. \n"
		"                                                                      \n";
	fprintf (stderr, "%s", doc);
}

int main (int argc, char *argv[])
{
	char m_name[512],vel_name[512],wav_name[512];
	struct SeisHeader *h_m=NULL;
	struct SeisHeader *h_vel=NULL;
	struct SeisHeader *h_wav=NULL;
	int nx,ny,nz,nt,ix,iz,nref;
	float **m=NULL,**vel=NULL,**wav=NULL,sx,sy,sz;
	float ox,dx,oy,dy,oz,dz,ot,dt,fmin,fmax;
	int ntraces;
	int padt,padx;
	bool pspi,verbose;
	struct SeisFileHeader fh;
	float kz_eps;
	
	if (!par_read_bool(argc,argv,"pspi",&pspi)) pspi = true;
	if (!par_read_bool(argc,argv,"verbose",&verbose)) verbose = false;
	if (!par_read_string(argc,argv,"m", m_name)) { docs (); exit (1); }
	if (!par_read_string(argc,argv,"vel", vel_name)) { docs (); exit (1); }
	if (!par_read_string(argc,argv,"wav", wav_name)) { docs (); exit (1); }
	if (!par_read_float(argc,argv,"sz",&sz)) sz = 0;
	if (!par_read_float(argc,argv,"sx",&sx)) sx = 0;
	if (!par_read_float(argc,argv,"sy",&sy)) sy = 0;	
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

	m = alloc2float(nz,nx*ny);
	for (iz=0;iz<nz;iz++) for (ix=0;ix<nx*ny;ix++) m[ix][iz] = 0.0;
	h_m = allocSeisHeader(nx*ny);
	for (ix=0;ix<nx*ny;ix++) h_m[ix] = h_vel[ix];	       
	illum(m,
	    wav,
	    nt,ot,dt,
	    nx,ox,dx,
	    ny,oy,dy,
	    sx,sy,
	    nz,oz,dz,
	    sz,
	    vel,nref,
	    fmin,fmax,
	    padt,padx,
	    pspi,
	    verbose,
	    kz_eps);
	for (ix=0;ix<nx*ny;ix++) h_m[ix].mx = ox + ix*dx;
	for (ix=0;ix<nx*ny;ix++) h_m[ix].my = oy;
	for (ix=0;ix<nx*ny;ix++) h_m[ix].imx = (int) (h_m[ix].mx - ox)/dx;
	for (ix=0;ix<nx*ny;ix++) h_m[ix].imy = (int) (h_m[ix].my - oy)/dy;
	InitFileHeader(&fh);
	fh.n1 = nz; fh.o1 = oz; fh.d1 = dz;
	fh.n2 = nx; fh.o2 = ox; fh.d2 = dx;
	fh.n3 = ny; fh.o3 = oy; fh.d3 = dy;
	SeisWrite(m_name,m,h_m,&fh);
	free2float(m); 
	free2float(vel); 
	free2float(wav); 
	freeSeisHeader(h_m);
	freeSeisHeader(h_vel);
	freeSeisHeader(h_wav);
	return 0;
}
