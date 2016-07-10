#include "seismic.h"
#include "wem_suboff.h"

void docs () {
	char * doc =
		"                                                                      \n"
		" Common shot wave equation migration for 2D-1C acoustic, isotropic    \n"
		" data.                                                                \n"
		"                                                                      \n";
	fprintf (stderr, "%s", doc);
}

int main (int argc, char *argv[])
{
	char uz_name[512];
	char mpp_name[512];
	char vp_name[512],wav_name[512];
	struct SeisHeader *h_uz=NULL;
	struct SeisHeader *h_mpp=NULL;
	struct SeisHeader *h_vp=NULL;
	struct SeisHeader *h_wav=NULL;
	int nx,nz,nt,ix,iz,it,nref,ihx;
	float **uz=NULL,**mpp=NULL,**vp=NULL,**wav=NULL,sx,sz,gz;
	float ox,dx,oz,dz,ot,dt,fmin,fmax;
	int ntraces;
	int padt,padx;
	bool adj,pspi,verbose;
	struct SeisFileHeader fh;
	int Lw;
	float kz_eps;
	int nhx;
	float ohx,dhx;

	if (!par_read_bool(argc,argv,"adj",&adj)) adj = true;
	if (!par_read_bool(argc,argv,"pspi",&pspi)) pspi = true;
	if (!par_read_bool(argc,argv,"verbose",&verbose)) verbose = false;
	if (!par_read_string(argc,argv,"uz", uz_name)) { docs (); exit (1); }
	if (!par_read_string(argc,argv,"mpp", mpp_name)) { docs (); exit (1); }
	if (!par_read_string(argc,argv,"vp", vp_name)) { docs (); exit (1); }
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
	if (!par_read_int(argc,argv,"Lw",&Lw)) Lw = 300;
	if (!par_read_float(argc,argv,"kz_eps",&kz_eps)) kz_eps = 0.0005;
	
	// get dimensions from velocity (nz,oz,dz,nx,ox,dx) and wavelet (nt) files
	InitFileHeader(&fh);
	ReadFileHeader(wav_name,&fh);
	nt = fh.n1;
	ReadFileHeader(vp_name,&fh);
	nz = fh.n1; ntraces = fh.n2*fh.n3*fh.n4*fh.n5;

	h_wav = allocSeisHeader(1);
	wav = alloc2float(nt,1); 		
	SeisRead(wav_name,wav,h_wav,&fh);
	ot = fh.o1;
	dt = fh.d1;
	h_vp = allocSeisHeader(ntraces);
	vp = alloc2float(nz,ntraces);
	SeisRead(vp_name,vp,h_vp,&fh); 	
	nx = h_vp[ntraces-1].imx - h_vp[0].imx + 1;
	oz = fh.o1;
	dz = fh.d1;
    	ox = h_vp[0].mx;
    	dx = h_vp[1].mx - h_vp[0].mx;

	if (adj){
		h_uz = allocSeisHeader(nx);
		uz = alloc2float(nt,nx);
		SeisRead(uz_name,uz,h_uz,&fh);
	}        
	else{
		h_mpp = allocSeisHeader(nx*nhx);
		mpp = alloc2float(nz,nx*nhx);
		SeisRead(mpp_name,mpp,h_mpp,&fh);
	}
	if (adj){
		mpp = alloc2float(nz,nx*nhx);
		for (iz=0;iz<nz;iz++) for (ix=0;ix<nx*nhx;ix++) mpp[ix][iz] = 0.0;
		h_mpp = allocSeisHeader(nx*nhx);
		for (ix=0;ix<nx;ix++){
			for (ihx=0;ihx<nhx;ihx++){ 
				h_mpp[ix*nhx + ihx] = h_vp[ix];
				h_mpp[ix*nhx + ihx].ihx = ihx;
				h_mpp[ix*nhx + ihx].hx = ihx*dhx + ohx;
			}
		}
	}
	else{
		h_uz = allocSeisHeader(nx);
		uz = alloc2float(nt,nx);
		for (it=0;it<nt;it++) for (ix=0;ix<nx;ix++) uz[ix][it] = 0.0;
		for (ix=0;ix<nx;ix++) h_uz[ix] = h_mpp[ix];
		for (ix=0;ix<nx;ix++) h_uz[ix].n1 = nt;
		for (ix=0;ix<nx;ix++) h_uz[ix].d1 = dt;
		for (ix=0;ix<nx;ix++) h_uz[ix].o1 = ot;
		for (ix=0;ix<nx;ix++) h_uz[ix].sx = sx;
		for (ix=0;ix<nx;ix++) h_uz[ix].gx = h_mpp[ix].mx;
		for (ix=0;ix<nx;ix++) h_uz[ix].hx = h_mpp[ix].mx - sx;
	}

	wem_suboff(uz,
	     mpp, 
	     wav,
	     nt,ot,dt, 
	     nx,ox,dx,
	     nhx,ohx,dhx,
	     sx,
	     nz,oz,dz,gz,sz,
	     vp,nref, 
	     fmin,fmax,
	     padt,padx,
	     adj,pspi,verbose,
	     kz_eps);
	if (adj){
		InitFileHeader(&fh);
		fh.n1 = nz; fh.o1 = oz; fh.d1 = dz;
		fh.n2 = nhx; fh.o2 = ohx; fh.d2 = dhx;
		fh.n3 = nx; fh.o3 = ox; fh.d3 = dx;
		fh.n4 = 1; fh.o4 = 0.; fh.d4 = 1.;
		fh.n5 = 1; fh.o5 = 0.; fh.d5 = 1.;
		for (ix=0;ix<nx*nhx;ix++) h_mpp[ix].trid = 1;
		for (ix=0;ix<nx*nhx;ix++) h_mpp[ix].mx = h_uz[ix].gx;
		for (ix=0;ix<nx*nhx;ix++) h_mpp[ix].imx = (int) (h_mpp[ix].mx - ox)/dx;
		SeisWrite(mpp_name,mpp,h_mpp,&fh);
	}
	else{
		InitFileHeader(&fh);
		fh.n1 = nt; fh.o1 = ot; fh.d1 = dt;
		fh.n2 = nx; fh.o2 = ox; fh.d2 = dx;
		fh.n3 = 1; fh.o3 = 0; fh.d3 = 1;
		SeisWrite(uz_name,uz,h_uz,&fh);
	}
	free2float(uz); 
	free2float(mpp); 
	free2float(vp); 
	free2float(wav); 
	freeSeisHeader(h_uz);
	freeSeisHeader(h_mpp);
	freeSeisHeader(h_vp);
	freeSeisHeader(h_wav);
	return 0;
}
