#include "seismic.h"
#include "ewem.h"
void ewem(float **ux, float **uy, float **uz,
	  float **mpp, float **mps1, float **mps2,
	  float **wav,
	  int nt, float ot, float dt, 
	  int nmx,float omx, float dmx,
	  int nmy,float omy, float dmy,
	  float sx,float sy,
	  int nz, float oz, float dz, float gz, float sz,
	  float **velp, float **vels, int nref, int Lw, 
	  float fmin, float fmax,
	  int padt, int padx,
	  bool adj, bool pspi, bool verbose,
	  float kz_eps)
/*< elastic wave equation depth migration operator. >*/
{
	int iz,ix,imx,imy,igx,igy,ik,iw,it,nw,nkx,nky,ntfft;
	float dw,dkx,dky;
	int ifmin,ifmax;
	float *d_t;
	complex *d_w;
	complex **ux_g_wx,**uy_g_wx,**uz_g_wx,**u_s_wx;
	fftwf_complex *a,*b;
	int *n;
	fftwf_plan p1,p2;
	float *po_p,**pd_p;
	float *po_s,**pd_s;
	float progress;
	int ithread,nthread;
	float max_source;
	float **mpp_threads,**mps1_threads,**mps2_threads;
	float **vpref,vpmin,vpmax,vp;
	int **ipref1,**ipref2;
	float **vsref,vsmin,vsmax,vs;
	int **isref1,**isref2;
	int iref;
	
	/* decompose slowness into layer average, and layer purturbation */
	po_p = alloc1float(nz); 
	pd_p = alloc2float(nz,nmx*nmy); 
	for (iz=0;iz<nz;iz++){
		po_p[iz] = 0.;
		for (ix=0;ix<nmx*nmy;ix++) po_p[iz] += velp[ix][iz];
		po_p[iz] /= (float) nmx*nmy;
		po_p[iz]  = 1./po_p[iz];
		for (ix=0;ix<nmx*nmy;ix++) pd_p[ix][iz] = 1.0/velp[ix][iz] - po_p[iz];
	}
	po_s = alloc1float(nz); 
	pd_s = alloc2float(nz,nmx*nmy); 
	for (iz=0;iz<nz;iz++){
		po_s[iz] = 0.;
		for (ix=0;ix<nmx*nmy;ix++) po_s[iz] += vels[ix][iz];
		po_s[iz] /= (float) nmx*nmy;
		po_s[iz]  = 1./po_s[iz];
		for (ix=0;ix<nmx*nmy;ix++) pd_s[ix][iz] = 1.0/vels[ix][iz] - po_s[iz];
	}

	// P wave reference velocities
	/****************************************************************************************/
	/* generate reference velocities for each depth step */
	vpref = alloc2float(nz,nref); /* reference velocities for each layer */
	ipref1 = alloc2int(nz,nmx*nmy); /* index of nearest lower reference velocity for each subsurface point */
	ipref2 = alloc2int(nz,nmx*nmy); /* index of nearest upper reference velocity for each subsurface point */
	for (iz=0;iz<nz;iz++){
		vpmin=velp[0][iz];
		for (ix=0;ix<nmx;ix++) if (velp[ix][iz] < vpmin) vpmin = velp[ix][iz];
		vpmax=velp[nmx-1][iz];
		for (ix=0;ix<nmx*nmy;ix++) if (velp[ix][iz] > vpmax) vpmax = velp[ix][iz];

		for (iref=0;iref<nref;iref++) vpref[iref][iz] = vpmin + (float) iref*(vpmax-vpmin)/((float) nref-1);
		for (ix=0;ix<nmx*nmy;ix++){
			vp = velp[ix][iz];
			if (vpmax>vpmin+10){
				iref = (int) (nref-1)*(vp-vpmin)/(vpmax-vpmin);
				ipref1[ix][iz] = iref;
				ipref2[ix][iz] = iref+1;
				if (iref>nref-2){
					ipref1[ix][iz] = nref-1;
					ipref2[ix][iz] = nref-1;
				}
			}
			else{
				ipref1[ix][iz] = 0;
				ipref2[ix][iz] = 0;
			}
		}
	}
	/****************************************************************************************/

	// S wave reference velocities
	/****************************************************************************************/
	/* generate reference velocities for each depth step */
	vsref = alloc2float(nz,nref); /* reference velocities for each layer */
	isref1 = alloc2int(nz,nmx*nmy); /* index of nearest lower reference velocity for each subsurface point */
	isref2 = alloc2int(nz,nmx*nmy); /* index of nearest upper reference velocity for each subsurface point */
	for (iz=0;iz<nz;iz++){
		vsmin=vels[0][iz];
		for (ix=0;ix<nmx;ix++) if (vels[ix][iz] < vsmin) vsmin = vels[ix][iz];
		vsmax=vels[nmx-1][iz];
		for (ix=0;ix<nmx*nmy;ix++) if (vels[ix][iz] > vsmax) vsmax = vels[ix][iz];

		for (iref=0;iref<nref;iref++) vsref[iref][iz] = vsmin + (float) iref*(vsmax-vsmin)/((float) nref-1);
		for (ix=0;ix<nmx*nmy;ix++){
			vs = vels[ix][iz];
			if (vsmax>vsmin+10){
				iref = (int) (nref-1)*(vs-vsmin)/(vsmax-vsmin);
				isref1[ix][iz] = iref;
				isref2[ix][iz] = iref+1;
				if (iref>nref-2){
					isref1[ix][iz] = nref-1;
					isref2[ix][iz] = nref-1;
				}
			}
			else{
				isref1[ix][iz] = 0;
				isref2[ix][iz] = 0;
			}
		}
	}
	/****************************************************************************************/
		
	if (adj){
		for (ix=0;ix<nmx*nmy;ix++) for (iz=0;iz<nz;iz++) mpp[ix][iz] = 0.;
		for (ix=0;ix<nmx*nmy;ix++) for (iz=0;iz<nz;iz++) mps1[ix][iz] = 0.;
		for (ix=0;ix<nmx*nmy;ix++) for (iz=0;iz<nz;iz++) mps2[ix][iz] = 0.;
	}
	else{
		for (ix=0;ix<nmx*nmy;ix++) for (it=0;it<nt;it++) ux[ix][it] = 0.;
		for (ix=0;ix<nmx*nmy;ix++) for (it=0;it<nt;it++) uy[ix][it] = 0.;
		for (ix=0;ix<nmx*nmy;ix++) for (it=0;it<nt;it++) uz[ix][it] = 0.;
	}
	ntfft = (int) 2*truncf(padt*((float) nt)/2);
	nw = (int) truncf(ntfft/2)+1;
	nkx = nmx > 1 ? padx*nmx : 1;
	nky = nmy > 1 ? padx*nmy : 1;
	dkx = 2*PI/((float) nkx)/dmx;
	dky = 2*PI/((float) nky)/dmy;
	dw = 2*PI/((float) ntfft)/dt;
	if(fmax*dt*ntfft+1<nw) ifmax = trunc(fmax*dt*ntfft)+1;
	else ifmax = nw;
	if(fmin*dt*ntfft+1<ifmax) ifmin = trunc(fmin*dt*ntfft);
	else ifmin = 0;
	ux_g_wx = alloc2complex(nw,nmx*nmy);
	uy_g_wx = alloc2complex(nw,nmx*nmy);
	uz_g_wx = alloc2complex(nw,nmx*nmy);
	u_s_wx = alloc2complex(nw,nmx*nmy);
	d_t = alloc1float(nt);
	d_w = alloc1complex(nw);
	for (it=0;it<nt;it++)  d_t[it] = 0.;  
	for (iw=0;iw<nw;iw++)  d_w[iw] = 0.; 

	/* set up fftw plans and pass them to the OMP region of the code */
	a  = fftwf_malloc(sizeof(fftwf_complex) * nkx*nky);
	b  = fftwf_malloc(sizeof(fftwf_complex) * nkx*nky);
	n = alloc1int(2); 
	n[0] = nkx;
	n[1] = nky;
	p1 = fftwf_plan_dft(2, n, a, a, FFTW_FORWARD, FFTW_MEASURE);
	p2 = fftwf_plan_dft(2, n, b, b, FFTW_BACKWARD, FFTW_MEASURE);
	for (ik=0;ik<nkx*nky;ik++){
		a[ik] = 0.;
		b[ik] = 0.;
	} 
	fftwf_execute_dft(p1,a,a);
	fftwf_execute_dft(p2,b,b);

	/**********************************************************************/
	igx = (int) (sx - omx)/dmx; /*position to inject source in x-dir*/
	igy = (int) (sy - omy)/dmy; /*position to inject source in y-dir*/
	/* source wavefield*/
	for (ix=0;ix<nmx*nmy;ix++) for (iw=0;iw<nw;iw++) u_s_wx[ix][iw] = 0.;
	for (it=0;it<nt;it++) d_t[it] = wav[0][it];
	f_op(d_w,d_t,nw,nt,1); /* d_t to d_w */
	for (iw=0;iw<nw;iw++) u_s_wx[igx*nmy + igy][iw] = d_w[iw];
	/* receiver wavefield*/
	if (adj){
		for (ix=0;ix<nmx*nmy;ix++){
			// x component
			for (it=0;it<nt;it++) d_t[it] = ux[ix][it];
			f_op(d_w,d_t,nw,nt,1); /* d_t to d_w */
			for (iw=0;iw<ifmin;iw++) ux_g_wx[ix][iw] = 0.;
			for (iw=ifmin;iw<ifmax;iw++) ux_g_wx[ix][iw] = d_w[iw];
			for (iw=ifmax;iw<nw;iw++) ux_g_wx[ix][iw] = 0.;
			// y component
			for (it=0;it<nt;it++) d_t[it] = uy[ix][it];
			f_op(d_w,d_t,nw,nt,1); /* d_t to d_w */
			for (iw=0;iw<ifmin;iw++) uy_g_wx[ix][iw] = 0.;
			for (iw=ifmin;iw<ifmax;iw++) uy_g_wx[ix][iw] = d_w[iw];
			for (iw=ifmax;iw<nw;iw++) uy_g_wx[ix][iw] = 0.;
			// z component
			for (it=0;it<nt;it++) d_t[it] = uz[ix][it];
			f_op(d_w,d_t,nw,nt,1); /* d_t to d_w */
			for (iw=0;iw<ifmin;iw++) uz_g_wx[ix][iw] = 0.;
			for (iw=ifmin;iw<ifmax;iw++) uz_g_wx[ix][iw] = d_w[iw];
			for (iw=ifmax;iw<nw;iw++) uz_g_wx[ix][iw] = 0.;
		}
	}
	else{
		for (ix=0;ix<nmx*nmy;ix++){
			for (iw=0;iw<nw;iw++){
				ux_g_wx[ix][iw] = 0.;
				uy_g_wx[ix][iw] = 0.;
				uz_g_wx[ix][iw] = 0.;
			}
		}
	}
	max_source = 0.;
	for (it=0;it<nt;it++) if (max_source < fabsf(wav[0][it])/sqrtf((float) ntfft)) max_source = fabsf(wav[0][it])/sqrtf((float) ntfft);

	nthread = omp_thread_count();
	//fprintf(stderr,"nthread=%d\n",nthread);
	if (adj){
		mpp_threads = alloc2float(nz,nmx*nmy*nthread);
		mps1_threads = alloc2float(nz,nmx*nmy*nthread);
		mps2_threads = alloc2float(nz,nmx*nmy*nthread);
		for (imx=0;imx<nmx;imx++){
			for (imy=0;imy<nmy;imy++){
				for (ithread=0;ithread<nthread;ithread++){
					for (iz=0;iz<nz;iz++){
						mpp_threads[imx*nmy*nthread + imy*nthread + ithread][iz] = 0.;
						mps1_threads[imx*nmy*nthread + imy*nthread + ithread][iz] = 0.;
						mps2_threads[imx*nmy*nthread + imy*nthread + ithread][iz] = 0.;
					}
				}
			}
		}
	}
	else{
		mpp_threads = alloc2float(nz,nmx*nmy);
		mps1_threads = alloc2float(nz,nmx*nmy);
		mps2_threads = alloc2float(nz,nmx*nmy);
		for (imx=0;imx<nmx;imx++){
			for (imy=0;imy<nmy;imy++){
				for (iz=0;iz<nz;iz++){
					mpp_threads[imx*nmy + imy][iz] = mpp[imx*nmy + imy][iz];
					mps1_threads[imx*nmy + imy][iz] = mps1[imx*nmy + imy][iz];
					mps2_threads[imx*nmy + imy][iz] = mps2[imx*nmy + imy][iz];
				}
			}
		}
	}

	progress = 0.;
#pragma omp parallel for private(iw) shared(mpp_threads,mps1_threads,mps2_threads,ux_g_wx,uy_g_wx,uz_g_wx,u_s_wx)
	for (iw=ifmin;iw<ifmax;iw++){ 
		progress += 1./((float) ifmax - ifmin);
		if (verbose) progress_msg(progress);
		elastic_extrap1f(mpp_threads,mps1_threads,mps2_threads,
				ux_g_wx,uy_g_wx,uz_g_wx,
				u_s_wx,
				max_source,iw,nw,ifmax,ntfft,
				dw,dkx,dky,nkx,nky,
				nz,oz,dz,gz,sz,nmx,omx,dmx,nmy,omy,dmy,
				nthread,
				velp,po_p,pd_p,vels,po_s,pd_s,
				vpref,ipref1,ipref2,nref,
				vsref,isref1,isref2,
				Lw,
				p1,p2,adj,pspi,verbose,kz_eps);
	}
	if (verbose) fprintf(stderr,"\n");

	if (adj){
		// reduction over parallel axis 
		for (imx=0;imx<nmx;imx++){ 
			for (imy=0;imy<nmy;imy++){
				for (ithread=0;ithread<nthread;ithread++){
					for (iz=0;iz<nz;iz++){
						mpp[imx*nmy + imy][iz] += mpp_threads[imx*nmy*nthread + imy*nthread + ithread][iz];
						mps1[imx*nmy + imy][iz] += mps1_threads[imx*nmy*nthread + imy*nthread + ithread][iz];
						mps2[imx*nmy + imy][iz] += mps2_threads[imx*nmy*nthread + imy*nthread + ithread][iz];
					}
				}
			}
		}
	}
	else{
		for (ix=0;ix<nmx*nmy;ix++){
			// x component
			for (iw=0;iw<ifmin;iw++) d_w[iw] = 0.;
			for (iw=ifmin;iw<ifmax;iw++) d_w[iw] = ux_g_wx[ix][iw];
			for (iw=ifmax;iw<nw;iw++) d_w[iw] = 0.;
			f_op(d_w,d_t,nw,nt,0); /* d_w to d_t */
			for (it=0;it<nt;it++) ux[ix][it] = d_t[it];
			// y component
			for (iw=0;iw<ifmin;iw++) d_w[iw] = 0.;
			for (iw=ifmin;iw<ifmax;iw++) d_w[iw] = uy_g_wx[ix][iw];
			for (iw=ifmax;iw<nw;iw++) d_w[iw] = 0.;
			f_op(d_w,d_t,nw,nt,0); /* d_w to d_t */
			for (it=0;it<nt;it++) uy[ix][it] = d_t[it];
			// x component
			for (iw=0;iw<ifmin;iw++) d_w[iw] = 0.;
			for (iw=ifmin;iw<ifmax;iw++) d_w[iw] = uz_g_wx[ix][iw];
			for (iw=ifmax;iw<nw;iw++) d_w[iw] = 0.;
			f_op(d_w,d_t,nw,nt,0); /* d_w to d_t */
			for (it=0;it<nt;it++) uz[ix][it] = d_t[it];
		}
	}

	free1int(n); 
	fftwf_free(a);
	fftwf_free(b);
	fftwf_destroy_plan(p1);
	fftwf_destroy_plan(p2);
	free1float(d_t);
	free1complex(d_w);
	free2complex(ux_g_wx);
	free2complex(uy_g_wx);
	free2complex(uz_g_wx);
	free2complex(u_s_wx);
	free1float(po_p);
	free2float(pd_p);
	free1float(po_s);
	free2float(pd_s);
	free2float(vpref);
	free2int(ipref1);
	free2int(ipref2);
	free2float(mpp_threads);
	free2float(mps1_threads);
	free2float(mps2_threads);

	return;
} 

void elastic_extrap1f(float **mpp, float **mps1, float **mps2,
		complex **ux_g_wx, complex **uy_g_wx, complex **uz_g_wx, 
		complex **u_s_wx,
		float max_source, int iw, int nw,int ifmax,int ntfft,float dw,float dkx,float dky,int nkx,int nky,
		int nz, float oz, float dz, float gz, float sz,
		int nmx,float omx, float dmx,
		int nmy,float omy, float dmy,
		int nthread,
		float **vp,float *po_p,float **pd_p,
		float **vs,float *po_s,float **pd_s,
		float **vpref, int **ipref1, int **ipref2, int nref,
		float **vsref, int **isref1, int **isref2,
		int Lw,
		fftwf_plan p1,fftwf_plan p2,
		bool adj, bool pspi, bool verbose,
		float kz_eps)
/*< extrapolate 1 frequency >*/
{
	float w,factor,z;
	int iz,ix,imx,imy,ithread;
	complex *ux_xg,*uy_xg,*uz_xg;
	complex *up_xg,*us1_xg,*us2_xg;
	complex *up_xs;
	complex **smig;

	ithread = omp_get_thread_num(); 
	//fprintf(stderr,"ithread=%d\n",ithread);

	ux_xg = alloc1complex(nmx*nmy);
	uy_xg = alloc1complex(nmx*nmy);
	uz_xg = alloc1complex(nmx*nmy);
	up_xg = alloc1complex(nmx*nmy);
	us1_xg = alloc1complex(nmx*nmy);
	us2_xg = alloc1complex(nmx*nmy);
	up_xs = alloc1complex(nmx*nmy);
	smig = alloc2complex(nz,nmx*nmy);

	for (ix=0;ix<nmx*nmy;ix++) ux_xg[ix] = 0.;
	for (ix=0;ix<nmx*nmy;ix++) uy_xg[ix] = 0.;
	for (ix=0;ix<nmx*nmy;ix++) uz_xg[ix] = 0.;
	for (ix=0;ix<nmx*nmy;ix++) up_xg[ix] = 0.;
	for (ix=0;ix<nmx*nmy;ix++) us1_xg[ix] = 0.;
	for (ix=0;ix<nmx*nmy;ix++) us2_xg[ix] = 0.;
	for (ix=0;ix<nmx*nmy;ix++) up_xs[ix] = 0.;
	if (iw==0) factor = 1.;
	else factor = 2.;

	w = iw*dw;
	for (ix=0;ix<nmx*nmy;ix++){ 
		up_xs[ix] = u_s_wx[ix][iw]/sqrtf((float) ntfft);
		ux_xg[ix] = ux_g_wx[ix][iw]/sqrtf((float) ntfft); // *pow(w,2)
		uy_xg[ix] = uy_g_wx[ix][iw]/sqrtf((float) ntfft); // *pow(w,2)
		uz_xg[ix] = uz_g_wx[ix][iw]/sqrtf((float) ntfft); // *pow(w,2)

	}
	for (ix=0;ix<nmx*nmy;ix++) up_xs[ix] = u_s_wx[ix][iw]/sqrtf((float) ntfft);
	for (iz=0;iz<nz;iz++){ // extrapolate source wavefield 
		z = oz + dz*iz;
		if (z >= sz){
			if (pspi) pspiop(up_xs,w,dkx,dky,nkx,nky,nmx,omx,dmx,nmy,omy,dmy,dz,iz,vp,po_p,pd_p,vpref,ipref1,ipref2,nref,p1,p2,true,true,verbose,kz_eps);
			else 	  ssop(up_xs,w,dkx,dky,nkx,nky,nmx,omx,dmx,nmy,omy,dmy,dz,iz,vp,po_p,pd_p,p1,p2,true,true,verbose,kz_eps);
			for (ix=0;ix<nmx*nmy;ix++) smig[ix][iz]  = up_xs[ix]/max_source;
		}
		else{
			for (ix=0;ix<nmx*nmy;ix++) smig[ix][iz] = 0.;
		}
	}
	if (adj){
		for (iz=0;iz<nz;iz++){ // extrapolate receiver wavefield
			z = oz + dz*iz;
			if (z >= gz){
				elastic_separate_3d(ux_xg,uy_xg,uz_xg,up_xg,us1_xg,us2_xg,w,dkx,nkx,nmx,omx,dmx,dky,nky,nmy,omy,dmy,1./po_p[iz],1./po_s[iz],p1,p2,true,adj);
				//elastic_separate_3d_window(ux_xg,uy_xg,uz_xg,up_xg,us1_xg,us2_xg,w,dkx,nkx,nmx,omx,dmx,dky,nky,nmy,omy,dmy,vp,vs,iz,Lw,p1,p2,true,adj);
				if (pspi){ 
					pspiop(up_xg,w,dkx,dky,nkx,nky,nmx,omx,dmx,nmy,omy,dmy,dz,iz,vp,po_p,pd_p,vpref,ipref1,ipref2,nref,p1,p2,true,false,verbose,kz_eps);
					pspiop(us1_xg,w,dkx,dky,nkx,nky,nmx,omx,dmx,nmy,omy,dmy,dz,iz,vs,po_s,pd_s,vsref,isref1,isref2,nref,p1,p2,true,false,verbose,kz_eps);
					pspiop(us2_xg,w,dkx,dky,nkx,nky,nmx,omx,dmx,nmy,omy,dmy,dz,iz,vs,po_s,pd_s,vsref,isref1,isref2,nref,p1,p2,true,false,verbose,kz_eps);
				}
				else{
					ssop(up_xg,w,dkx,dky,nkx,nky,nmx,omx,dmx,nmy,omy,dmy,dz,iz,vp,po_p,pd_p,p1,p2,true,false,verbose,kz_eps);
					ssop(us1_xg,w,dkx,dky,nkx,nky,nmx,omx,dmx,nmy,omy,dmy,dz,iz,vs,po_s,pd_s,p1,p2,true,false,verbose,kz_eps);
					ssop(us2_xg,w,dkx,dky,nkx,nky,nmx,omx,dmx,nmy,omy,dmy,dz,iz,vs,po_s,pd_s,p1,p2,true,false,verbose,kz_eps);
				}
				elastic_separate_3d(ux_xg,uy_xg,uz_xg,up_xg,us1_xg,us2_xg,w,dkx,nkx,nmx,omx,dmx,dky,nky,nmy,omy,dmy,1./po_p[iz],1./po_s[iz],p1,p2,false,adj);
				//elastic_separate_3d_window(ux_xg,uy_xg,uz_xg,up_xg,us1_xg,us2_xg,w,dkx,nkx,nmx,omx,dmx,dky,nky,nmy,omy,dmy,vp,vs,iz,Lw,p1,p2,false,adj);
				for (imx=0;imx<nmx;imx++){ 
					for (imy=0;imy<nmy;imy++){
						mpp[imx*nmy*nthread + imy*nthread + ithread][iz]  += factor*crealf(conjf(smig[imx*nmy + imy][iz])*up_xg[imx*nmy + imy]);
						mps1[imx*nmy*nthread + imy*nthread + ithread][iz] += factor*crealf(conjf(smig[imx*nmy + imy][iz])*us1_xg[imx*nmy + imy]);
						mps2[imx*nmy*nthread + imy*nthread + ithread][iz] += factor*crealf(conjf(smig[imx*nmy + imy][iz])*us2_xg[imx*nmy + imy]);
					}
				}
			}
		}
	}
	else{
		for (iz=nz-1;iz>=0;iz--){ // extrapolate receiver wavefield 
			z = oz + dz*iz;
			if (z >= gz){
				elastic_separate_3d(ux_xg,uy_xg,uz_xg,up_xg,us1_xg,us2_xg,w,dkx,nkx,nmx,omx,dmx,dky,nky,nmy,omy,dmy,1./po_p[iz],1./po_s[iz],p1,p2,true,adj);
				//elastic_separate_3d_window(ux_xg,uy_xg,uz_xg,up_xg,us1_xg,us2_xg,w,dkx,nkx,nmx,omx,dmx,dky,nky,nmy,omy,dmy,vp,vs,iz,Lw,p1,p2,true,adj);
				for (ix=0;ix<nmx*nmy;ix++) up_xg[ix]  += smig[ix][iz]*mpp[ix][iz];
				for (ix=0;ix<nmx*nmy;ix++) us1_xg[ix] += smig[ix][iz]*mps1[ix][iz];
				for (ix=0;ix<nmx*nmy;ix++) us2_xg[ix] += smig[ix][iz]*mps2[ix][iz];
				if (pspi){ 
					pspiop(up_xg,w,dkx,dky,nkx,nky,nmx,omx,dmx,nmy,omy,dmy,dz,iz,vp,po_p,pd_p,vpref,ipref1,ipref2,nref,p1,p2,false,false,verbose,kz_eps);
					pspiop(us1_xg,w,dkx,dky,nkx,nky,nmx,omx,dmx,nmy,omy,dmy,dz,iz,vs,po_s,pd_s,vsref,isref1,isref2,nref,p1,p2,false,false,verbose,kz_eps);
					pspiop(us2_xg,w,dkx,dky,nkx,nky,nmx,omx,dmx,nmy,omy,dmy,dz,iz,vs,po_s,pd_s,vsref,isref1,isref2,nref,p1,p2,false,false,verbose,kz_eps);
				}
				else{
					ssop(up_xg,w,dkx,dky,nkx,nky,nmx,omx,dmx,nmy,omy,dmy,dz,iz,vp,po_p,pd_p,p1,p2,false,false,verbose,kz_eps);
					ssop(us1_xg,w,dkx,dky,nkx,nky,nmx,omx,dmx,nmy,omy,dmy,dz,iz,vs,po_s,pd_s,p1,p2,false,false,verbose,kz_eps);
					ssop(us2_xg,w,dkx,dky,nkx,nky,nmx,omx,dmx,nmy,omy,dmy,dz,iz,vs,po_s,pd_s,p1,p2,false,false,verbose,kz_eps);
				}
				elastic_separate_3d(ux_xg,uy_xg,uz_xg,up_xg,us1_xg,us2_xg,w,dkx,nkx,nmx,omx,dmx,dky,nky,nmy,omy,dmy,1./po_p[iz],1./po_s[iz],p1,p2,false,adj);
				//elastic_separate_3d_window(ux_xg,uy_xg,uz_xg,up_xg,us1_xg,us2_xg,w,dkx,nkx,nmx,omx,dmx,dky,nky,nmy,omy,dmy,vp,vs,iz,Lw,p1,p2,false,adj);
			}
		}
		for (ix=0;ix<nmx*nmy;ix++){
			ux_g_wx[ix][iw] = ux_xg[ix]/sqrtf((float) ntfft); // *pow(w,2)
			uy_g_wx[ix][iw] = uy_xg[ix]/sqrtf((float) ntfft); // *pow(w,2)
			uz_g_wx[ix][iw] = uz_xg[ix]/sqrtf((float) ntfft); // *pow(w,2)
		}
	}

	free1complex(ux_xg);
	free1complex(uy_xg);
	free1complex(uz_xg);
	free1complex(up_xg);
	free1complex(us1_xg);
	free1complex(us2_xg);
	free1complex(up_xs);
	free2complex(smig);

	return;
}

void ssop(complex *d_x,
		float w,float dkx,float dky,int nkx,int nky,int nmx,float omx,float dmx,int nmy,float omy,float dmy,float dz,int iz,
		float **v,float *po,float **pd,
		fftwf_plan p1,fftwf_plan p2,
		bool adj, 
		bool src,
		bool verbose,
		float kz_eps)
{
	float kx,ky;
	complex L;
	int ik,ikx,iky,imx,imy; 
	complex *d_k;
	fftwf_complex *a,*b;
	complex w2,s;	

	int lmx,lmy;
	if (nmx>100) lmx=30;
	else lmx=0;
	if (nmy>100) lmy=30;
	else lmy=0;
	a  = fftwf_malloc(sizeof(fftwf_complex) * nkx*nky);
	b  = fftwf_malloc(sizeof(fftwf_complex) * nkx*nky);
	d_k = alloc1complex(nkx*nky);
	w2 = kz_eps + w*I;
	w2 *= w2;
	w2 *= po[iz];
	w2 *= po[iz];
	if (adj){
		for(imx=0; imx<nkx;imx++){ 
			for(imy=0; imy<nky;imy++){ 
				if (imx < nmx && imy < nmy) a[imx*nky + imy] = d_x[imx*nmy + imy];
				else a[imx*nky + imy] = 0.;
			}
		}
		fftwf_execute_dft(p1,a,a); 
		for (ikx=0;ikx<nkx;ikx++){
			if (ikx<= (int) nkx/2) kx = (float) dkx*ikx;
			else                   kx = -((float) dkx*nkx - dkx*ikx);
			for (iky=0;iky<nky;iky++){
				if (iky<= (int) nky/2) ky = (float) dky*iky;
				else                   ky = -((float) dky*nky - dky*iky);
				s = w2 + powf(kx,2) + powf(ky,2);
				L = cexp(-csqrtf(s)*dz);
				if (src) d_k[ikx*nky + iky] = ((complex) a[ikx*nky + iky])*L/sqrtf((float) nkx*nky);
				else     d_k[ikx*nky + iky] = ((complex) a[ikx*nky + iky])*conjf(L)/sqrtf((float) nkx*nky);   
			}
		}
		for(ik=0; ik<nkx*nky;ik++) b[ik] = (fftwf_complex) d_k[ik];
		fftwf_execute_dft(p2,b,b);
		for(imx=0; imx<nkx;imx++){ 
			for(imy=0; imy<nky;imy++){ 
				if (imx < nmx && imy < nmy){
					L = cexpf(I*w*pd[imx*nmy + imy][iz]*dz);
					if (src) d_x[imx*nmy + imy] = ((complex) b[imx*nky + imy])*conjf(L)/sqrtf((float) nkx*nky); // SS operator
					else     d_x[imx*nmy + imy] = ((complex) b[imx*nky + imy])*L/sqrtf((float) nkx*nky); // SS operator
				}
			}
		}
		boundary_condition(d_x,nmx,lmx,nmy,lmy);    
	}
	else{
		boundary_condition(d_x,nmx,lmx,nmy,lmy);    
		for(imx=0; imx<nkx;imx++){ 
			for(imy=0; imy<nky;imy++){ 
				if (imx < nmx && imy < nmy){
					L = cexpf(I*w*pd[imx*nmy + imy][iz]*dz);  
					a[imx*nky + imy] = d_x[imx*nmy + imy]*conjf(L); // SS operator
				}
				else a[imx*nky + imy] = 0.0;
			}
		}
		fftwf_execute_dft(p1,a,a); 
		for (ikx=0;ikx<nkx;ikx++){
			if (ikx<= (int) nkx/2) kx = (float) dkx*ikx;
			else                   kx = -((float) dkx*nkx - dkx*ikx);
			for (iky=0;iky<nky;iky++){
				if (iky<= (int) nky/2) ky = (float) dky*iky;
				else                   ky = -((float) dky*nky - dky*iky);
				s = w2 + powf(kx,2) + powf(ky,2);
				L = cexp(-csqrtf(s)*dz);
				d_k[ikx*nky + iky] = ((complex) a[ikx*nky + iky])*L/sqrtf((float) nkx*nky);
			}
		}
		for(ik=0; ik<nkx*nky;ik++) b[ik] = (fftwf_complex) d_k[ik];
		fftwf_execute_dft(p2,b,b);
		for(imx=0; imx<nkx;imx++){ 
			for(imy=0; imy<nky;imy++){ 
				if (imx < nmx && imy < nmy){ 
					d_x[imx*nmy + imy] = ((complex) b[imx*nky + imy])/sqrtf((float) nkx*nky);
				}
			}
		}
	}
	free1complex(d_k);
	fftwf_free(a);
	fftwf_free(b);

	return;
}

void pspiop(complex *d_x,
		float w,float dkx,float dky,int nkx,int nky,
		int nmx,float omx,float dmx,int nmy,float omy,float dmy,
		float dz,int iz,
		float **vel,float *po,float **pd,
		float **vref, int **iref1, int **iref2, int nref,
		fftwf_plan p1,fftwf_plan p2,
		bool adj, 
		bool src,
		bool verbose,
		float kz_eps)
{


	float kx,ky;
	complex L;
	int ik,ikx,iky,imx,imy,ix; 
	complex *d_k;
	fftwf_complex *a,*b;
	int lmx,lmy;
	complex **dref;
	complex w2,s;	
	int iref;
	float v,vref1,vref2,aa,pd_ref;
	if (nmx>100) lmx=30;
	else lmx=0;
	if (nmy>100) lmy=30;
	else lmy=0;
	a  = fftwf_malloc(sizeof(fftwf_complex) * nkx*nky);
	b  = fftwf_malloc(sizeof(fftwf_complex) * nkx*nky);
	dref = alloc2complex(nref,nmx*nmy);
	d_k = alloc1complex(nkx*nky);

	if (adj){
		for(imx=0; imx<nkx;imx++){ 
			for(imy=0; imy<nky;imy++){ 
				if (imx < nmx && imy < nmy) a[imx*nky + imy] = d_x[imx*nmy + imy];
				else a[imx*nky + imy] = 0.;
			}
		}
		fftwf_execute_dft(p1,a,a);
		for (iref=0;iref<nref;iref++){
			w2 = kz_eps + w*I;
			w2 *= w2;
			w2 /= vref[iref][iz];
			w2 /= vref[iref][iz];
			for (ikx=0;ikx<nkx;ikx++){
				if (ikx<= (int) nkx/2) kx = (float) dkx*ikx;
				else                   kx = -((float) dkx*nkx - dkx*ikx);
				for (iky=0;iky<nky;iky++){
					if (iky<= (int) nky/2) ky = (float) dky*iky;
					else                   ky = -((float) dky*nky - dky*iky);
					s = w2 + powf(kx,2) + powf(ky,2);
					L = cexp(-csqrtf(s)*dz);
					if (src) d_k[ikx*nky + iky] = ((complex) a[ikx*nky + iky])*L/sqrtf((float) nkx*nky);
					else     d_k[ikx*nky + iky] = ((complex) a[ikx*nky + iky])*conjf(L)/sqrtf((float) nkx*nky);
				}
			}
			for(ik=0; ik<nkx*nky;ik++) b[ik] = (fftwf_complex) d_k[ik];
			fftwf_execute_dft(p2,b,b);
			for(imx=0; imx<nkx;imx++){ 
				for(imy=0; imy<nky;imy++){ 
					if (imx < nmx && imy < nmy){
						pd_ref = (1.0/vel[imx*nmy + imy][iz]) - (1.0/vref[iref][iz]);
						L = cexpf(I*w*pd_ref*dz);
						if (src) dref[imx*nmy + imy][iref] = ((complex) b[imx*nky + imy])*conjf(L)/sqrtf((float) nkx*nky);
						else     dref[imx*nmy + imy][iref] = ((complex) b[imx*nky + imy])*L/sqrtf((float) nkx*nky);
					}
				}
			}
		}
		for (ix=0;ix<nmx*nmy;ix++){
			v = vel[ix][iz];
			vref1 = vref[iref1[ix][iz]][iz];
			vref2 = vref[iref2[ix][iz]][iz];
			aa = linear_interp(vref1,vref2,v);
			d_x[ix] = aa*dref[ix][iref1[ix][iz]] + (1.0-aa)*dref[ix][iref2[ix][iz]];
		}
		boundary_condition(d_x,nmx,lmx,nmy,lmy);
	}
	else{
		boundary_condition(d_x,nmx,lmx,nmy,lmy);    
		for (ix=0;ix<nmx*nmy;ix++) for (iref=0;iref<nref;iref++) dref[ix][iref] = 0.0;
		for (ix=0;ix<nmx*nmy;ix++){
			v = vel[ix][iz];
			vref1 = vref[iref1[ix][iz]][iz];
			vref2 = vref[iref2[ix][iz]][iz];
			aa = linear_interp(vref1,vref2,v);
			dref[ix][iref1[ix][iz]] = aa*d_x[ix];
			dref[ix][iref2[ix][iz]] = (1.0-aa)*d_x[ix];
		}
		for (ix=0;ix<nmx*nmy;ix++) d_x[ix] = 0.;
		for (iref=0;iref<nref;iref++){
			w2 = kz_eps + w*I;
			w2 *= w2;
			w2 /= vref[iref][iz];
			w2 /= vref[iref][iz];
			for(imx=0; imx<nkx;imx++){
				for(imy=0; imy<nky;imy++){
					if (imx < nmx && imy < nmy){ 
						pd_ref = (1.0/vel[imx*nmy + imy][iz]) - (1.0/vref[iref][iz]);
						L = cexpf(I*w*pd_ref*dz);  
						a[imx*nky + imy] = dref[imx*nmy + imy][iref]*conjf(L);
					}
					else a[imx*nky + imy] = 0.;
				}
			}
			fftwf_execute_dft(p1,a,a);
			for (ikx=0;ikx<nkx;ikx++){
				if (ikx<= (int) nkx/2) kx = (float) dkx*ikx;
				else                   kx = -((float) dkx*nkx - dkx*ikx);
				for (iky=0;iky<nky;iky++){
					if (iky<= (int) nky/2) ky = (float) dky*iky;
					else                   ky = -((float) dky*nky - dky*iky);
					s = w2 + powf(kx,2) + powf(ky,2);
					L = cexp(-csqrtf(s)*dz);
					d_k[ikx*nky + iky] = ((complex) a[ikx*nky + iky])*L/sqrtf((float) nkx*nky);
				}
			}
			for(ik=0; ik<nkx*nky;ik++) b[ik] = (fftwf_complex) d_k[ik];
			fftwf_execute_dft(p2,b,b);
			for(imx=0; imx<nkx;imx++){ 
				for(imy=0; imy<nky;imy++){ 
					if (imx < nmx && imy < nmy){
						d_x[imx*nmy + imy] += ((complex) b[imx*nky + imy])/sqrtf((float) nkx*nky);
					}
				}
			}
		}
	}
	free1complex(d_k);
	fftwf_free(a);
	fftwf_free(b);
	free2complex(dref);

	return;
}

float linear_interp(float x1,float x2,float x)
	/*< adjoint of linear interpolation between two points. aa is a scalar between 1 and 0 that is applied as aa*y1 and (1-aa)*y2 >*/
{	
	float aa;
	aa = x2 - x1 > 1.0 ? (x2-x)/(x2-x1) : 0.0;
	return  aa;
}

void f_op(complex *m,float *d,int nw,int nt,bool adj)
{
	fftwf_complex *out1a,*in1b;
	float *in1a,*out1b;
	int ntfft,it,iw;
	fftwf_plan p1a,p1b;

	ntfft = (nw-1)*2;

	if (adj){ /* data --> model */
		out1a = fftwf_malloc(sizeof(fftwf_complex) * nw);
		in1a = alloc1float(ntfft);
		p1a = fftwf_plan_dft_r2c_1d(ntfft, in1a, (fftwf_complex*)out1a, FFTW_ESTIMATE);
		for(it=0;it<nt;it++) in1a[it] = d[it];
		for(it=nt;it<ntfft;it++) in1a[it] = 0.;
		fftwf_execute(p1a); 
		for(iw=0;iw<nw;iw++) m[iw] = out1a[iw];
		fftwf_destroy_plan(p1a);
		fftwf_free(in1a); fftwf_free(out1a);
	}

	else{ /* model --> data */
		out1b = alloc1float(ntfft);
		in1b = fftwf_malloc(sizeof(fftwf_complex) * ntfft);
		p1b = fftwf_plan_dft_c2r_1d(ntfft, (fftwf_complex*)in1b, out1b, FFTW_ESTIMATE);
		for(iw=0;iw<nw;iw++) in1b[iw] = m[iw];
		for(iw=nw;iw<ntfft;iw++) in1b[iw] = 0.;
		fftwf_execute(p1b); 
		for(it=0;it<nt;it++) d[it] = out1b[it];
		fftwf_destroy_plan(p1b);
		fftwf_free(in1b); fftwf_free(out1b);
	}

	return;
}

void progress_msg(float progress)
{ 
	fprintf(stderr,"\r[%6.2f%% complete]      ",progress*100);
	return;
}

float signf(float a)
	/*< sign of a float >*/
{
	float b;
	if (a>0.)      b = 1.;
	else if (a<0.) b =-1.;
	else          b = 0.;
	return b;
}

float signfnonzero(float a)
	/*< sign of a float, if a==0 then gives a value of 1. >*/
{
	float b;
	if (a>=0.)      b = 1.;
	else b =-1.;
	return b;
}

int compare (const void * a, const void * b)
{
	float fa = *(const float*) a;
	float fb = *(const float*) b;
	return (fa > fb) - (fa < fb);
}

int omp_thread_count() {
	int n = 0;
#pragma omp parallel reduction(+:n)
	n += 1;
	return n;
}

void boundary_condition(complex *d_x,int nmx,int lmx,int nmy,int lmy)
{
	int imx,imy;
	float tmx,tmy; 
	tmx = 1.;tmy = 1.;
	for (imx=0;imx<nmx;imx++){
		if (imx>=0   && imx<lmx) tmx = expf(-powf(0.015*((float) lmx - imx),2.));
		if (imx>=lmx && imx<=nmx-lmx) tmx = 1.;
		if (imx>nmx-lmx && imx<nmx) tmx = expf(-powf(0.015*((float) imx - nmx + lmx),2.));
		for (imy=0;imy<nmy;imy++){
			if (imy>=0   && imy<lmy) tmy = expf(-powf(0.015*((float) lmy - imy),2.));
			if (imy>=lmy && imy<=nmy-lmy) tmy = 1.;
			if (imy>nmy-lmy && imy<nmy) tmy = expf(-powf(0.015*((float) imy - nmy + lmy),2.));
			d_x[imx*nmy + imy] *= tmx*tmy;
		}
	}  

	return;
}

void elastic_separate_3d(complex *ux, complex *uy, complex *uz,
		complex *up, complex *us1, complex *us2,
		float w, 
		float dkx, int nkx, int nmx, float omx, float dmx,
		float dky, int nky, int nmy, float omy, float dmy,
		float vp,float vs,
		fftwf_plan p1,fftwf_plan p2,
		bool sep, bool adj)
{
	// elastic wavefield decomposition / recomposition with using the average medium 
	// parameters for a given layer. 

	int imx,imy,ikx,iky,ik;
	fftwf_complex *a,*b;  
	complex *ux_k,*uz_k,*up_k,*us2_k;
	float kx,kxp,kxs,kzp,kzs,sp,ss,norm,norm_p,norm_s;
	ux_k = alloc1complex(nkx*nky);
	//uy_k = alloc1complex(nkx*nky);
	uz_k = alloc1complex(nkx*nky);
	up_k = alloc1complex(nkx*nky);
	//us1_k = alloc1complex(nkx*nky);
	us2_k = alloc1complex(nkx*nky);
	a  = fftwf_malloc(sizeof(fftwf_complex) * nkx*nky);
	b  = fftwf_malloc(sizeof(fftwf_complex) * nkx*nky);
	if (sep){ /* separation of wavefield components to wavefield potentials */
		// x-component
		for(imx=0;imx<nkx;imx++){ 
			for(imy=0;imy<nky;imy++){ 
				a[imx*nky + imy] = (imx < nmx && imy < nmy) ? ux[imx*nmy + imy] : 0.;
			}
		}
		fftwf_execute_dft(p1,a,a);
		for(ik=0;ik<nkx*nky;ik++) ux_k[ik] = a[ik]/sqrtf(nkx*nky);
		// y-component
		//for(imx=0;imx<nkx;imx++){ 
		//	for(imy=0;imy<nky;imy++){ 
		//		a[imx*nky + imy] = (imx < nmx && imy < nmy) ? uy[imx*nmy + imy] : 0.;
		//	}
		//}
		//fftwf_execute_dft(p1,a,a);
		//for(ik=0;ik<nkx*nky;ik++) uy_k[ik] = a[ik]/sqrtf(nkx*nky);
		// z-component
		for(imx=0;imx<nkx;imx++){ 
			for(imy=0;imy<nky;imy++){ 
				a[imx*nky + imy] = (imx < nmx && imy < nmy) ? uz[imx*nmy + imy] : 0.;
			}
		}
		fftwf_execute_dft(p1,a,a);
		for(ik=0;ik<nkx*nky;ik++) uz_k[ik] = a[ik]/sqrtf(nkx*nky); 
		for (ikx=0;ikx<nkx;ikx++){
			kx = ikx<nkx/2. ? dkx*ikx : -(dkx*nkx - dkx*ikx);
			for (iky=0;iky<nky;iky++){
				sp = w*w/(vp*vp) - kx*kx;
				ss = w*w/(vs*vs) - kx*kx;
				kzp = sp > 0. ? sqrtf(sp) : 0.;
				kzs = ss > 0. ? sqrtf(ss) : 0.;
				norm_p = sqrtf(kx*kx + kzp*kzp);
				norm_s = sqrtf(kx*kx + kzs*kzs);
				kxp = kx/norm_p;
				kxs = kx/norm_s;
				kzp = kzp/norm_p;
				kzs = kzs/norm_s;
				norm = kxp*kxs + kzp*kzs;
				if (norm >= 0.01 && w >= 10.){
					if (!adj){
						up_k[ikx*nky + iky]  = ( kxs*ux_k[ikx*nky + iky] + kzs*uz_k[ikx*nky + iky])/norm;  
						us2_k[ikx*nky + iky] = (-kzp*ux_k[ikx*nky + iky] + kxp*uz_k[ikx*nky + iky])/norm;
					}
					else{
						up_k[ikx*nky + iky]  =  kxp*ux_k[ikx*nky + iky] + kzp*uz_k[ikx*nky + iky];  
						us2_k[ikx*nky + iky] = -kzs*ux_k[ikx*nky + iky] + kxs*uz_k[ikx*nky + iky];					
					}
				}
				else{
					up_k[ikx*nky + iky]  = uz_k[ikx*nky + iky];  
					us2_k[ikx*nky + iky] = ux_k[ikx*nky + iky];
				}
			}
		}
		for (ik=0;ik<nkx*nky;ik++) b[ik] = up_k[ik];
		fftwf_execute_dft(p2,b,b);
		for(imx=0; imx<nkx;imx++){ 
			for(imy=0; imy<nky;imy++){ 
				if (imx < nmx && imy < nmy) up[imx*nmy + imy] = b[imx*nky + imy]/sqrtf(nkx*nky);
			}
		}      
		//for (ik=0;ik<nkx*nky;ik++) b[ik] = us1_k[ik];
		//fftwf_execute_dft(p2,b,b);
		//for(imx=0; imx<nkx;imx++){ 
		//	for(imy=0; imy<nky;imy++){ 
		//		if (imx < nmx && imy < nmy) us1[imx*nmy + imy] = b[imx*nky + imy]/sqrtf(nkx*nky);
		//	}
		//}      
		for (ik=0;ik<nkx*nky;ik++) b[ik] = us2_k[ik];
		fftwf_execute_dft(p2,b,b);
		for(imx=0; imx<nkx;imx++){ 
			for(imy=0; imy<nky;imy++){ 
				if (imx < nmx && imy < nmy) us2[imx*nmy + imy] = b[imx*nky + imy]/sqrtf(nkx*nky);
			}
		}      
	}
	else { /* combination of wavefield potentials to wavefield components */
		// p-component
		for(imx=0;imx<nkx;imx++){ 
			for(imy=0;imy<nky;imy++){ 
				a[imx*nky + imy] = (imx < nmx && imy < nmy) ? up[imx*nmy + imy] : 0.;
			}
		}
		fftwf_execute_dft(p1,a,a);
		for(ik=0;ik<nkx*nky;ik++) up_k[ik] = a[ik]/sqrtf(nkx*nky);
		// s1-component
		//for(imx=0;imx<nkx;imx++){ 
		//	for(imy=0;imy<nky;imy++){ 
		//		a[imx*nky + imy] = (imx < nmx && imy < nmy) ? us1[imx*nmy + imy] : 0.;
		//	}
		//}
		//fftwf_execute_dft(p1,a,a);
		//for(ik=0;ik<nkx*nky;ik++) us1_k[ik] = a[ik]/sqrtf(nkx*nky);
		// s2-component
		for(imx=0;imx<nkx;imx++){ 
			for(imy=0;imy<nky;imy++){ 
				a[imx*nky + imy] = (imx < nmx && imy < nmy) ? us2[imx*nmy + imy] : 0.;
			}
		}
		fftwf_execute_dft(p1,a,a);
		for(ik=0;ik<nkx*nky;ik++) us2_k[ik] = a[ik]/sqrtf(nkx*nky); 
		for (ikx=0;ikx<nkx;ikx++){
			kx = ikx<nkx/2. ? dkx*ikx : -(dkx*nkx - dkx*ikx);
			for (iky=0;iky<nky;iky++){
				sp = w*w/(vp*vp) - kx*kx;
				ss = w*w/(vs*vs) - kx*kx;
				kzp = sp > 0. ? sqrtf(sp) : 0.;
				kzs = ss > 0. ? sqrtf(ss) : 0.;
				norm_p = sqrtf(kx*kx + kzp*kzp);
				norm_s = sqrtf(kx*kx + kzs*kzs);
				kxp = kx/norm_p;
				kxs = kx/norm_s;
				kzp = kzp/norm_p;
				kzs = kzs/norm_s;
				norm = kxp*kxs + kzp*kzs;
				if (norm >= 0.01 && w >= 10.){
					if (!adj){
						ux_k[ikx*nky + iky] = kxp*up_k[ikx*nky + iky] - kzs*us2_k[ikx*nky + iky];
						uz_k[ikx*nky + iky] = kzp*up_k[ikx*nky + iky] + kxs*us2_k[ikx*nky + iky];
					}
					else{
						ux_k[ikx*nky + iky] = (kxs*up_k[ikx*nky + iky] - kzp*us2_k[ikx*nky + iky])/norm;
						uz_k[ikx*nky + iky] = (kzs*up_k[ikx*nky + iky] + kxp*us2_k[ikx*nky + iky])/norm;
					}
				}
				else{
					ux_k[ikx*nky + iky] = us2_k[ikx*nky + iky];
					uz_k[ikx*nky + iky] = up_k[ikx*nky + iky];  
				}
			}
		}
		for (ik=0;ik<nkx*nky;ik++) b[ik] = ux_k[ik];
		fftwf_execute_dft(p2,b,b);
		for(imx=0; imx<nkx;imx++){ 
			for(imy=0; imy<nky;imy++){ 
				if (imx < nmx && imy < nmy) ux[imx*nmy + imy] = b[imx*nky + imy]/sqrtf(nkx*nky);
			}
		}      
		//for (ik=0;ik<nkx*nky;ik++) b[ik] = uy_k[ik];
		//fftwf_execute_dft(p2,b,b);
		//for(imx=0; imx<nkx;imx++){ 
		//	for(imy=0; imy<nky;imy++){ 
		//		if (imx < nmx && imy < nmy) uy[imx*nmy + imy] = b[imx*nky + imy]/sqrtf(nkx*nky);
		//	}
		//}      
		for (ik=0;ik<nkx*nky;ik++) b[ik] = uz_k[ik];
		fftwf_execute_dft(p2,b,b);
		for(imx=0; imx<nkx;imx++){ 
			for(imy=0; imy<nky;imy++){ 
				if (imx < nmx && imy < nmy) uz[imx*nmy + imy] = b[imx*nky + imy]/sqrtf(nkx*nky);
			}
		}      
	}

	fftwf_free(a);
	fftwf_free(b);
	free1complex(ux_k);
	//free1complex(uy_k);
	free1complex(uz_k);
	free1complex(up_k);
	//free1complex(us1_k);
	free1complex(us2_k);

	return;
}

void elastic_separate_3d_window(complex *ux, complex *uy, complex *uz,
		complex *up, complex *us1, complex *us2,
		float w, 
		float dkx, int nkx, int nmx, float omx, float dmx,
		float dky, int nky, int nmy, float omy, float dmy,
		float **vp,float **vs, int iz, int Lw,
		fftwf_plan p1,fftwf_plan p2,
		bool sep, bool adj)
{
	int ix,Nw,iw,iwi,iwf;
	float vp0,vs0;
	complex *W1 = alloc1complex(nmx);
	complex *W2 = alloc1complex(nmx);
	complex *W3 = alloc1complex(nmx);
	for (ix=0;ix<nmx;ix++) W1[ix] = 0.;
	for (ix=0;ix<nmx;ix++) W2[ix] = 0.;
	for (ix=0;ix<nmx;ix++) W3[ix] = 0.;
	Nw = (int) nmx/Lw + 1;
	Nw = (Nw-1)*Lw >= nmx ? Nw - 1 : Nw;
	//fprintf(stderr,"nmx=%d,Lw=%d,Nw=%d\n",nmx,Lw,Nw);
	for (iw=0;iw<Nw;iw++){
		iwi = iw*Lw;
		iwf = iwi + Lw < nmx ? iwi + Lw : nmx;
		vp0 = 0.0; vs0 = 0.0;
		for (ix=iwi;ix<iwf;ix++){
			vp0 += vp[ix][iz];
			vs0 += vs[ix][iz];
		}
		vp0 /= (iwf-iwi+1);
		vs0 /= (iwf-iwi+1);
		//fprintf(stderr,"iwi=%d,iwf=%d,vp0=%f,vs0=%f\n",iwi,iwf,vp0,vs0);
		if (sep){
			elastic_separate_3d(ux,uy,uz,W1,W2,W3,w,dkx,nkx,nmx,omx,dmx,dky,nky,nmy,omy,dmy,vp0,vs0,p1,p2,true,adj);
			for (ix=iwi;ix<iwf;ix++){
				up[ix]  = W1[ix];	
				us1[ix] = W2[ix];	
				us2[ix] = W3[ix];	
			}
		}
		else{
			elastic_separate_3d(W1,W2,W3,up,us1,us2,w,dkx,nkx,nmx,omx,dmx,dky,nky,nmy,omy,dmy,vp0,vs0,p1,p2,false,adj);
			for (ix=iwi;ix<iwf;ix++){
				ux[ix] = W1[ix];	
				uy[ix] = W2[ix];	
				uz[ix] = W3[ix];	
			}
		}
	}

	free1complex(W1);
	free1complex(W2);
	free1complex(W3);
	
	return;
}

