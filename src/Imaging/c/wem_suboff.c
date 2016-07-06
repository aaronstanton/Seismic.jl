#include "seismic.h"
#include "wem_suboff.h"
void wem_suboff(float **d, float **m, float **wav,
		int nt, float ot, float dt, 
		int nmx,float omx, float dmx,
		int nhx,float ohx, float dhx,
		float sx,
		int nz, float oz, float dz, float gz, float sz,
		float **vel, int nref, float fmin, float fmax,
		int padt, int padx,
		bool adj, bool pspi, bool verbose,
		float kz_eps)
/*< wave equation depth migration operator. >*/
{
	int iz,ix,imx,ihx,igx,ik,iw,it,nw,nkx,ntfft;
	float dw,dkx;
	int ifmin,ifmax;
	float *d_t;
	complex *d_w,**d_g_wx,**d_s_wx;
	fftwf_complex *a,*b;
	int *n;
	fftwf_plan p1,p2;
	float *po,**pd;
	float progress;
	int ithread,nthread;
	float **m_threads;
	float **vref,vmin,vmax,v;
	int **iref1,**iref2;
	int iref;

	/* decompose slowness into layer average, and layer purturbation */
	po = alloc1float(nz); 
	pd = alloc2float(nz,nmx); 
	for (iz=0;iz<nz;iz++){
		po[iz] = 0.;
		for (ix=0;ix<nmx;ix++) po[iz] += vel[ix][iz];
		po[iz] /= (float) nmx;
		po[iz]  = 1./po[iz];
		for (ix=0;ix<nmx;ix++) pd[ix][iz] = 1.0/vel[ix][iz] - po[iz];
	}

	/****************************************************************************************/
	/* generate reference velocities for each depth step */
	vref = alloc2float(nz,nref); /* reference velocities for each layer */
	iref1 = alloc2int(nz,nmx); /* index of nearest lower reference velocity for each subsurface point */
	iref2 = alloc2int(nz,nmx); /* index of nearest upper reference velocity for each subsurface point */
	for (iz=0;iz<nz;iz++){
		vmin=vel[0][iz];
		for (ix=0;ix<nmx;ix++) if (vel[ix][iz] < vmin) vmin = vel[ix][iz];
		vmax=vel[nmx-1][iz];
		for (ix=0;ix<nmx;ix++) if (vel[ix][iz] > vmax) vmax = vel[ix][iz];

		for (iref=0;iref<nref;iref++) vref[iref][iz] = vmin + (float) iref*(vmax-vmin)/((float) nref-1);
		for (ix=0;ix<nmx;ix++){
			v = vel[ix][iz];
			if (vmax>vmin+10){
				iref = (int) (nref-1)*(v-vmin)/(vmax-vmin);
				iref1[ix][iz] = iref;
				iref2[ix][iz] = iref+1;
				if (iref>nref-2){
					iref1[ix][iz] = nref-1;
					iref2[ix][iz] = nref-1;
				}
			}
			else{
				iref1[ix][iz] = 0;
				iref2[ix][iz] = 0;
			}
		}
	}
	/****************************************************************************************/

	if (adj){
		for (ix=0;ix<nmx*nhx;ix++) for (iz=0;iz<nz;iz++) m[ix][iz] = 0.;
	}
	else{
		for (ix=0;ix<nmx;ix++) for (it=0;it<nt;it++) d[ix][it] = 0.;
	}
	ntfft = (int) 2*padt*((float) nt)/2;
	nw = (int) ntfft/2 + 1;
	nkx = nmx > 1 ? padx*nmx : 1;
	dkx = 2.0*PI/((float) nkx)/dmx;
	dw = 2.0*PI/((float) ntfft)/dt;

	if(fmax*dt*ntfft+1<nw) ifmax = (int) fmax*dt*ntfft + 1;
	else ifmax = nw;
	if(fmin*dt*ntfft+1<ifmax) ifmin = (int) fmin*dt*ntfft;
	else ifmin = 0;
	d_g_wx = alloc2complex(nw,nmx);
	d_s_wx = alloc2complex(nw,nmx);
	d_t = alloc1float(nt);
	d_w = alloc1complex(nw);
	for (it=0;it<nt;it++)  d_t[it] = 0.;  
	for (iw=0;iw<nw;iw++)  d_w[iw] = 0.; 

	/* set up fftw plans and pass them to the OMP region of the code */
	a  = fftwf_malloc(sizeof(fftwf_complex) * nkx);
	b  = fftwf_malloc(sizeof(fftwf_complex) * nkx);
	n = alloc1int(1); 
	n[0] = nkx;
	p1 = fftwf_plan_dft(1, n, a, a, FFTW_FORWARD, FFTW_ESTIMATE);
	p2 = fftwf_plan_dft(1, n, b, b, FFTW_BACKWARD, FFTW_ESTIMATE);
	for (ik=0;ik<nkx;ik++){
		a[ik] = 0.;
		b[ik] = 0.;
	} 
	fftwf_execute_dft(p1,a,a);
	fftwf_execute_dft(p2,b,b);
	/**********************************************************************/
	igx = (int) (sx - omx)/dmx; /*position to inject source in x-dir*/
	/* source wavefield*/
	for (ix=0;ix<nmx;ix++) for (iw=0;iw<nw;iw++) d_s_wx[ix][iw] = 0.;
	for (it=0;it<nt;it++) d_t[it] = wav[0][it];
	f_op(d_w,d_t,nw,nt,1); /* d_t to d_w */
	for (iw=0;iw<nw;iw++) d_s_wx[igx][iw] = d_w[iw];
	/* receiver wavefield*/
	if (adj){
		for (ix=0;ix<nmx;ix++){
			for (it=0;it<nt;it++) d_t[it] = d[ix][it];
			f_op(d_w,d_t,nw,nt,1); /* d_t to d_w */
			for (iw=0;iw<ifmin;iw++) d_g_wx[ix][iw] = 0.;
			for (iw=ifmin;iw<ifmax;iw++){ 
				d_g_wx[ix][iw] = d_w[iw];
			}
			for (iw=ifmax;iw<nw;iw++) d_g_wx[ix][iw] = 0.;
		}
	}
	else{
		for (ix=0;ix<nmx;ix++){
			for (iw=0;iw<nw;iw++){
				d_g_wx[ix][iw] = 0.;
			}
		}
	}
	nthread = omp_thread_count();
	if (adj){
		m_threads = alloc2float(nz,nmx*nhx*nthread);
		for (imx=0;imx<nmx;imx++){
			for (ihx=0;ihx<nhx;ihx++){
				for (ithread=0;ithread<nthread;ithread++){
					for (iz=0;iz<nz;iz++){
						m_threads[imx*nhx*nthread + ihx*nthread + ithread][iz] = 0.;
					}
				}
			}
		}
	}
	else{
		m_threads = alloc2float(nz,nmx*nhx);
		for (imx=0;imx<nmx;imx++){
			for (ihx=0;ihx<nhx;ihx++){
				for (iz=0;iz<nz;iz++){
					m_threads[imx*nhx + ihx][iz] = m[imx*nhx + ihx][iz];
				}
			}
		}
	}
	progress = 0.;
#pragma omp parallel for private(iw) shared(m_threads,d_g_wx,d_s_wx)
	for (iw=ifmin;iw<ifmax;iw++){ 
		progress += 1./((float) ifmax - ifmin);
		if (verbose) progress_msg(progress);
		extrap1f(m_threads,d_g_wx,d_s_wx,iw,ifmax,nw,ifmax,ntfft,dw,dkx,nkx,nz,oz,dz,gz,sz,nmx,omx,dmx,nhx,ohx,dhx,nthread,vel,po,pd,vref,iref1,iref2,nref,p1,p2,adj,pspi,verbose,kz_eps);
	}
	if (verbose) fprintf(stderr,"\n");
	if (adj){
		// reduction over parallel axis 
		for (imx=0;imx<nmx;imx++) for (ihx=0;ihx<nhx;ihx++) for (ithread=0;ithread<nthread;ithread++) for (iz=0;iz<nz;iz++) m[imx*nhx + ihx][iz] += m_threads[imx*nhx*nthread + ihx*nthread + ithread][iz];
	}
	else{
		for (ix=0;ix<nmx;ix++){
			for (iw=0;iw<ifmin;iw++) d_w[iw] = 0.;
			for (iw=ifmin;iw<ifmax;iw++){ 
				d_w[iw] = d_g_wx[ix][iw];
			}
			for (iw=ifmax;iw<nw;iw++) d_w[iw] = 0.;
			f_op(d_w,d_t,nw,nt,0); /* d_w to d_t */
			for (it=0;it<nt;it++) d[ix][it] = d_t[it];
		}
	}

	free1int(n); 
	fftwf_free(a);
	fftwf_free(b);  
	fftwf_destroy_plan(p1);
	fftwf_destroy_plan(p2);  
	free1float(d_t);
	free1complex(d_w);
	free2float(m_threads);
	free2complex(d_g_wx);
	free2complex(d_s_wx);
	free2float(vref);
	free2int(iref1);
	free2int(iref2);
	return;
} 

void extrap1f(float **m,complex **d_g_wx, complex **d_s_wx,
		int iw, int ang_iw_max, int nw,int ifmax,int ntfft,float dw,float dkx,int nkx,
		int nz, float oz, float dz, float gz, float sz,
		int nmx,float omx, float dmx,
		int nhx,float ohx, float dhx,
		int nthread,
		float **v,float *po,float **pd,
		float **vref, int **iref1, int **iref2, int nref,
		fftwf_plan p1,fftwf_plan p2,
		bool adj, bool pspi, bool verbose,
		float kz_eps)   
/*< extrapolate 1 frequency >*/
{
	float w,factor,z,hx,sx,gx;
	int iz,ix,imx,ihx,isx,igx,ithread;
	complex *d_xg,*d_xs,**smig;
	ithread = omp_get_thread_num(); 
	d_xg = alloc1complex(nmx);
	d_xs = alloc1complex(nmx);
	for (ix=0;ix<nmx;ix++) d_xg[ix] = 0.;
	for (ix=0;ix<nmx;ix++) d_xs[ix] = 0.;
	if (iw==0) factor = 1.;
	else factor = 2.;
	w = iw*dw;
	if (adj){
		for (ix=0;ix<nmx;ix++){ 
			d_xs[ix] = d_s_wx[ix][iw]/sqrtf((float) ntfft);
			d_xg[ix] = powf(w,2)*d_g_wx[ix][iw]/sqrtf((float) ntfft); //powf(w,2)*
		}
		for (iz=0;iz<nz;iz++){ // extrapolate source and receiver wavefields
			z = oz + dz*iz;
			if (z >= sz){
				if (pspi) pspiop(d_xs,w,dkx,nkx,nmx,omx,dmx,dz,iz,v,po,pd,vref,iref1,iref2,nref,p1,p2,true,true,verbose,kz_eps);
				else ssop(d_xs,w,dkx,nkx,nmx,omx,dmx,dz,iz,v,po,pd,p1,p2,true,true,verbose,kz_eps);
			} 
			if (z >= gz){
				if (pspi) pspiop(d_xg,w,dkx,nkx,nmx,omx,dmx,dz,iz,v,po,pd,vref,iref1,iref2,nref,p1,p2,true,false,verbose,kz_eps);
				else ssop(d_xg,w,dkx,nkx,nmx,omx,dmx,dz,iz,v,po,pd,p1,p2,true,false,verbose,kz_eps);
				for (imx=0;imx<nmx;imx++){ 
					for (ihx=0;ihx<nhx;ihx++){
						hx = ihx*dhx + ohx;
						sx = (imx*dmx + omx) - hx;
						gx = (imx*dmx + omx) + hx;
						isx = (int) truncf((sx - omx)/dmx);
						igx = (int) truncf((gx - omx)/dmx);
						if (isx >=0 && isx < nmx && igx >=0 && igx < nmx){
							m[imx*nhx*nthread + ihx*nthread + ithread][iz] += factor*crealf(d_xg[igx]*conjf(d_xs[isx]));
						}
					}
				}
			}
		}
	}
	else{
		smig = alloc2complex(nz,nmx);
		for (ix=0;ix<nmx;ix++) d_xs[ix] = d_s_wx[ix][iw]/sqrtf((float) ntfft);
		for (iz=0;iz<nz;iz++){ // extrapolate source wavefield 
			z = oz + dz*iz;
			if (z >= sz){
				if (pspi) pspiop(d_xs,w,dkx,nkx,nmx,omx,dmx,dz,iz,v,po,pd,vref,iref1,iref2,nref,p1,p2,true,true,verbose,kz_eps);
				else ssop(d_xs,w,dkx,nkx,nmx,omx,dmx,dz,iz,v,po,pd,p1,p2,true,true,verbose,kz_eps);
				for (ix=0;ix<nmx;ix++) smig[ix][iz] = d_xs[ix];
			}
			else{
				for (ix=0;ix<nmx;ix++) smig[ix][iz] = 0.;
			}
		}
		for (ix=0;ix<nmx;ix++) d_xg[ix] = 0.;
		for (iz=nz-1;iz>=0;iz--){ // extrapolate receiver wavefield 
			z = oz + dz*iz;
			if (z >= gz){
				for (imx=0;imx<nmx;imx++){ 
					for (ihx=0;ihx<nhx;ihx++){
						hx = ihx*dhx + ohx;
						sx = (imx*dmx + omx) - hx;
						gx = (imx*dmx + omx) + hx;
						isx = (int) truncf((sx - omx)/dmx);
						igx = (int) truncf((gx - omx)/dmx);
						if (isx >=0 && isx < nmx && igx >=0 && igx < nmx){
							d_xg[igx] = d_xg[igx] + smig[isx][iz]*m[imx*nhx*nthread + ihx*nthread + ithread][iz];
						}
					}
				}
			}
			if (pspi) pspiop(d_xg,w,dkx,nkx,nmx,omx,dmx,dz,iz,v,po,pd,vref,iref1,iref2,nref,p1,p2,false,false,verbose,kz_eps);
			else ssop(d_xg,w,dkx,nkx,nmx,omx,dmx,dz,iz,v,po,pd,p1,p2,false,false,verbose,kz_eps);
		}
		for (ix=0;ix<nmx;ix++){
			d_g_wx[ix][iw] = powf(w,2)*d_xg[ix]/sqrtf((float) ntfft);// powf(w,2)*			
		}
		free2complex(smig);
	}

	free1complex(d_xs);
	free1complex(d_xg);

	return;
}

void ssop(complex *d_x,
		float w,float dkx,int nkx,int nmx,float omx,float dmx,float dz,int iz,
		float **v,float *po,float **pd,
		fftwf_plan p1,fftwf_plan p2,
		bool adj, 
		bool src,
		bool verbose,
		float kz_eps)
{
	float kx;
	complex L;
	int ik,ikx,imx; 
	complex *d_k;
	fftwf_complex *a,*b;
	complex w2,s;	

	int lmx;
	if (nmx>100) lmx=30;
	else lmx=0;
	a  = fftwf_malloc(sizeof(fftwf_complex) * nkx);
	b  = fftwf_malloc(sizeof(fftwf_complex) * nkx);
	d_k = alloc1complex(nkx);
	w2 = kz_eps + w*I;
	w2 *= w2;
	w2 *= po[iz];
	w2 *= po[iz];
	if (adj){
		for(imx=0; imx<nkx;imx++){ 
			if (imx < nmx) a[imx] = d_x[imx];
			else a[imx] = 0.;
		}
		fftwf_execute_dft(p1,a,a); 
		for (ikx=0;ikx<nkx;ikx++){
			if (ikx<= (int) nkx/2) kx = (float) dkx*ikx;
			else                   kx = -((float) dkx*nkx - dkx*ikx);
			s = w2 + powf(kx,2);
			L = cexp(-csqrtf(s)*dz);
			if (src) d_k[ikx] = ((complex) a[ikx])*L/sqrtf((float) nkx);
			else     d_k[ikx] = ((complex) a[ikx])*conjf(L)/sqrtf((float) nkx);   
		}
		for(ik=0; ik<nkx;ik++) b[ik] = (fftwf_complex) d_k[ik];
		fftwf_execute_dft(p2,b,b);
		for(imx=0; imx<nkx;imx++){ 
			if (imx < nmx){
				L = cexpf(I*w*pd[imx][iz]*dz);
				if (src) d_x[imx] = ((complex) b[imx])*conjf(L)/sqrtf((float) nkx); // SS operator
				else     d_x[imx] = ((complex) b[imx])*L/sqrtf((float) nkx); // SS operator
			}
		}
		boundary_condition(d_x,nmx,lmx);    
	}
	else{
		boundary_condition(d_x,nmx,lmx);    
		for(imx=0; imx<nkx;imx++){ 
			if (imx < nmx){
				L = cexpf(I*w*pd[imx][iz]*dz);  
				a[imx] = d_x[imx]*conjf(L); // SS operator
			}
			else a[imx] = 0.0;
		}
		fftwf_execute_dft(p1,a,a); 
		for (ikx=0;ikx<nkx;ikx++){
			if (ikx<= (int) nkx/2) kx = (float) dkx*ikx;
			else                   kx = -((float) dkx*nkx - dkx*ikx);
			s = w2 + powf(kx,2);
			L = cexp(-csqrtf(s)*dz);
			d_k[ikx] = ((complex) a[ikx])*L/sqrtf((float) nkx);
		}
		for(ik=0; ik<nkx;ik++) b[ik] = (fftwf_complex) d_k[ik];
		fftwf_execute_dft(p2,b,b);
		for(imx=0; imx<nkx;imx++){ 
			if (imx < nmx){ 
				d_x[imx] = ((complex) b[imx])/sqrtf((float) nkx);
			}
		}
	}
	free1complex(d_k);
	fftwf_free(a);
	fftwf_free(b);

	return;
}

void pspiop(complex *d_x,
		float w,float dkx,int nkx,
		int nmx,float omx,float dmx,
		float dz,int iz,
		float **vel,float *po,float **pd,
		float **vref, int **iref1, int **iref2, int nref,
		fftwf_plan p1,fftwf_plan p2,
		bool adj, 
		bool src,
		bool verbose,
		float kz_eps)
{


	float kx;
	complex L;
	int ik,ikx,imx,ix; 
	complex *d_k;
	fftwf_complex *a,*b;
	int lmx;
	complex **dref;
	complex w2,s;	
	int iref;
	float v,vref1,vref2,aa,pd_ref;
	if (nmx>100) lmx=30;
	else lmx=0;
	a  = fftwf_malloc(sizeof(fftwf_complex) * nkx);
	b  = fftwf_malloc(sizeof(fftwf_complex) * nkx);
	dref = alloc2complex(nref,nmx);
	d_k = alloc1complex(nkx);

	if (adj){
		for(imx=0; imx<nkx;imx++){ 
			if (imx < nmx) a[imx] = d_x[imx];
			else a[imx] = 0.;
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
				s = w2 + powf(kx,2);
				L = cexp(-csqrtf(s)*dz);
				if (src) d_k[ikx] = ((complex) a[ikx])*L/sqrtf((float) nkx);
				else     d_k[ikx] = ((complex) a[ikx])*conjf(L)/sqrtf((float) nkx);
			}
			for(ik=0; ik<nkx;ik++) b[ik] = (fftwf_complex) d_k[ik];
			fftwf_execute_dft(p2,b,b);
			for(imx=0; imx<nkx;imx++){ 
				if (imx < nmx){
					pd_ref = (1.0/vel[imx][iz]) - (1.0/vref[iref][iz]);
					L = cexpf(I*w*pd_ref*dz);
					if (src) dref[imx][iref] = ((complex) b[imx])*conjf(L)/sqrtf((float) nkx);
					else     dref[imx][iref] = ((complex) b[imx])*L/sqrtf((float) nkx);
				}
			}
		}
		for (ix=0;ix<nmx;ix++){
			v = vel[ix][iz];
			vref1 = vref[iref1[ix][iz]][iz];
			vref2 = vref[iref2[ix][iz]][iz];
			aa = linear_interp(vref1,vref2,v);
			d_x[ix] = aa*dref[ix][iref1[ix][iz]] + (1.0-aa)*dref[ix][iref2[ix][iz]];
		}
		boundary_condition(d_x,nmx,lmx);
	}
	else{
		boundary_condition(d_x,nmx,lmx);    
		for (ix=0;ix<nmx;ix++) for (iref=0;iref<nref;iref++) dref[ix][iref] = 0.0;
		for (ix=0;ix<nmx;ix++){
			v = vel[ix][iz];
			vref1 = vref[iref1[ix][iz]][iz];
			vref2 = vref[iref2[ix][iz]][iz];
			aa = linear_interp(vref1,vref2,v);
			dref[ix][iref1[ix][iz]] = aa*d_x[ix];
			dref[ix][iref2[ix][iz]] = (1.0-aa)*d_x[ix];
		}
		for (ix=0;ix<nmx;ix++) d_x[ix] = 0.;
		for (iref=0;iref<nref;iref++){
			w2 = 0.5f + w*I;
			w2 *= w2;
			w2 /= vref[iref][iz];
			w2 /= vref[iref][iz];
			for(imx=0; imx<nkx;imx++){
				if (imx < nmx){ 
					pd_ref = (1.0/vel[imx][iz]) - (1.0/vref[iref][iz]);
					L = cexpf(I*w*pd_ref*dz);  
					a[imx] = dref[imx][iref]*conjf(L);
				}
				else a[imx] = 0.;
			}
			fftwf_execute_dft(p1,a,a);
			for (ikx=0;ikx<nkx;ikx++){
				if (ikx<= (int) nkx/2) kx = (float) dkx*ikx;
				else                   kx = -((float) dkx*nkx - dkx*ikx);
				s = w2 + powf(kx,2);
				L = cexp(-csqrtf(s)*dz);
				d_k[ikx] = ((complex) a[ikx])*L/sqrtf((float) nkx);
			}
			for(ik=0; ik<nkx;ik++) b[ik] = (fftwf_complex) d_k[ik];
			fftwf_execute_dft(p2,b,b);
			for(imx=0; imx<nkx;imx++){ 
				if (imx < nmx){
					d_x[imx] += ((complex) b[imx])/sqrtf((float) nkx);
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

void boundary_condition(complex *d_x,int nmx,int lmx)
{
	int imx;
	float tmx; 
	tmx = 1.;
	for (imx=0;imx<nmx;imx++){
		if (imx>=0   && imx<lmx) tmx = expf(-powf(0.015*((float) lmx - imx),2.));
		if (imx>=lmx && imx<=nmx-lmx) tmx = 1.;
		if (imx>nmx-lmx && imx<nmx) tmx = expf(-powf(0.015*((float) imx - nmx + lmx),2.));
		d_x[imx] *= tmx;
	}  

	return;
}

float signf1(float a)
	/*< sign of a float, always has an amplitude of 1 >*/
{
	float b;
	if (a>=0.) b = 1.;
	else b =-1.;
	return b;
}
