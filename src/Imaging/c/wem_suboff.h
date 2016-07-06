#ifndef _wem_suboff_h_
#define _wem_suboff_h_
void wem_suboff(float **d, float **m, float **wav,
		int nt, float ot, float dt, 
		int nmx,float omx, float dmx,
		int nhx,float ohx, float dhx,
		float sx,
		int nz, float oz, float dz, float gz, float sz,
		float **vel, int nref, float fmin, float fmax,
		int padt, int padx,
		bool adj, bool pspi, bool verbose,
		float kz_eps);
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
              float kz_eps);
void ssop(complex *d_x,
          float w,float dkx,int nkx,int nmx,float omx,float dmx,float dz,int iz,
          float **v,float *po,float **pd,
          fftwf_plan p1,fftwf_plan p2,
          bool adj, 
          bool src,
          bool verbose,
          float kz_eps);
void pspiop(complex *d_x,
		float w,float dkx,int nkx,int nmx,float omx,float dmx,float dz,int iz,
		float **vel,float *po,float **pd,
		float **vref, int **iref1, int **iref2, int nref,
		fftwf_plan p1,fftwf_plan p2,
		bool adj, 
		bool src,
		bool verbose,
        float kz_eps);
float linear_interp(float x1,float x2,float x);
void f_op(complex *m,float *d,int nw,int nt,bool adj);
void progress_msg(float progress);
float signf(float a);
int compare (const void * a, const void * b);
int omp_thread_count();
void boundary_condition(complex *d_x,int nmx,int lmx);
float signf1(float a);            
#endif
