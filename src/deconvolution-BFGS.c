#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multimin.h>

#define MaxS 35000
#define Nx 6
#define Npar 2

double flow[35000];
double tm[35000];
long nf;
double tmax;

void fgetstr(FILE *f, char *s)
{
 long n=0;
 while(n<MaxS && !feof(f)) 
  {
   s[n]=fgetc(f);
   if(s[n]=='\n') break; 
   if(!feof(f)) n++;
  }
 s[n]='\0';
}

long read_hf(char *fname)
{
 char s[MaxS];
 FILE *f;
 float t,T,hf,h,nh;

 f=fopen(fname,"r");
 fgetstr(f,s);
 nf=0;
 while(!feof(f))
  {
   fgetstr(f,s);
   if(s==NULL) continue;
   sscanf(s,"%f,%f,%f,%f,%f",&t,&T,&hf,&h,&nh);
   tmax=t;flow[nf]=nh;tm[nf]=t;
   nf++;
  }
 fclose(f);
// printf("n=%ld   tmax=%f\n",nf,tmax);
 return nf;
}

double test_kinetic(double k1, double A01, double C01, double H1, double k2, double A02, double C02, double H2)
{
 double m1, A1, C1, dC1, P1;
 double m2, A2, C2, dC2, P2;
 double t,dt,Psum,dQ;
 double Q;
 long n;

 m1=1.0;m2=1.0;

// printf("t, k1=%f k2=%f A01=%f A02=%f C01=%f C02=%f, 0, 0, Hf1, Hf2, Hf\n",k1,k2,A01,A02,C01,C02);
 C1=C01;
 C2=C02;
 dt=tmax/nf;
 Q=0;n=0;
// printf("dt=%f\n",dt);
 for(t=0;t<tmax;t+=dt)
  {
    A1=A01-m1*(C1-C01);
    dC1=k1*(A01-m1*(C1-C01))*C1;
    P1=-1000.0*H1*dC1/dt;
    A2=A02-m2*(C2-C02);
    dC2=k2*(A02-m2*(C2-C02))*C2;
    P2=-1000.0*H2*dC2/dt;
    Psum=P1+P2;
//    printf("%.1f, %f, %f, %f, %f, %f, %f, %f, %f, %f\n",t,dC1/dt,dC2/dt,A1,A2,C1,C2,Hf1,Hf2,Hf);
    C1+=dC1;
    C2+=dC2;
    if(t<tm[0]) continue;
    dQ=flow[n]-Psum;
    Q+=dQ*dQ;
//    printf("t=%f %f\n",t,tm[n]);
    n++;
  }
 return Q/n;
}

double my_f (const gsl_vector *v, void *params) {
  double x, y;
  double *p = (double *)params;
  double k1, k2, C01, C02, H1, H2, A01, A02;
  double r;

  k1 = gsl_vector_get(v, 0);
  k2 = gsl_vector_get(v, 1);
  C01 = gsl_vector_get(v, 2);
  C02 = gsl_vector_get(v, 3);
  H1 = gsl_vector_get(v, 4);
  H2 = gsl_vector_get(v, 5);
  A01 = p[0];
  A02 = p[1];

  r = test_kinetic(k1, A01, C01, H1, k2, A02, C02, H2);
//  printf("*** Q=%e  k1=%e  k2=%e  C01=%e  C02=%e  H1=%e  H2=%e\n", r, k1, k2, C01, C02, H1, H2);
  return r;
}

void my_df (const gsl_vector *v, void *params, gsl_vector *df) {
  double f0, f1, f2, dx;
  double *p = (double *)params;
  gsl_vector *vn;
  long i;

  vn = gsl_vector_alloc (Nx);
  for(i=0; i<Nx; i++) gsl_vector_set(vn, i, gsl_vector_get(v, i) );
  f0 = my_f(v, params);
//  printf("\n");
  for(i=0; i<Nx; i++) {
//    dx = 0.001 + fabs( gsl_vector_get(v, i) )*0.01;
    dx = fabs( gsl_vector_get(v, i) )*0.00001;
//    if(dx==0.0) dx=1e-8;
    gsl_vector_set(vn, i, gsl_vector_get(v, i)-dx/2 );
    f1 = my_f(vn, params);
    gsl_vector_set(vn, i, gsl_vector_get(v, i)+dx/2 );
    f2 = my_f(vn, params);
    gsl_vector_set(df, i, (f2-f1)/dx );
//    printf("df[%ld]=%e (%f +- %f) (%e, %e)\n", i, (f2-f1)/dx, gsl_vector_get(v, i), dx, f2, f1 );
    gsl_vector_set(vn, i, gsl_vector_get(v, i) );
  }
  gsl_vector_free (vn);
}

void my_fdf (const gsl_vector *x, void *params, double *f, gsl_vector *df) {
  *f = my_f(x, params);
  my_df(x, params, df);
}

void randomize() {
 struct timeval tv;
 struct timezone tz;

 gettimeofday(&tv,&tz);
 srandom(tv.tv_usec);
}

long rnd(long nmax) {return(rand()%nmax);}
double frnd() {return(1.0*rnd(1000000L)/1000000L);}

int minimizer_BFGS (double *rx, double *par) {
  size_t iter = 0;
  int status;
  long i;

  randomize();
  const gsl_multimin_fdfminimizer_type *T;
  gsl_multimin_fdfminimizer *s;

  gsl_vector *x;
  gsl_multimin_function_fdf my_func;

  my_func.n = Nx;
  my_func.f = my_f;
  my_func.df = my_df;
  my_func.fdf = my_fdf;
  my_func.params = par;

  x = gsl_vector_alloc (Nx);
  for(i=0;i<Nx;i++) gsl_vector_set (x, i, rx[i]*(1.25-0.5*frnd()) );
//  for(i=0;i<Nx;i++) gsl_vector_set (x, i, rx[i]);

  T = gsl_multimin_fdfminimizer_conjugate_fr;
  s = gsl_multimin_fdfminimizer_alloc (T, Nx);

  gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.0002, 1e-8);

//  printf ("*** %5d [", iter);
//  for(i=0;i<Nx;i++) printf (" %.9f", gsl_vector_get (s->x, i));
//  printf (" ] Q=%e\n", s->f);

  do {
      iter++;
      status = gsl_multimin_fdfminimizer_iterate (s);

      if (status) break;

      status = gsl_multimin_test_gradient (s->gradient, 1e-9);

      if (status == GSL_SUCCESS) printf ("Minimum found at:\n");

//      printf ("*** %5d [", iter);
//      for(i=0;i<Nx;i++) printf (" %.9f", gsl_vector_get (s->x, i));
//      printf (" ] Q=%e\n", s->f);

  } while (status == GSL_CONTINUE && iter < 10000);


  printf ("*** %5d [", iter);
  for(i=0;i<Nx;i++) printf (" %.9f", gsl_vector_get (s->x, i));
  printf (" ] Q=%e\n", s->f);

  gsl_multimin_fdfminimizer_free (s);
  gsl_vector_free (x);

  return 0;
}

int main(int argc, char *argv[])
{
  double k1avg,k2avg,C01avg,C02avg,H1avg,H2avg,Q,k1,k2,C01,C02,H1,H2;
  double Qmin,k1opt,k2opt,C01opt,C02opt,H1opt,H2opt;
  int j;
  long n;
  char *fname;

  double rx[Nx];
  double par[Npar] = { 0.0035185, 0.0035185 }; // A01, A02 start concentration

  for(j=1;j<argc;j++) {
    if(strcmp(argv[j],"-k1")==0) {k1avg=atof(argv[j+1]);j++;continue;}
    if(strcmp(argv[j],"-k2")==0) {k2avg=atof(argv[j+1]);j++;continue;}
    if(strcmp(argv[j],"-C01")==0) {C01avg=atof(argv[j+1]);j++;continue;}
    if(strcmp(argv[j],"-C02")==0) {C02avg=atof(argv[j+1]);j++;continue;}
    if(strcmp(argv[j],"-H1")==0) {H1avg=atof(argv[j+1]);j++;continue;}
    if(strcmp(argv[j],"-H2")==0) {H2avg=atof(argv[j+1]);j++;continue;}
    if(strcmp(argv[j],"-n")==0) {n=atol(argv[j+1]);j++;continue;}
    if(strcmp(argv[j],"-f")==0) {fname=strdup(argv[j+1]);j++;continue;}
  }
  read_hf(fname);

  rx[0]=k1avg;
  rx[1]=k2avg;
  rx[2]=C01avg;
  rx[3]=C02avg;
  rx[4]=H1avg;
  rx[5]=H2avg;
  minimizer_BFGS(rx, par);
}
