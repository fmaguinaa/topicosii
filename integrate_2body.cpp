#include <cstdlib>
#include <iostream>
#include <cmath>

#include "utils.h"
#include "integrate_ode.h"

using namespace std;

/*
problema de dos cuerpos
      r' = v
      v' = -r/r^3
*/

#define dimension 4


void leapfrog_step(double* yout, double* y, double t0, double h,
	void (*f)(double*, double*, double)) {
  long size=v_getsize(yout);
  double* rhalf=vector(size);
  double* a=vector(size);
  double t=t0;
  v_copy(yout, y);

  // r_{1/2} = r_0 + h/2*v
  rhalf[1]=y[1]+h/2.0*y[3];     
  rhalf[2]=y[2]+h/2.0*y[4];   
  rhalf[3]=0;  
  rhalf[4]=0;  
  
  // v_1 = v_0 + h*F(r_{1/2})
  f(a, rhalf, t0+h/2);
  yout[3]=y[3]+h*a[3];
  yout[4]=y[4]+h*a[4];
  
  // r_1 = r_{1/2} + h/2 * v_1
  yout[1]=rhalf[1]+h/2.0*yout[3];
  yout[2]=rhalf[2]+h/2.0*yout[4];
  
  v_free(rhalf);
  v_free(a);
}


void f(double* yout, double* y, double x) {
  double r=sqrt(y[1]*y[1]+y[2]*y[2]);   
  yout[1]=y[3];
  yout[2]=y[4];
  yout[3]=-y[1]/(r*r*r);
  yout[4]=-y[2]/(r*r*r);
}

void angular_momentum(double* j, double* pos) {
     double* r=vector(3);
     double* p=vector(3);
     
     r[1]=pos[1]; r[2]=pos[2]; r[3]=0;
     p[1]=pos[3]; p[2]=pos[4]; p[3]=0;
     
     v_cross(j, r, p);
     v_free(r);
     v_free(p);
}

void runge_lenz(double* e, double* pos) {
     double* r=vector(3);
     double* p=vector(3);
     double* j=vector(3);
     
     r[1]=pos[1]; r[2]=pos[2]; r[3]=0;
     p[1]=pos[3]; p[2]=pos[4]; p[3]=0;
     
     angular_momentum(j, pos);
     v_cross(e, p, j);
     v_normalize(r);
     v_multconstself(r, -1.0);
     v_addself(e,r);
     
     v_free(r);
     v_free(p);
     v_free(j);
}

double energy(double* pos) {
  return 0.5*(pos[3]*pos[3]+pos[4]*pos[4])-1/sqrt(pos[1]*pos[1]+pos[2]*pos[2]);
}

int main(int argc, char *argv[])
{
    double* y0=vector(dimension);
    double* yout2=vector(dimension);
    double* yout4=vector(dimension);
    double* youte=vector(dimension);
    double* youtl=vector(dimension);
    double h=0.01;
    double tmax=100;
    double E0, E;
    int count=0;
    int cmax=20;
    y0[1]=1.0;
    y0[2]=0;
    y0[3]=0;
    y0[4]=1.05;
    v_copy(yout2, y0);
    v_copy(yout4, y0);
    v_copy(youte, y0);
    v_copy(youtl, y0);
    E0=energy(y0);
    cout.precision(15);
    for(double t=0; t<=tmax; t+=h) {
      count++;
      if (count>cmax) cout<<t<<"\t";

      v_copy(y0, youte);
      euler_step(youte, y0, t, h, f);
      E=energy(youte);
      if (count>cmax) cout <<youte[1]<<"\t"<<youte[2]<<"\t"<<abs((E-E0)/E0)<<"\t";

      v_copy(y0, yout2);
      rk2_step(yout2, y0, t, h, f);
      E=energy(yout2);
      if (count>cmax) cout <<yout2[1]<<"\t"<<yout2[2]<<"\t"<<abs((E-E0)/E0)<<"\t";

      v_copy(y0, yout4);
      rk4_step(yout4, y0, t, h, f);
      E=energy(yout4);
      if (count>cmax) cout <<yout4[1]<<"\t"<<yout4[2]<<"\t"<<abs((E-E0)/E0)<<"\t";

      v_copy(y0, youtl);
      leapfrog_step(youtl, y0, t, h, f);
      E=energy(youtl);
      if (count>cmax) cout <<youtl[1]<<"\t"<<youtl[2]<<"\t"<<abs((E-E0)/E0)<<"\t";
       
      if (count>cmax) {cout<<endl; count=0;}
    }                
    v_free(y0);
    v_free(yout2);
    v_free(yout4);
    v_free(youte);
    v_free(youtl);
    return EXIT_SUCCESS;
}
