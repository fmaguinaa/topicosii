//#include"integrate_ode.h"

void euler(double* yout, double* y0, double a, double b, double h,
	void (*f)(double*, double*, double)) {
  long N=(long)abs((b-a)/h);
  long size=v_getsize(yout);
  double* fo=vector(size);
  double x=a;
  v_copy(yout, y0);
  for (int i=1; i<=N; i++) {
     //cout<<x<<",\t"<<yout[1]<<endl; 
     
     // yout = yout + h*f(yout, x)
     f(fo, yout, x); 
     v_multconstself(fo, h);
     v_addself(yout, fo);
     
     // x=x+h
     x=x+h; 
  }
  //cout<<x<<",\t"<<yout[1]<<endl; 
  v_free(fo);
}

void euler_step(double* yout, double* y0, double a, double h,
	void (*f)(double*, double*, double)) {
  long size=v_getsize(yout);
  double* fo=vector(size);
  double x=a;
  v_copy(yout, y0);

     // yout = yout + h*f(yout, x)
     f(fo, yout, x); 
     v_multconstself(fo, h);
     v_addself(yout, fo);

  v_free(fo);
}

void rk2(double* yout, double* y0, double a, double b, double h,
	void (*f)(double*, double*, double)) {
  long N=(long)abs((b-a)/h);
  long size=v_getsize(yout);
  double* fo=vector(size);
  double* k1=vector(size);
  double* k2=vector(size);
  double* yd=vector(size);
  double x=a, xd;
  v_copy(yout, y0);
  for (int i=1; i<=N; i++) {
     //cout<<x<<",\t"<<yout[1]<<endl; 
     
     // k1 = f(x,y)
     f(k1, yout, x); 
     
     // yd = y+ (h/2)*k1;    xd=x+(h/2)
     v_multconstself(k1, h/2.0);  
     v_add(yd, yout, k1);
     xd=x+h/2.0;
     
     // k2 = f(xd, yd)
     f(k2,yd,xd);
     
     // yout=yout + h * k2
     v_multconstself(k2, h);  
     v_addself(yout, k2);
     x=x+h; 
  }
  //cout<<x<<",\t"<<yout[1]<<endl; 
  v_free(yd);
  v_free(fo);
  v_free(k1);
  v_free(k2);
}

void rk2_step(double* yout, double* y0, double a, double h,
	void (*f)(double*, double*, double)) {
  long size=v_getsize(yout);
  double* fo=vector(size);
  double* k1=vector(size);
  double* k2=vector(size);
  double* yd=vector(size);
  double x=a, xd;
  v_copy(yout, y0);

     // k1 = f(x,y)
     f(k1, yout, x); 
     
     // yd = y+ (h/2)*k1;    xd=x+(h/2)
     v_multconstself(k1, h/2.0);  
     v_add(yd, yout, k1);
     xd=x+h/2.0;
     
     // k2 = f(xd, yd)
     f(k2,yd,xd);
     
     // yout=yout + h * k2
     v_multconstself(k2, h); 
     v_addself(yout, k2);

  v_free(yd);
  v_free(fo);
  v_free(k1);
  v_free(k2);
}

void rk4(double* yout, double* y0, double a, double b, double h,
	void (*f)(double*, double*, double)) {
  long N=(long)abs((b-a)/h);
  long size=v_getsize(yout);
  double* fo=vector(size);
  double* k1=vector(size);
  double* k2=vector(size);
  double* k3=vector(size);
  double* k=vector(size);
  double* k4=vector(size);
  double* yd=vector(size);
  double x=a, xd;
  v_copy(yout, y0);
  for (int i=1; i<=N; i++) {
     //cout<<x<<",\t"<<yout[1]<<endl; 
     
     // k1 = f(x,y)
     f(k1, yout, x); 
     v_copy(k,k1); // k=k1 -> Zwischenspeichern, da für Rechnung k1..k3 verändert werden muss
     
     // k2 = f(xd, yd)
     //    wobei: yd = y+ (h/2)*k1;    xd=x+(h/2)
     v_multconstself(k, h/2.0);  
     v_add(yd, yout, k);
     xd=x+h/2.0;
     f(k2,yd,xd);
     v_copy(k,k2); // k=k2
     
     // k3 = f(xd, yd)
     //    wobei: yd = y+ (h/2)*k2;    xd=x+(h/2)
     v_multconstself(k, h/2.0);  
     v_add(yd, yout, k);
     xd=x+h/2.0;
     f(k3,yd,xd);
     v_copy(k,k3); // k=k1
     
     // k4 = f(xd, yd)
     //    wobei: yd = y+ h*k3;    xd=x+h
     v_multconstself(k, h);  
     v_add(yd, yout, k);
     xd=x+h;
     f(k4,yd,xd);
    
     // yout=yout + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
     v_add(k, k1, k4);
     v_multconstself(k2, 2);
     v_multconstself(k3, 2);
     v_addself(k, k2);
     v_addself(k, k3);
     v_multconstself(k, h/6.0);  
     v_addself(yout, k);
     x=x+h; 
  }
  //cout<<x<<",\t"<<yout[1]<<endl; 
  v_free(yd);
  v_free(fo);
  v_free(k1);
  v_free(k2);
  v_free(k3);
  v_free(k4);
  v_free(k);
}


void rk4_step(double* yout, double* y0, double a, double h,
	void (*f)(double*, double*, double)) {
  long size=v_getsize(yout);
  double* fo=vector(size);
  double* k1=vector(size);
  double* k2=vector(size);
  double* k3=vector(size);
  double* k=vector(size);
  double* k4=vector(size);
  double* yd=vector(size);
  double x=a, xd;
  v_copy(yout, y0);

     // k1 = f(x,y)
     f(k1, yout, x); 
     v_copy(k,k1); // k=k1 -> Zwischenspeichern, da für Rechnung k1..k3 verändert werden muss
     
     // k2 = f(xd, yd)
     //    wobei: yd = y+ (h/2)*k1;    xd=x+(h/2)
     v_multconstself(k, h/2.0);  
     v_add(yd, yout, k);
     xd=x+h/2.0;
     f(k2,yd,xd);
     v_copy(k,k2); // k=k2
     
     // k3 = f(xd, yd)
     //    wobei: yd = y+ (h/2)*k2;    xd=x+(h/2)
     v_multconstself(k, h/2.0);  
     v_add(yd, yout, k);
     xd=x+h/2.0;
     f(k3,yd,xd);
     v_copy(k,k3); // k=k1
     
     // k4 = f(xd, yd)
     //    wobei: yd = y+ h*k3;    xd=x+h
     v_multconstself(k, h);  
     v_add(yd, yout, k);
     xd=x+h;
     f(k4,yd,xd);
    
     // yout=yout + (h/6) * (k1 + 2*k2 + 2*k3 + k4)
     v_add(k, k1, k4);
     v_multconstself(k2, 2);
     v_multconstself(k3, 2);
     v_addself(k, k2);
     v_addself(k, k3);
     v_multconstself(k, h/6.0);  
     v_addself(yout, k);

  v_free(yd);
  v_free(fo);
  v_free(k1);
  v_free(k2);
  v_free(k3);
  v_free(k4);
  v_free(k);
}
