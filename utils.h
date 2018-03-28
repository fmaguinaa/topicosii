//#include "utils.h"

double *vector(long size) {
	double *v;
	v=(double *)calloc((size+1), sizeof(double));
	((long*)v)[0]=size;
	return v;
}

double **matrix(long size) {
	// Speicher für Zeilen und Spalten allozieren
	double** m= (double **) calloc(size+1, sizeof(double*));
	
    // Speicher für Zeilen allozieren
    for (int i=0; i<=size; i++) {
        m[i]= (double *) calloc(size+1, sizeof(double));
    }
    
    ((long*)m[0])[0]=size;
    return m;
}

double **matrix_unity(long size) {
  double** m=matrix(size);
  for (int i=1; i<=size; i++) {
      m[i][i]=1;
  }
  return m;
}

long m_getsize(double** m) {
  return ((long*)m[0])[0];    
}

long v_getsize(double* v) {
  return ((long*)v)[0];    
}


void m_free(double **m)
{
  long size=m_getsize(m);   
  for (int i=0; i<=size; i++) {
    free(m[i]);
  }
  free(m);
}

void v_free(double *v)
{
	free(v);
}



void m_print(double** a) {
  long size=m_getsize(a);   
  for (int i=1; i<=size; i++) {
    for (int j=1; j<=size; j++) {
//      cout<<a[i][j]<<"\t";
    }
//    cout<<endl;
  }
//  cout<<endl;
}

void v_print(double* a) {
  long size=v_getsize(a);   
//  cout<<"(";
  for (int i=1; i<=size; i++) {
//      cout<<a[i];
//      if (i!=size) cout<<",\t";
  }
//  cout<<")"<<endl<<endl;
}

double** m_mult(double** a, double** b) {
  long size=m_getsize(a);   
  double** res=matrix(size);
  
  for (int i=1; i<=size; i++) {
    for (int j=1; j<=size; j++) {
      res[i][j]=0;
      for (int k=1; k<=size; k++) {
        res[i][j] = res[i][j] + a[i][k]*b[k][j];
      }	
    }
  }
  
  return res;
}

void m_mult(double** a, double** b, double** c) {
  long size=m_getsize(b);   
  
  for (int i=1; i<=size; i++) {
    for (int j=1; j<=size; j++) {
      a[i][j]=0;
      for (int k=1; k<=size; k++) {
        a[i][j] = a[i][j] + b[i][k]*c[k][j];
      }	
    }
  }
}

void m_multself(double** b, double** c) {
  long size=m_getsize(b);   
  double** a=matrix(size);
  for (int i=1; i<=size; i++) {
    for (int j=1; j<=size; j++) {
      a[i][j]=0;
      for (int k=1; k<=size; k++) {
        a[i][j] = a[i][j] + b[i][k]*c[k][j];
      }	
    }
  }
//  m_copy(a,b);
//  m_free(a);
}

double** m_multconst(double** a, double c){
  long size=m_getsize(a);   
  double** res=matrix(size);
  
  for (int i=1; i<=size; i++) {
    for (int j=1; j<=size; j++) {
      res[i][j]=c*a[i][j];
    }
  }
  
  return res;
}

void m_multconst(double** a, double** b, double c){
  long size=m_getsize(b);   
  
  for (int i=1; i<=size; i++) {
    for (int j=1; j<=size; j++) {
      a[i][j]=c*b[i][j];
    }
  }
}

void m_multconstself(double** a, double c){
  long size=m_getsize(a);   
  
  for (int i=1; i<=size; i++) {
    for (int j=1; j<=size; j++) {
      a[i][j]=c*a[i][j];
    }
  }
}

double** m_add(double** a, double** b){
  long size=m_getsize(a);   
  double** res=matrix(size);
  
  for (int i=1; i<=size; i++) {
    for (int j=1; j<=size; j++) {
      res[i][j]=a[i][j]+b[i][j];
    }
  }
  
  return res;
}

void m_add(double** a, double** b, double** c){
  long size=m_getsize(b);   
  
  for (int i=1; i<=size; i++) {
    for (int j=1; j<=size; j++) {
      a[i][j]=b[i][j]+c[i][j];
    }
  }
}

void m_addself(double** a, double** c){
  long size=m_getsize(a);   
  
  for (int i=1; i<=size; i++) {
    for (int j=1; j<=size; j++) {
      a[i][j]=a[i][j]+c[i][j];
    }
  }
}

double* v_add(double* a, double* b){
  long size=v_getsize(a);   
  double* res=vector(size);
  
  for (int i=1; i<=size; i++) {
      res[i]=a[i]+b[i];
  }
  
  return res;
}

void v_add(double* a, double* b, double* c){
  long size=v_getsize(a);   
  
  for (int i=1; i<=size; i++) {
      a[i]=b[i]+c[i];
  }
}

void v_addself(double* a, double* c){  
  long size=v_getsize(a);   
  
  for (int i=1; i<=size; i++) {
      a[i]=a[i]+c[i];
  }
}


double** m_copy(double** a) {
  long size=m_getsize(a);   
  double** res=matrix(size);
  
  for (int i=1; i<=size; i++) {
    for (int j=1; j<=size; j++) {
      res[i][j]=a[i][j];
    }
  }
  
  return res;
}


double* v_copy(double* a) {
  long size=v_getsize(a);   
  double* res=vector(size);
  
  for (int i=1; i<=size; i++) {
      res[i]=a[i];
  }
  
  return res;
}


void m_copy(double** a, double** b) {
  long size=m_getsize(a);   
  
  for (int i=1; i<=size; i++) {
    for (int j=1; j<=size; j++) {
      a[i][j]=b[i][j];
    }
  }
}


void v_copy(double* a, double* b) {
  long size=v_getsize(a);   
  
  for (int i=1; i<=size; i++) {
      a[i]=b[i];
  }
}


double* vm_multr(double** a, double* b) {
  long size=v_getsize(b);   
  double* res=vector(size); 
  for (int i=1; i<=size; i++) {
    res[i]=0;  
    for (int j=1; j<=size; j++){
      res[i]=res[i]+a[i][j]*b[j];  
    }    
  }     
  return res;
} 

void vm_multr(double* a, double** b, double* c) {
  long size=m_getsize(b);   
  for (int i=1; i<=size; i++) {
    a[i]=0;  
    for (int j=1; j<=size; j++){
      a[i]=a[i]+b[i][j]*c[j];  
    }    
  }     
}

double* vm_multl(double* a, double** b) {
  long size=v_getsize(a);   
  double* res=vector(size); 
  for (int i=1; i<=size; i++) {
    res[i]=0;  
    for (int j=1; j<=size; j++){
      res[i]=res[i]+a[j]*b[j][i];  
    }    
  }     
  return res;
}

void vm_multl(double* a, double* b, double** c) {
  long size=v_getsize(b);   
  for (int i=1; i<=size; i++) {
    a[i]=0;  
    for (int j=1; j<=size; j++){
      a[i]=a[i]+b[j]*c[j][i];  
    }    
  }     
}

double* v_multconst(double* a, double c) {
  long size=v_getsize(a);   
  double* res=vector(size);
  
  for (int i=1; i<=size; i++) {
      res[i]=a[i]*c;
  }
  
  return res;
}

void v_multconst(double* a, double* b, double c) {
  long size=v_getsize(b);   
  
  for (int i=1; i<=size; i++) {
      a[i]=b[i]*c;
  }
}

void v_multconstself(double* a, double c) { 
  long size=v_getsize(a);   
  
  for (int i=1; i<=size; i++) {
      a[i]=a[i]*c;
  }
}

double v_skalar(double* a, double* b) {
  long size=v_getsize(b); 
  double s=0;  
  
  for (int i=1; i<=size; i++) {
      s=s+a[i]*b[i];
  }
  return s;
}


void v_crossself(double* a, double* b){
  long size=v_getsize(a);   
  double* res=vector(size);
  res[1]=a[2]*b[3]-a[3]*b[2]; 
  res[2]=a[3]*b[1]-a[1]*b[3]; 
  res[3]=a[1]*b[2]-a[2]*b[1]; 
  v_copy(a,res);
  v_free(res);
}

double* v_cross(double* a, double* b){
  long size=v_getsize(a);   
  double* res=vector(size);
  res[1]=a[2]*b[3]-a[3]*b[2]; 
  res[2]=a[3]*b[1]-a[1]*b[3]; 
  res[3]=a[1]*b[2]-a[2]*b[1]; 
  return res;
}

void v_cross(double* a, double* b, double* c){
  long size=v_getsize(a);   
  double* res=vector(size);
  a[1]=b[2]*c[3]-b[3]*c[2]; 
  a[2]=b[3]*c[1]-b[1]*c[3]; 
  a[3]=b[1]*c[2]-b[2]*c[1]; 
}

double v_abs(double* a) {
  long size=v_getsize(a);   
  double ab=0;
  for (int i=1; i<=size; i++)
    { ab+=a[i]*a[i]; }
  return sqrt(ab);
}
void v_normalize(double* a){
  long size=v_getsize(a);   
  double ab=v_abs(a);
  for (int i=1; i<=size; i++)
    { a[i]=a[i]/ab; } 
}


