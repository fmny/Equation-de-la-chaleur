//même programme que Q5impli2.cpp sauf que la condition initiale n'est plus 
//u(x,0)=0 mais u(x,0)=sin(pi*x/L)+0.25*sin(10*pi*x/L)
//on compare la solution numerique a la solution analytique
//calulée pour un mur d'épaisseur finie

#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

double F(double x,double L){ return sin(M_PI*x/L)+0.25*sin(10*M_PI*x/L);}

double us=0;
const double L = 1.;
double kappa=1.;
const int M= 50; 
const int itermax=20000;  
double h= L/M, h2=h*h;
double deltaT=h*h/(2*kappa);  
int compteur=0;




//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//Représentation de la solution analytique
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  

//onde=A[k]*exp(-(k*PI/L)²)phik[x]
//sommedesondes=sum(A[k])   k=1,..,20



void solanalytique ( int M1, int itermax1 )

{
  ofstream graph2("Q52");
  double *A=new double [20];
  double *onde=new double[21];
  double *sommedesondes=new double[M1+1];
  A[1]=1.; A[10]=1./4;  

  for ( int t=0;t<=itermax1;t++ )
    {
      for ( int i=0;i<=M1;i++ )
	{
	  for ( int k=1;k<=20;k++ )
	    {
	      onde[k]=onde[k-1]+ 
		A[k]*exp((-k*M_PI*(1./L))*(k*M_PI*(1./L))*t*deltaT )	
		*sin(k*M_PI*i*h*(1./L)); 
	    }
	 
	  sommedesondes[i]=onde[20];

		  
	  if ( t%100 == 1 || (t<100 && t%10==1 )||(t<10))
	    {
	      graph2<<i*h<<"   "<<(1-i*h/L)*us+sommedesondes[i]<<endl; 
	    }
	}
      
    }
}



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//Début de l'algorithme de decomposition LU
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

void  LU ( double *a, double *b, double *c, double *X ,double *f ) 

{
  double *bstar=new double[M+1];
  double *cstar=new double[M+1];
  double *Y=new double[M+1];
   
 
  bstar[1]=b[1];
  cstar[1]=c[1]/b[1];

  for ( int k=2; k<=M; k++ )
    {      bstar[k]=b[k]-a[k]*cstar[k-1];
      cstar[k]=c[k]/bstar[k];    }

  Y[1]=f[1]/bstar[1];
  for ( int k=2; k<=M; k++ )
    {    Y[k]=(f[k]-a[k]*Y[k-1])/bstar[k]; }

  X[M]=Y[M];
  for ( int k=M-1; k>=1; k-- )
    {    X[k]=Y[k]-cstar[k]*X[k+1] ;    }
}

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //fin de l'algorithme de  décomposition LU
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




int main()
{
  ofstream graph("Q51");
  double *a=new double[M+2];
  double *b=new double[M+2];
  double *c=new double[M+2];
  double *f=new double[M+2];
  double *Y=new double[M+2];
  double *X=new double[M+2];
  
  //condition aux limites u(0,t)=us ; u(x,0)=0 ; 
  
  //solution initiale juste
 for ( int i=1;i<=M+1; i++ )
   {   f[i]=F((i-1)*h,L) ;    X[i]=us;  }
 
 //ok
  b[1]=1.;   
  for ( int i=2;i<=M;i++ )
    { b[i]=2*kappa*deltaT/h2+1; }
  b[M+1]=1.;      
  
  //ok
  c[1]=0.;
  for ( int i=2;i<M+1;i++ )
    { c[i]=-kappa*deltaT/h2; }
  
  //ok
  a[M+1]=0.;
  for ( int i=2;i<M+1;i++ )
    { a[i]=-kappa*deltaT/h2 ; }
  
  //
  for ( int i=0;i<=M;i++)
    {  graph<<i*h<<"   "<<F(i*h,L)<<endl; }

  LU(a,b,c,X,f);

  for ( int t=0; t<=itermax; t++ )
    {    
      f[1]=us;   f[M+1]=0;  X[M+1]=0;    X[1]=us;
            
      for ( int i=2;i<=M; i++ )
	{  f[i]=X[i] ; }
      
      for ( int i=1;i<=M+1; i++ )
	{ Y[i]=X[i]; }
      //On stocke les valeurs de la fonction a l'instant precedent
      //pour permettre le calcul de l'erreur      
 
      LU (a,b,c,X,f);
    
      for ( int i=0;i<=M; i++ )
	{
	  if ( t%100 == 1 || (t<100 && t%10==1 )||(t<10))
	    graph<<i*h<<"   "<<X[i+1]<<endl;
	}
      
      
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      // Calcul de l'erreur pour tester la vitesse
      //de convergence
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      double err=0;
      for (int i=1; i<= M+1; i=i+1)     
	{ err += abs(X[i]-Y[i]);}
      //         cout << iter << " erreur =  " << err << endl;
      if(err < 1e-10) break; 
      
      for ( int i=1;i<=M+1;i++ )
	{Y[i]=X[i];}
      compteur=compteur+1;
      
    }
  cout<<"convergence du schéma implicite en "<<compteur<<" itérations."<<endl;

  

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //solution analytique pour un domaine fini
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  solanalytique(M,itermax);
  
 
}
	
