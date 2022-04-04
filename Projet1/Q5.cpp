//reponse a la question Q5.3)

#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;

double F(double x){ return 0;}

double us=1.;
const double L = 1.;
double kappa=1.;
const int M= 50; 
const int itermax=1000; 
double h= L/M, h2=h*h;
double deltaT=h*h/(2*kappa);  
int compteur=0;
const int gv=100000; 



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//Représentation de la solution analytique
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  

  //onde=A[k]*exp(-(k*PI/L)²)phik[x]
  //sommedesondes=sum(A[k])   k=1,..,20



void solanalytique ( int M1, int itermax1 )

{
  ofstream graph2("Q52");//domaine fini
  double *onde=new double[21];
  double *sommedesondes=new double[M1+1];
  
  for ( int t=0;t<=itermax1;t++ )
    {
      for ( int i=0;i<=M1;i++ )
	{
	  for ( int k=1;k<=20;k++ )
	    {
	      onde[k]=onde[k-1]+ 
		( -2*us*(1./(k*M_PI)))*
		exp((-k*M_PI*(1./L))*(k*M_PI*(1./L))*t*deltaT )	
		*sin(k*M_PI*i*h*(1./L)); 
	    }
	 
	  sommedesondes[i]=onde[20];
		  
	  if ( t%20 == 0)
	    {
	      graph2<<i*h<<"   "<<(1-i*h/L)*us+sommedesondes[i]<<endl; 
	    }
	}
      
    }
}


  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //solution analytique pour le cas d'un mur de longueur infini
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
void solanalytique2 ( int M2 ,int itermax2 )
{
  ofstream graph3("Q53");//domaine infini
  double *H=new double[M2+1];
  
  for ( int i=0;i<=M2;i++ )
    {  graph3<<i*h<<"    "<<0<<endl; }
  
  for ( int t=1;t<=itermax2;t++ )
    {
      for ( int i=0;i<=M2;i++ )
	{  H[i]=(1-erf(i*h/(2*sqrt(kappa*t*deltaT))))*us; }
      
      for ( int i=0;i<=M2;i++ )
	{      if ( t%20 == 0)
	  { 
	    graph3<<i*h<<"   "<<H[i]<<endl; 
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
  
  double *a=new double[M+1];
  double *b=new double[M+1];
  double *c=new double[M+1];
  double *f=new double[M+1];
  double *Y=new double[M+1];
  double *X=new double[M+2];
  ofstream graph("Q51");
  //condition aux limites u(0,t)=us ; u(x,0)=0 ; 
  
 for ( int i=1;i<=M; i++ )
   {   f[i]=0 ;   }//u(x,0)=0
 
 for ( int i=0;i<=M;i++ )
   { X[i]=us; }



  b[1]=1.;   
  for ( int i=2;i<M;i++ )
    { b[i]=2*kappa*deltaT/h2+1; }
  b[M]=1.;      
  
  c[1]=0.;
  for ( int i=2;i<M;i++ )
    { c[i]=-kappa*deltaT/h2; }
  
  a[M]=0.;
  for ( int i=2;i<M;i++ )
    { a[i]=-kappa*deltaT/h2 ; }
  
  for ( int i=0;i<=M;i++)
    {  graph<<i*h<<"   "<<F(i*h)<<endl; }

  LU(a,b,c,X,f);

  for ( int iter=0; iter<gv; iter++ )
    {    
      f[1]=us;   f[M]=0;  X[M]=0;    X[1]=us;
      //On stocke les valeurs de la fonction a l'instant precedent
      //pour permettre le calcul de l'erreur
      
      for ( int i=2;i<M; i++ )
	{  f[i]=X[i] ; }
      
      for ( int i=1;i<=M; i++ )
	{ Y[i]=X[i]; }
       
      LU (a,b,c,X,f);

      //      graph<<0<<"   "<<us<<endl;  
   
      for ( int i=0;i<M; i++ )
	{ if (iter%20 == 0  )
	  graph<<i*h<<"   "<<X[i+1]<<endl;
	}
     
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      // Calcul de l'erreur pour tester la vitesse
      //de convergence
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      double err=0;
      for (int i=1; i<= M; i=i+1)     
	{ err += abs(X[i]-Y[i]);}
      // cout << iter << " err " << err << endl;
      if(err < 1e-10) break; 
      
      for ( int i=1;i<M+1;i++ )
	{Y[i]=X[i];}
      compteur=compteur+1;
      
    }
  cout<<"convergence du schéma implicite en "<<compteur<<" itérations."<<endl;
  //cout<<"erreur = "<<err<<endl;
  
  // for (int i=0;i<=M;i++ )
  //{ cout

  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //solution analytique pour un domaine fini et infini
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  solanalytique(M,itermax);
  solanalytique2(M,itermax);
  
}
	
