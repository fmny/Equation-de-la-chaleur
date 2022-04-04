//solution numérique qui semble coller avce la solution analytique


#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;
typedef double R;


double kappa=1.,L=1.,us=0.;

int M=50,itermax=10000;


double h=L/M;
double deltaT=h*h/(2*kappa);

R max (R a,R b)
{
  if (a<b) return b;
  else return a;
}


 //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //Représentation de la solution analytique
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
  //onde=A[k]*exp(-(k*PI/L)²)phik[x]
  //sommedesondes=sum(A[k])   k=1,..,20
  
void solanalytique (int M1,int itermax1)  
{
  double *onde=new double[21];
  double *sommedesondes=new double[M1+1];
  double *A=new double [21];
 
  for ( int i=0;i<20;i++ )
    A[i]=0;//-2*us/(i*M_PI);

 A[1]=1.; A[10]=1./4;

  ofstream graph2("Q42");
  

  for ( int t=0;t<=itermax1;t++ )
    {
      for ( int i=0;i<=M1;i++ )
	{

	   for ( int k=1;k<=20;k++ )
	     {
	       onde[k]=onde[k-1]+ 
		 (A[k])*
		 exp((-k*M_PI*(1./L))*(k*M_PI*(1./L))*t*deltaT )	
		 *sin(k*M_PI*i*h*(1./L)); 
	     }
	  
	  sommedesondes[i]=onde[20];
	  
	  //  cout<<sommedesondes[i]<<endl;
	  
	  if ( t%200 == 0 || (t<100 && t%5==0 ))
	    {
	      graph2<<i*h<<"   "<<(1-i*h/L)*us+sommedesondes[i]<<endl; 
	    }
	}
      
      
    }
}
  
  


int main()
{
  double *u = new double[M+1];
  double *v = new double[M+1];
  //double erreur;
  
  //solution numerique
  ofstream graph("Q41");  
  ofstream graph0("amortissement");  
  double amort,amortmax;
  int compteur;
  compteur=0;
  amortmax=0;

  for ( int i=0 ; i<=M ; i++ )
	{
	  u[i]=sin(M_PI*i*h/L)+(1./4)*sin(10*M_PI*i*h/L) ;
	  graph<<i*h<<"   "<<u[i]<<endl;;    
	}

  for (int t=0;t<=itermax;t++)
    {      
      v[0]=us; v[M]=0;
      for ( int i=1 ; i<M ; i++ )
	{ v[i]=(u[i]+(kappa*deltaT/(h*h))*(u[i+1]-2*u[i]+u[i-1])); }
      

      //amort max=1.22041 a peu pres
      amort=0;
      for ( int i=0;i<M;i++)
	{ amort= max(amort,v[i]); }

      amortmax=max(amortmax,amort);

      //      cout<<"amormax="<<amortmax<<endl;

      if ( amort > amortmax/2 )
	{compteur=compteur+1;}


     if ( t< 400  )
	{  graph0<<t<<"    "<<amort<<endl;}

  if ( t%200 == 0 || (t<100 && t%5==0 ))	
   {
     for ( int i=0;i<=M;i++ )
       {  graph<<i*h<<"   "<<v[i]<<endl; }
   }	

      for ( int i=0;i<=M;i++ )
	u[i]=v[i];
    }

  cout<<"nombre d'itérations pour que l'amplitude de l'onde soit plus que divisée par deux="<<compteur<<endl;
  
 
  solanalytique (M,itermax);

}
