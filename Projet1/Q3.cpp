//Réponse à la question 3),mettre dans un fichier graphique
//Q31 et Q32 puis dans un autre Q31 et Q33


#include <iostream>
#include <cmath>
#include <fstream>

using namespace std;


double kappa=1.,L=1.,us=1.,u0=0.;

int M=50,itermax=50000;
double h=L/M;
double deltaT=h*h/(2*kappa);//deltaT=h*h/(2*kappa)
double theta0=0.;//a verifier que theta0=0;



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//Représentation de la solution analytique pour le cas du domaine fini
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  //onde=A[k]*exp(-(k*PI/L)²)phik[x]
  //sommedesondes=sum(A[k])   k=1,..,20

void solanalytique ( int M1, int itermax1 )
  
{
  ofstream graph2("Q32");
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

	  
	  if (( t<100)&&(t%5==0) ||   (t%200 == 0)  )
	    {
	      //graph2<<t*deltaT<<"    "<<i*h
	      //    <<"  "<<(1-i*h/L)*us+sommedesondes[i]<<endl; 	   
	    
	      graph2<<i*h<<"  "<<(1-i*h/L)*us+sommedesondes[i]<<endl; 	   
	    


	    }
	}
    }
}


  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  //solution analytique pour le cas d'un mur de longueur infini
  //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  
void solanalytique2 ( int M2 ,int itermax2 )
{

  ofstream graph3("Q33");
  double *H=new double[M2+1];

  for ( int i=0;i<=M2;i++ )
    {  H[i]=0; }

  for ( int i=0;i<=M2;i++ )
    {  graph3<<i*h<<"    "<<H[i]<<endl; }
   
  for ( int t=1;t<=itermax2;t++ )
    {
      for ( int i=0;i<=M2;i++ )
	{  H[i]=(1-erf(i*h*(1./(2*sqrt(kappa*t*deltaT)))))*us; }
      
      for ( int i=0;i<=M2;i++ )
	{   	  
	  if ( ( t<100)&&(t%5==0) ||   (t%200 == 0)  )
	  { graph3<<i*h<<"   "<<H[i]<<endl; }
	}
      
    }

}



//Définition de la fonction erf
/*
double erf ( double L)
{
  int P=100;
  double h=L/P;
  double *u=new double[P];
  double v;


  
  if ( L<=0 ) return 0;
  
  else
  u[0]=h*exp(-h*h);
  
  for ( int i=0; i<=P-2 ; i++ )
  {	u[i+1]=u[i]+h*exp((-(i+1)*h)*((i+1)*h)); }
  
  v=2./(sqrt(M_PI))*u[P-1]; 
  
  return v;
  }
*/



int main()
{
  
  double *u = new double[M+1];
  double *v = new double[M+1];
  int compteur=0;
  ofstream graph("Q31");  
  ofstream err1("er1");  
  ofstream err2("er2");  
  double erreur1,erreur2; 

  double *onde1=new double[21];
  double *sommedesondes1=new double[M+1];
  

  for ( int i=0 ; i<M ; i++ )
    { 
      graph<<i*h<<"   "<<u0<<endl;
      u[i]=0;	
    }
  u[0]=us;
  
  
  for (int t=1;t<=itermax;t++)
    {      
      for ( int i=1 ; i<M ; i++ )
	{  v[i]=u[i]+(kappa*deltaT/(h*h))*(u[i+1]-2*u[i]+u[i-1]); }	
      
      v[M]=theta0; //u(l,t)=theta0 pout tout t
      v[0]=us;  //u0(t)=us,u(l,t)=0 pour tout t
      
      /*   
      for ( int i=0;i<=M;i++ )
	{  erreur1+=(( (1-erf(i*h*(1./(2*sqrt(kappa*t*deltaT)))))*us )-v[i] )
*(( (1-erf(i*h*(1./(2*sqrt(kappa*t*deltaT)))))*us )-v[i] ); }

      for ( int i=0;i<=M;i++ )
	{
	  for ( int k=1;k<=20;k++ )
	    {
	      onde1[k]=onde1[k-1]+ 
		( -2*us*(1./(k*M_PI)))*
		exp((-k*M_PI*(1./L))*(k*M_PI*(1./L))*t*deltaT )	
		*sin(k*M_PI*i*h*(1./L)); 
	    }

	  sommedesondes1[i]=onde1[20];
	}
            
     
      for ( int i=0;i<=M;i++ )
	{  erreur2+=((1-i*h/L)*us+sommedesondes1[i]-v[i] )*
	     ((1-i*h/L)*us+sommedesondes1[i]-v[i] )  ; }
      */

      for ( int i=0 ; i<M ; i++ )
	{  
	  //	  if ( t% 200 == 0 )	
	    //  graph<<t*deltaT<<"       "<<i*h<<"   "<<v[i]<<endl;
	  
	  if ( ( t<100)&&(t%5==0)|| (t%200==0 ))
	    {  graph<<i*h<<"    "<<v[i]<<endl; }
	

	  /*       if (t<100)
	  {
	    err1<<t*deltaT<<"  "<<sqrt(erreur1)<<endl; 	   
	    err2<<t*deltaT<<"   "<<sqrt(erreur2)<<endl;
	  }
	  */

	}


   

      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      // Calcul de l'erreur pour tester la vitesse
      // de convergence
      //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      double err=0;
      for (int i=0; i<= M; i=i+1)     
	{ err += abs(u[i]-v[i]);}
      
      if(err < 1e-10) break; 
      
      for ( int i=1;i<=M-1;i++ )
	{u[i]=v[i];}
      compteur=compteur+1;
    }
  
  cout<<"convergence du schéma explicite en "<<compteur<<" itérations."<<endl;

  solanalytique(M,compteur);  
  solanalytique2(M,compteur);  


  //calcul de l'erreur entre la solution numerique
  //et les 2 solutions exactes

  /* 
  for ( int t=0;t<=itermax;t++ )
    {
      for ( int i=0;i<=M;i++ )
	{
	  
	  if ( ( t<100 )&&(t%5==0) ||   (t%200 == 0)  )
	    {
	      
	      //graph2<<t*deltaT<<"    "<<i*h
	      // <<"  "<<(1-i*h/L)*us+sommedesondes[i]<<endl; 	   
	    
	      graph2<<i*h<<"  "<< ( (1-i*h/L)*us+sommedesondes[i] )-v[i]-<<endl; 	   
      
	      
	    }
	}
    }
  
  */
  
  
}
