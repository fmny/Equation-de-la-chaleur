//Réponse à la question Q6


#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cstring>
#include "sfem.hpp"
#include "RNM.hpp"
#include "GC.hpp"
#include "gnuplot-iso.hpp"


double L=2., H=2., l=1., h=0.4, k1=1., k2=10., us=100., un=50.;
double xr=(-l/2); 
double yr=(H+h)/4;
int nr=6;
double q=1.;


//fonction pour sauvegarder la solution 
// et la tracer sous gnuplot (commande splot)

void savesplot(const Mesh & Th,char * fileplotname,KN_<R> &x)
{
  ofstream plot(fileplotname);
  for(int it=0;it<Th.nt;it++)
    plot << (R2) Th[it][0] << " " << x[Th(it,0)] <<  endl 
     << (R2) Th[it][1] << " " << x[Th(it,1)] << endl 
     << (R2) Th[it][2] << " " << x[Th(it,2)] << endl 
     << (R2) Th[it][0] << " " << x[Th(it,0)] << endl << endl << endl; 
}

// une classe matrice diagonale pour le preconditionneaur diagonal
template <class R>
class MatDiag: VirtualMatrice<R> { public:
  const KN_<R> d; // Adr du Vecteur diagonal
  typedef typename  VirtualMatrice<R>::plusAx plusAx;
  
  MatDiag(KN_<R> dd) :d(dd) {}
  
  void addMatMul(const  KN_<R>  & x, KN_<R> & Ax) const { 
    assert(x.N()==Ax.N() && x.N() == d.N() );
    Ax+=DotStar_KN_<R>(x,d);   // <=>  for(int i=0;i<x.N;i++) Ax[i] += x[i]*d[i];
  }
  
  plusAx operator*(const KN<R> &  x) const {
    return plusAx(this,x);}
};
 

class Mat_lap: public VirtualMatrice<R> { public:
  const Mesh & Th;
  typedef VirtualMatrice<R>::plusAx plusAx;
  
  Mat_lap(const Mesh & T) : Th(T) {}; 
  
  void addMatMul(const  KN_<R>  & x, KN_<R> & Ax) const{
  for (int k=0;k<Th.nt;k++)
    {
      const Triangle & K(Th[k]);
      double conductivite = (K.lab == 0)? k2 : k1;
      for (int il=0;il<=2;il++)
	{
	  int ig0 = Th(k,il);
	  int ig1 = Th(k,((il+1) %3));
	  int ig2 = Th(k,((il+2) %3));
	  
	      if( Th(ig0).lab != 1 && Th(ig0).lab != 3) 
		{		  
		  Ax[ig0] += conductivite * K.area * (x[ig0] * (K.H(il) , K.H(il)) + x[ig1] * (K.H(il) , K.H((il+1) %3)) + x[ig2] * (K.H(il) , K.H((il+2) %3)) );
		}
	}
    }
  
  }; 
  plusAx operator*(const KN<R> &  x) const {return plusAx(this,x);}
};



int recherchelabel ( R2 &P , Mesh &T )
{   
  int labelpoint;  
  double lambda0,lambda1,lambda2;
  Mesh &Th(T);
  int k=0;

  while ( k <= Th.nt ) 
    {
      if ( k == Th.nt ) { return 10000; }	 
      
      lambda0=0.; lambda1=0.; lambda2=0.;    
      
      Triangle & K(Th[k]);
      
      R2 i0(K[0]),i1(K[1]),i2(K[2]); //coordonnées des 3 sommets 
      
      double a,b,c,d,e,f,det;
      a=0.; b=0.; c=0.; d=0.; e=0; f=0; det=0.;
      
      a=i0.x-i2.x; b=i1.x-i2.x; c=i0.y-i2.y; d=i1.y-i2.y; e=i2.x; f=i2.y;
      
      det=a*d-b*c;	  
      
      if (  a != 0   )
	{
	      lambda1=(a*P.y-(P.x-e)*c-f*a)*(1./det);
	      lambda0=(P.x-lambda1*b-e)*(1./a);
	      lambda2=1-lambda0-lambda1;
	}
      
      else
	{
	  lambda1=(P.x-e)*(1./b);
	  lambda0=(P.y-lambda1*d-f)*(1./c);
	  lambda2=1-lambda0-lambda1;	 
	}
      
      
      //la condition de sortie fonctionne
      if   (  (  lambda0 >= 0  )&&(  lambda1 >=0 )&&( lambda2 >= 0 )
	      &&(  lambda0 <= 1  )&&(  lambda1 <=1 )&&( lambda2 <= 1 ) )  
	{ labelpoint=K.lab; return labelpoint; }
      
      else { k=k+1; }
      
    }
}

R U ( int l, R2 &P , Mesh &T)
{
  Mesh &Th(T);
 
  for ( int k=1;k<=nr;k++ )
    {
      if ( l == k )
	{ 
	  if ( recherchelabel(P,Th) == 14-k )
	    { return q; }
	  else { return 0.;}
	}
       
  }
  
  
} 



int main(int argc , char** argv )
{
  Mesh Th(argc>1? argv[1] : "four.msh");
  int N=Th.nv;
  
  // construction de b
  KNM<R> B(N,nr+1);
    B=0;

  // Assemblage de la matrice et du second membre 
  for (int k=0;k<Th.nt;k++)
    {
      Triangle & K(Th[k]);//Th[k]= triangle k
       
      R2 sA(K[0]),sB(K[1]),sC(K[2]);  //coord des sommets A,B,C du triangle k
      
      R2 gw[3]={K.H(0),K.H(1),K.H(2)};       //K.H(0)=gradient 
      
      R2 G_K((sA+sB+sC)/3.);  //G_K barycentre du triangle k
      // R intfgkwi= f(G_K)*K.area/3; //  $ \int_K f w_i = f(G_K) |K|/3 $
      //intfgkwi=valeur de f evaluée en le
      // barycentre du triangle k *aire du triangle/3 
      // int aux=0;

      for ( int il=0;il<=2;il++)
	{
	  int ig = Th(k,il);//Th(k,il)=sommet de numero local il du triangle k	 
	  Vertex pt(Th(ig));//nomme pt le sommet ig
	  
	  if (  ( pt.lab ==1) ||( pt.lab ==3 )  )
	    { 
	      for ( int k=0;k<=nr;k++ )
		{  B(ig,k)=0;	}
	     
	    }
	  else
	    {
	     
	      for ( int k=1;k<=nr;k++ )
		{  B(ig,k)+=K.area*U(k,G_K,Th)/3; }
	     
	    }
	  
	}      
    }
 


  KN<R> *xx=new KN<R>[nr+1](N);
 
 for (int k=0;k<Th.nv;k++)
    {
      Vertex & v = Th(k);
     
      if (v.lab==1) 	{ xx[0][k]=us;  }  
      if (v.lab==3) 	{ xx[0][k]=un;  }
     
      } 
     
  
  MatriceIdentite<R> I; // Matrice identite

  Mat_lap *AA=new Mat_lap[nr+1](Th);
  bool *converge=new bool[nr+1];

 
  for ( int k=0;k<=nr;k++ )
    { converge[k]=GradienConjugue(AA[k],I, B('.',k),xx[k],N,1e-10); }
 
  char filename0[256];
  char filename1[256];
  char filename2[256];
 
  for ( int k=0;k<=nr;k++ )
    {
      if (converge[k])
	{ // resultat pour gnuplot 
	
	  sprintf(filename0,"sol%d",k);
	  sprintf(filename1,"soll%d",k);
	  sprintf(filename2,"iso%d",k);

	  savesplot(Th,filename0,xx[k]);
	  
	  ofstream fsol(filename1);
	  
	  fsol<< xx[k]; 
	  gnuplot_iso(filename2, Th, xx[k] , 20);
	  
	  } 
      else 
	{
	  cerr << "Erreur : non convergence du gradient conjuge " << endl;
	  return 1; // pour que la commande unix retourne une erreur
	}
    }


 //Résolution du systeme Ax=b avec A(6*6)
 
 KNM<R> A(nr,nr);
 A=0;
 KN<R>  X(nr);
 X=0;
 KN<R>  BB(nr);
 BB=0;
 R Uopt=250.;
 MatriceIdentite<R> I2;



 for (int k=0;k<Th.nt;k++)
   {
     Triangle & K(Th[k]);//Th[k]= triangle k
     
     if ( K.lab == 0 )
       {
	 for ( int il=0;il<=2;il++)
	   {
	     int ig1 = Th(k,il);
	     //Th(k,il)=sommet de numero local il du triangle k	 
	     
	     Vertex pt(Th(ig1));//nomme pt le sommet ig
	   

	     for ( int i=0;i<nr;i++ )
	       {
		 BB[i]+=xx[i+1][ig1]*(Uopt-xx[0][ig1])*K.area/3;
		 for ( int j=0;j<nr;j++ )
		   {  A(i,j)+=xx[i+1][ig1]*xx[j+1][ig1]*K.area/3;   }
	       }
	   }
       }
    }
 
 
 GradienConjugue(A,I2,BB,X,6,1e-10);
  
 cout<<"les valeurs optimales obtenues numériquement sont:"<<endl; 
 for (int i=0;i<nr;i++ )
   
   {   cout<<"résistance "<<i+1<<"="<<X(i)<<endl; }

   
 return 0;
	 
}
