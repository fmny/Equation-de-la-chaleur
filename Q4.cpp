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

// Resolution de l'EDP: -Delta u =f sur Omega 
// avec les conditions de Dirichlet : u=g sur Gamma

double L=2., H=2., l=1., h=0.4, k1=1., k2=10., us=100., un=50.;
double xr=(-l/2); 
double yr=(H+h)/4;



//R g(const R2 & P){return P.x*P.x+2*P.y*P.y;}  
//R f(const R2 & ){return -6.;}

/*
R g (const R2 &P, Mesh &Th)
{
  for ( int k=0; k<Th.nv; k++ )
    { 
      for ( int il=0; il<3; il++ )
	Vertex &V=Th[k][il];
      if ( ( V.x == P.x ) && ( V.y == P.y ) )
	{
	  if ( V.lab == 1 ) { return un; }
	  if ( V.lab == 3  ) { return us; }
	  }
	  else return 0.; 
	  
	  
	  }
*/

R g(const R2 & P)
{ 
  if ( P.y == H/2  ) return un ;
  if ( P.y == -H/2 ) return us ;
      else return 0.;
}

//R g(const R2 & P){return P.x*P.x+2*P.y*P.y;}  // boundary condition

//R f(const R2 & P){return 0;} // right hand side  $ - \Delta g = -2 -4 $ 

R f (const R2 &P)
{ return 0.;}


// fonction pour sauvegarder la solution 
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
		  Ax[ig0] += conductivite * K.area * (x[ig0] * (K.H(il) , K.H(il)) 
							 + x[ig1] * (K.H(il) , K.H((il+1) %3)) +
							 x[ig2] * (K.H(il) , K.H((il+2) %3)) );
		}
	      
	}
    }
  
  }; 
  plusAx operator*(const KN<R> &  x) const {return plusAx(this,x);}
};
 

/*
//resolution du systame A.x=b avec A(3*3)
void SystemeOrdre3( KNM<R> A, R2 &P,KN<R> &u)
{

 MatriceIdentite<R> I;
 KN<R> b(3);
 KN<R> x(3);

 b[0]=P.x; b[1]=P.y; b[2]=1.;


 GradienConjugue(A,I,b,x,3,1e-10);

 u[0]=x[0]; u[1]=x[1]; u[2]=x[2];
}
*/


/*
int recherchelabel ( R2 &P , Mesh &T )
{
  int labelPoint;
  Mesh &Th(T);
  for ( int k=0;k<Th.nt;k++ )
    {
      Triangle & K(Th[k]);
      
    
      R2 i0(K[0]),i1(K[1]),i2(K[2]); //coordonnées des 3 sommets 
 
      KN<R> v(3); 
      KNM<R> A(3,3);
 
      A(0,0)=i0.x; A(0,1)=i1.x; A(0,2)=i2.x;
      A(1,0)=i0.y; A(1,1)=i1.y; A(1,2)=i2.y;
      A(2,0)=1.;   A(2,1)=1.;   A(2,2)=1.;
           
      SystemeOrdre3(A,P,v);



      if( ( 0 <= v[0])&&( 0 <= v[1]  )&&( 0 <= v[2]  ) )
	{  labelPoint=K.lab; return labelPoint;   }
    }
}
*/     

/*
int recherchelabel ( R2 &P , Mesh &Th )
{   
  int labelpoint;  
  double lambda0,lambda1,lambda2;
  
  while( (  lambda0 >= 0  )&&(  lambda1 >=0 )&&( lambda2 >= 0 ) )
    {
      for ( int k=0;k<=Th.nt;k++ )


	{
	if ( k == Th.nt ) { return 10000; }	  


	  Triangle & K(Th[k]);
	  
	  R2 i0(K[0]),i1(K[1]),i2(K[2]); //coordonnées des 3 sommets 
	  
	  
	  double a,b,c,d;
	  a=i0.x-i2.x;
	  b=i1.x-i2.x;
	  c=i0.y-i2.y;
	  d=i1.y-i2.y;
	  
	  if (  a != 0   )
	    {
	      lambda1=(P.y-(P.x-i2.x)/a-i2.y)*(1./(d-c*b/a));
	      lambda0=(P.x-lambda1*b-i2.x)*(1./a);
	      lambda2=1-lambda0-lambda1;
	    }
	  else 
	    {
	      lambda1=(P.x-i2.x)/b;
	      lambda0=(P.y-lambda1*d-i2.y)/c;
	      lambda2=1-lambda0-lambda1;	   
	    }
	  
	  if (  (lambda0 <= 1 )&&( lambda0 >=0)  )
	    cout<<"lambda="<<lambda0<<"  "<<lambda1<<"  "<<lambda2<<endl;
	  
	  labelpoint=K.lab;	  
	}
      
    }
  return labelpoint;
}
*/

/*
  int recherchelabel ( R2 &P , Mesh &Th )
{          
  int labelpoint;  
  double lambda0,lambda1,lambda2;

  while( (  lambda0 >= 0  )&&(  lambda1 >=0 )&&( lambda2 >= 0 ) )
   {
     for ( int k=0;k<Th.nt;k++ )
       
       {
	 Triangle & K(Th[k]);
	 
	 R2 i0(K[0]),i1(K[1]),i2(K[2]); //coordonnées des 3 sommets 
	 
	 cout<<"coord"<<i0.x<<"  "<<i0.y<<"   "<<i1.x<<"   "<<i1.y
	     <<"   "<<i2.x<<"  "<<i2.y<<endl;
	 
	 double lambda0,lambda1,lambda2;
	 
	 double a,b,c,d;
	 a=i0.x-i2.x;
	 b=i1.x-i2.x;
	 c=i0.y-i2.y;
	 d=i1.y-i2.y;
	 
	 {
	   if(  ( b*c-a*d != 0   )&&( a != 0) )
	     {
	       lambda1=(c*P.x-a*P.y-c*i2.x+a*i2.y)/(b*c-a*d)     ;
	       lambda0=(P.x-lambda1*b-i2.x)*(1./a);
	       lambda2=1-lambda0-lambda1;
	     }
	   if( ( b*c-a*d != 0 )&& ( a = 0 )&& (d != 0  ))
	     {
	       lambda1=(P.x-i2.x)/b;
	       lambda0=(P.y-lambda1*d-i2.y)/d;
	       lambda2=1-lambda0-lambda1;	   
	     }
	   
	   if( ( b*c-a*d != 0 )&&( a == 0 )&&( d == 0) )
	     {
	       lambda1=(P.x-i2.x)/b;
	       lambda1=(P.y-i2.y)/c;
	       lambda2=1-lambda1-lambda0;
	     }
	   else { lambda0=1000; lambda1=1000; lambda2=1000; }      
	 }    
	 
	 if (  (lambda0 <= 1 )&&( lambda0 >=0)  )
	   cout<<"lambda="<<lambda0<<"  "<<lambda1<<"  "<<lambda2<<endl;
	 labelpoint=K.lab;
       }
     {  return labelpoint; }
     //      else { cout<<b*c-a*d<<endl; return 10000; };  
   }
}
*/



int main(int argc , char** argv )
{
  Mesh Th(argc>1? argv[1] : "four.msh");
  int N=Th.nv;
  
  // construction de b
  KN<R> b(N);
  b=0;
  
  /*
  R2 W1(0.01,0);
  cout<<"label de ("<<W1.x<<","<<W1.y<<")="<<recherchelabel(W1,Th)<<endl;
  */

  // Assemblage de la matrice et du second membre 
  for (int k=0;k<Th.nt;k++)
    {
      Triangle & K(Th[k]);//Th[k]= triangle k
      //      cout<<"     "<<K.lab<<endl;
      
      //      cout<<"K.lab="<<K.lab<<endl;
      
      R2 sA(K[0]),sB(K[1]),sC(K[2]);  //coord des sommets A,B,C du triangle k
      
      R2 gw[3]={K.H(0),K.H(1),K.H(2)}; //que vaut K.H(i)??
      
      //K.H(0)=gradient ???
      
      //cout<<K.H(0)<<"      "<<K.H(1)<<"        "<<K.H(2)<<endl;
      R2 G_K((sA+sB+sC)/3.);  //G_K barycentre du triangle k
      // R intfgkwi= f(G_K)*K.area/3; //  $ \int_K f w_i = f(G_K) |K|/3 $
      //intfgkwi=valeur de f evaluée en le
      // barycentre du triangle k *aire du triangle/3 
      
      for ( int il=0;il<=2;il++)
	{
	  int ig = Th(k,il);  //Th(k,il)=sommet de numero local il du triangle k	 
	  Vertex pt(Th(ig));   //nomme pt le sommet ig
	  
	  if (  ( pt.lab ==1 || pt.lab ==3 ) )
	    { b[ig]=0; }
	  else
	    { b[ig]+= K.area*f(G_K)/3; }
	}
    }
  
  // construction de x
   KN<R> x(N);   
  x=0;  
  
  for (int k=0;k<Th.nv;k++)
    {
      Vertex & v = Th(k);
      if (v.lab==1) 	{x[k]=us;}  
      if (v.lab==3) 	{x[k]=un;}
    } 
  
  
  MatriceIdentite<R> C1; // Matrice identite

  Mat_lap A(Th);

  bool converge=GradienConjugue(A,C1, b,x,N,1e-10);       



  //recherche des labels des triangles
  /*
  for (int k=0;k<Th.nt;k++)
    {
      Triangle & K(Th[k]);
      if (  ( K.lab == 14 )||( K.lab == 6)||(K.lab == 0)  )
	{
	  for ( int il=0;il<3;il++ ) 
	    x[Th(k,il)]=1000;
	}
    }
  */
 
  if (converge)
    { // resultat pour gnuplot 
      savesplot(Th,"q2b.splot",x);
      ofstream fsol("x.sol");
      fsol << x; 
      gnuplot_iso("fileiso_q2b", Th, x , 25);
    } 
  else 
    {
      cerr << "Erreur : non convergence du gradient conjuge " << endl;
      return 1; // pour que la commande unix retourne une erreur
    }
 
  return 0;
}
