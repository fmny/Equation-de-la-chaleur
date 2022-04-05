#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include "sfem.hpp"
#include "RNM.hpp"
#include "GC.hpp"
#include "gnuplot-iso.hpp"


double H=2;
double h=0.4;
double l=1.0;
double xr=(-l/2); 
double yr=(H+h)/4;
double un=50.;
double us=100.;

/*
R g(const R2 & P)
{ 
  if ( P.y == H/2  ) return un ;
  if ( P.y == -H/2 ) return us ;
      else return 0.;
}
*/

R g(const R2 & P){return P.x*P.x+2*P.y*P.y;}  // boundary condition
R f(const R2 & P){return -6;} // right hand side  $ - \Delta g = -2 -4 $ 
/*
R f (const R2 &P)
{
  if ( sqrt(    (P.x-xr)*(P.x-xr) + (P.y-yr )* (P.y-yr)   )< (H-h)/20 ) 
     return 25000.; 
  else return 0.;
}
*/


void savesplot(const Mesh & Th,char * fileplotname,KN_<R> &x)
{
  ofstream plot(fileplotname);
  for(int it=0;it<Th.nt;it++)
    plot << (R2) Th[it][0] << " " << x[Th(it,0)] <<  endl 
	 << (R2) Th[it][1] << " " << x[Th(it,1)] << endl 
	 << (R2) Th[it][2] << " " << x[Th(it,2)] << endl 
      //	 << (R2) Th[it][0] << " " << x[Th(it,0)] << endl 
	 << endl << endl; 
}


// une classe matrice diagonale pour le preconditionneaur diagonal
template <class R>
class MatDiag: VirtualMatrice<R> { public:
  const KN_<R> d; // Adr du Vecteur diagonal
  typedef typename  VirtualMatrice<R>::plusAx plusAx;
  
  MatDiag(KN_<R> dd) :d(dd) {}

  void addMatMul(const  KN_<R>  & x, KN_<R> & Ax) const { 
    assert(x.N()==Ax.N() && x.N() == d.N() );
    Ax+=DotStar_KN_<R>(x,d);   // <=>  for(int i=0;i<x.n;i++) Ax[i] += x[i]*d[i];
  }
  
  plusAx operator*(const KN<R> &  x) const {
    return plusAx(this,x);}
};

int main(int argc , char** argv )
{
  Mesh Th(argc>1? argv[1] : "carre.msh");
  int N=Th.nv; 
  KN<R> b(Th.nv);  //  $ b[i] =  \int_\Omega  f w_i \sim \sum_{K\in\mathcal{T}_h} f(G_K)\int_K w_i $
  KNM<R> A(N,N); // Matrice pleine stupide mais pour l'exemple
  b=0; // second membre 
  A =0; // matrice 
  // Assemblage de la matrice et du second membre 
  for (int k=0;k<Th.nt;k++)
    {
      Triangle & K(Th[k]);//Th[k]=triabgle k
  
      cout<<"label"<<K.lab<<endl;

      //cout<<Th[k]<<endl;
      R2 sA(K[0]),sB(K[1]),sC(K[2]);  //coord des sommets A,B,C du triangle k
      
      //cout<<sA<<endl;

      R2 gw[3]={K.H(0),K.H(1),K.H(2)}; //que vaut K.H(i)??

      //K.H(0)=gradient ???

      //cout<<K.H(0)<<"      "<<K.H(1)<<"        "<<K.H(2)<<endl;
      R2 G_K((sA+sB+sC)/3.);  //G_K barycentre du triangle k
      R intfgkwi= f(G_K)*K.area/3; //  $ \int_K f w_i = f(G_K) |K|/3 $
      //intfgkwi=valeur de f en le barycentre du triangle k *aire du triangle/3 

      ////////////////////////////////////////////

      for (int il=0;il<3;++il)
	{ 
	  int i=Th(K[il]); //  numero de sommet il du triangle k
	  b[i] += intfgkwi;  // Assemblage  du second membre 
	  for (int jl=0;jl<3;++jl) // Assemblage de la matrice 
	    { 
	      int j=Th(K[jl]); //  sommet jl du triangle k   
	      R a_K_ij= (gw[il],gw[jl])*K.area;
	      A(i,j) += a_K_ij;
	    }
	}
    }

  ////////////////////////////////////////////


  R tgv = 1e30; // tres grande valeur pour la penalisation exacte.
  
  KN<R> x(N),xe(N);   
  x=0;  //  donnee initial,   il faut mettre les conditions aux limites.
  
  // Prise en compte des condition au limite de diriclet par penalisation exact.
  KN<R> D1(N); // un vecteur pour stocker la 1/diagonale de la matrice A
  for (int i=0;i<Th.nv;i++)
    { 
      Vertex & Si(Th(i)); 
      xe[i]=g(Si);
      if (Si.onGamma()) 
	{
	  A(i,i) = tgv;  
	  b[i]= tgv*g(Si);  
	  x[i] = g(Si); // la donne initial verifie les CL (mieux)
	}
      D1[i]=1./A(i,i);  //  Construction du vecteur 1/diagonal
    }
  
  MatDiag<R> C(D1); // Matrice de preconditionnement (Dii)^-1
  bool converge=GradienConjugue(A,C, b,x,N,1e-10);       
  //necessite la matrice A en argument
  
  if (converge)
    { // a file for gnuplot 
      savesplot(Th,"x.splot",x);
      savesplot(Th,"xe.splot",xe);
      xe -= x;
      cout << " diff x-xe :  min =" << xe.min() << " , max = " << xe.max() << endl;
      ofstream fsol("x.sol");
      fsol << x; 

    } 
  else 
    {
      cerr << "Erreur no converge du gradient conjuge " << endl;
      return 1; // pour que la commande unix retour une erreur
    }
  gnuplot_iso ("fileiso",Th,x,20);
  return 0;
}


