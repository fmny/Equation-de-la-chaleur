//Réponse à la question Q5c)


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



R g(const R2 & P)
{ 
  if ( P.y == H/2  ) return un ;
  if ( P.y == -H/2 ) return us ;
      else return 0.;
}

//R g(const R2 & P){return P.x*P.x+2*P.y*P.y;}  // boundary condition

//R f(const R2 & P){return 0;} // right hand side  $ - \Delta g = -2 -4 $ 

R f (const R2 &P)
{
  if( ( ( sqrt(    (P.x-xr)*(P.x-xr) + (P.y-yr )* (P.y-yr)   )< (H-h)/20 ) )
    ||( (sqrt(    (P.x)*(P.x) + (P.y-yr )* (P.y-yr)   )< (H-h)/20 ) )
      ||((sqrt(    (P.x+xr)*(P.x+xr) + (P.y-yr )* (P.y-yr)   )< (H-h)/20 ) )
  ||((sqrt(    (P.x-xr)*(P.x-xr) + (P.y+yr )* (P.y+yr)   )< (H-h)/20 ) )
  ||((sqrt(    (P.x)*(P.x) + (P.y+yr )* (P.y+yr)   )< (H-h)/20 )) 
||((sqrt(    (P.x+xr)*(P.x+xr) + (P.y+yr )* (P.y+yr)   )< (H-h)/20 )) 
  )

 return 25000.;
  else return 0.;
}

R f1 (const R2 &P)
{
  return 0;
}


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
		  Ax[ig0] += conductivite * K.area * (x[ig0] * (K.H(il) , K.H(il)) + x[ig1] * (K.H(il) , K.H((il+1) %3)) + x[ig2] * (K.H(il) , K.H((il+2) %3)) );
		}
	      
	}
    }
  
  }; 
  plusAx operator*(const KN<R> &  x) const {return plusAx(this,x);}
};



int main(int argc , char** argv )
{
  Mesh Th(argc>1? argv[1] : "four.msh");
  int N=Th.nv;
  
  // construction de b
  KN<R> b(N);
  b=0; 
  

  // Assemblage de la matrice et du second membre 
  for (int k=0;k<Th.nt;k++)
    {
      Triangle & K(Th[k]);//Th[k]= triangle k
       
      R2 sA(K[0]),sB(K[1]),sC(K[2]);  //coord des sommets A,B,C du triangle k
      
      R2 gw[3]={K.H(0),K.H(1),K.H(2)}; //que vaut K.H(i)??
      
      //K.H(0)=gradient 
      
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
  
  // construction de x qui contient la solution du pb avec une resistance qui 
  //chauffe et un=50 et un=100
  //x1=sol avec une res et un=us=0

   KN<R> x(N);   
   x=0; 
 
 for (int k=0;k<Th.nv;k++)
    {
      Vertex & v = Th(k);
      if (v.lab==1) 	{x[k]=us;  }  
      if (v.lab==3) 	{x[k]=un;  }
    } 
     

 
  
  MatriceIdentite<R> C1; // Matrice identite

  Mat_lap A(Th);

  bool converge =GradienConjugue(A,C1, b,x,N,1e-10);       
 
 
  if (converge)
    { // resultat pour gnuplot 
      savesplot(Th,"x",x);
      ofstream fsol("x.sol");
      fsol << x; 
      gnuplot_iso("iso", Th, x , 20);
    } 
  else 
    {
      cerr << "Erreur : non convergence du gradient conjuge " << endl;
      return 1; // pour que la commande unix retourne une erreur
    }
 return 0;
 
}
