#include <string> 
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include "fonction0.hpp"
#include "sfem.hpp"
#include "RNM.hpp"
#include "GC.hpp"
#include "gnuplot-iso.hpp"



using namespace std;
//la grammaire des expressions:
//  exp = term | term '+' exp | '-' exp
//  term =  factor | factor '*' term  | factor '/' term  
//  factor = NUM  | '('exp')' |  FUNC'('exp')'| VAR |'-'factor|'+'factor ; 
//  


extern double * ptx,*pty,*Rr;  // pour le lien avec l'algebre de fonction (a suivre)
double *ptx,*pty,*Rr;
double H=2.;
double  h=0.4;
int nr=6;
double us=100.;
double un=50.;
double k1=1.;
double k2=10.;
double q1=25000.;




void savesplot(const Mesh & Th,char * fileplotname,KN_<R> &x)
{
  ofstream plot(fileplotname);
  for(int it=0;it<Th.nt;it++)
    plot << (R2) Th[it][0] << " " << x[Th(it,0)] <<  endl 
     << (R2) Th[it][1] << " " << x[Th(it,1)] << endl 
     << (R2) Th[it][2] << " " << x[Th(it,2)] << endl 
     << (R2) Th[it][0] << " " << x[Th(it,0)] << endl << endl << endl; 
}



// changement de const Fonction0 en Fonction0
double fr ( Fonction0 & f, R2 P ,R2 C )
{
  R2 X=P-C;
  *ptx=X.x;
  *pty=X.y;
  return f();
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
		  Ax[ig0] += conductivite * K.area * 
		    (x[ig0] * (K.H(il) , K.H(il)) 
		     + x[ig1] * (K.H(il) , K.H((il+1) %3)) 
		     + x[ig2] * (K.H(il) , K.H((il+2) %3)) );
		}
	}
    }
  
  }; 
  plusAx operator*(const KN<R> &  x) const {return plusAx(this,x);}
};


class Expf {  public: 
  typedef double Real;
  typedef Fonction0  Fonc ;
private:
  int sym;
  static const int NUM=257; // codage pour le terminal NUM (nombre)
  static const int FUNC=258;// idem pour les fonctions
  static const int VAR=259; // idem pour les variables
  
  istream & cin;
  Real valeur, valeuro;
  Real (*f1)(Real), (*f1o)(Real); //  si sym == FUNC
  Real *pv, *pvo;     		// Si sym == VAR
  
public: 
  Expf(istream & f) : cin(f) { NextSym();} // lit le 1er symbol 
  void fatal(string s){ cerr <<  s  << endl;exit(1);}  // erreur fatale
 
  void NextSym() {
    pvo=pv;		// on sauve les valeurs associées au symbole 
    f1o=f1;
    valeuro=valeur; 
    
    int c=cin.peek();
    if (isdigit(c) || c == '.') { sym=NUM; cin >> valeur; 
    //   cout << char(c) << " "<< valeur << " " << cin.eof() << " " <<  cin.good() <<endl;
     if(!cin.eof()) assert( cin.good());}
    else if (c=='+' || c=='*' || c=='(' || c==')'  || c==EOF) sym=cin.get();
    else if (c=='-' || c=='/' || c=='^' ) sym=cin.get();
    else if (isspace(c)) {sym=' '; cin.get();}
    else if (isalpha(c)) { 
     string buf;
     while (isalnum(cin.peek()))
       buf += cin.get();
     if ( buf == "exp" ) { sym=FUNC; f1=exp; }
     else if ( buf == "cos" ) { sym=FUNC; f1=cos;}
     else if ( buf == "x" ) { sym=VAR; pv=ptx;}
     else if ( buf == "pi" ) { sym=NUM; valeur=M_PI;}
     else if ( buf == "y" ) { sym=VAR; pv=pty;}
     else if ( buf == "R" ) { sym=VAR; pv=Rr; }//valeur = (H-h)/20
     else {
       cerr << " mot inconnue '" << buf <<"'" << endl; exit(1);}
     cout << "Id " << buf << " " << sym << endl;
    }   
    else {  cerr << char(c) << " (" << c<<") "; fatal("caractere invalide");}
  }
  
  void Match(int c) { 
    if( c!=sym ) {
      if( c==NUM ) cerr <<" un nombre "; else cerr << "'"<<char(c)<<"'" ; 
      fatal(" etait attendu"); }
    else NextSym();
 }
  
  bool ifSym(int c) { if (c==sym) { NextSym();return true;} else return false;} 
  
  bool factor(Fonction0 & v) {   
   if(ifSym(NUM)) {v=valeuro;   return true;}
   if(ifSym('-')) {
     if (!expr(v) ) fatal("dans factor on attendait une expression");
		v = -v;return true;}
   if(ifSym('+')) {
     if (!expr(v) ) fatal("factor: On attendait une expression apres +");
     return true;}
   else if (ifSym('(')) { 
     if (!expr(v) ) fatal(" (exp) : On attendait une expression ");
     Match(')');  
     return true;}
   else if (ifSym(VAR)){ v = pvo; return true; }
   else if (ifSym(FUNC)) {
     Real (*ff)(Real ) = f1o;
     Match('(');
     if (!expr(v) ) fatal(" On attendait une expression ");
      v=compose(ff,v); Match(')');    
      return true;}
   else return false;
  }
  
  bool term(Fonction0 & v) { 
    if ( ! factor(v) ) return false;
    else if (ifSym('*'))
      { 
	Fonction0 vv(0.);
	if (!term(vv) ) fatal(" On attendait un terme ");
	v=v*vv;}
    else if (ifSym('/')) 
      { 
	Fonction0 vv(0.);
	if (!term(vv) ) fatal(" On attendait un terme ");
	v=v/vv;}
    return true;
  }
  
  bool expr(Fonction0 & v) { 
    if ( !term(v) ) { cout << " expression fausse " << endl;return false;}
    else if (ifSym('+')) 
      {
	Fonction0 vv(0.);
	if (!expr(vv) ) fatal(" On attendait une expression ");
	v+=vv;}
    else if (ifSym('-')) 
      { 
	Fonction0 vv(0.);
	if (!expr(vv) ) fatal(" On attendait une expression ");
	v-=vv;}
    return true;
  }
};




int main(int argc,char **argv)
{
  double x,y,r;
 
  //ajout du fichier lecture.cpp au fichier interprete
  string psi;
  double cx,cy,temperaturecible;
 
  int cas;
   
  ifstream file("project.txt");

  assert(file);
  // lecture "  ....  "  complique mais universel 
  // recherche de la premiere "
  while( file.get()!='"' )
    {   assert(file.good());  }
    psi = ""; 
    // recherche de la derniere  "
    while( file.peek()!='"' )
      {
	psi += file.get();  // ajoute du caractere a la chaine 
      }
    file.get(); // eat " 
    //file2.get();
    assert(file.good());    // fichier ok ?
    // fin lecture "  ....  "  complique mais universel 
    file >> cx >> cy >> temperaturecible >> cas >> r;
    assert(file.good());    
   
     
 
    //psi=chaine de caractere (cx,cy)=coord du centre des res. 
    //tempcible=temperature
    //cout << psi << endl;
    //cout << cx << " " << cy << " " << temperaturecible<< " " << cas << endl;
    

    // un pseudo fichier string
    istringstream filepsi(psi);
    
    //  impression du pseudo fichier
    while ( filepsi.peek() != EOF)
    cout << (char) filepsi.get() ;
    cout << endl; 
    
    ifstream file2("project2.txt"); 
    
    Expf::Fonc v(0.);
    ptx=&x;
    pty=&y;
    Rr=&r;   
        
    //on sauvegarde psi(x,y) dans une fonction0
    //k indice de la resistance allumée
    //puis faire pour tout sommet du maillage:
    //pour k=1,..,nr
    //f(P,k)=q_k*fr(f0,R2 sommet,C_k)
    
    
    
    Mesh Th(argc>1? argv[1] : "four.msh");
    
    int N=Th.nv;
    int M=Th.nt;    

    //matrice contenant les valeurs de f sur les sommets
    KN<R> *U=new KN<R>[nr](N);

   
    R2 *C=new R2[nr]; 
    R2 P;
    
    Expf myexp(file2);
    

    
    C[0].x=-cx;    C[0].y=cy;
    C[1].x=0;      C[1].y=cy;
    C[2].x=cx;     C[2].y=cy;
    C[3].x=-cx;    C[3].y=-cy;
    C[4].x=0;      C[4].y=-cy;
    C[5].x=cx;     C[5].y=-cy;
    
    
    cout<<"cx="<<cx<<endl;
    cout<<"cy="<<cy<<endl;
    cout<<"Rr="<<r<<endl;    
    
    while (true) {
      
      if (  myexp.expr(v) ) 
	{
	  
	  for (int k=0;k<N;k++)
	    {	 
	      //definition du point courant P qui va décrire 
	      //tous les sommets de la triangulation
	      x=Th(k).x;
	      y=Th(k).y;
	      P.x=x; P.y=y;
	      

	      for (int t=0;t<nr;t++ )
		U[t][k]=fr(v,P,C[t]);//definition de la fonction f
	    }
	
	  // cout << "\n Autre exp : ";   
	  //myexp.Match(' '); 
	}
      else  break;
    }
    cout << " On a fini "<< endl;
    
    savesplot(Th,"1",U[0]);
    savesplot(Th,"2",U[1]);
    savesplot(Th,"3",U[2]);
    savesplot(Th,"4",U[3]);
    savesplot(Th,"5",U[4]);
    savesplot(Th,"6",U[5]);
    
    
    
      KN<R> somme(N);
      somme=0;
      
      for (int k=0;k<nr;k++)
      somme+=U[k];
      
      savesplot(Th,"7",somme);
    
   
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

	
	
	for ( int il=0;il<=2;il++)
	{
	  int ig = Th(k,il);//Th(k,il)=sommet global du triangle de numero 
	  //local il du triangle k	 
	  
	  Vertex pt(Th(ig));//nomme pt le sommet de numero global ig
	  
	  if (  ( pt.lab ==1) ||( pt.lab ==3 )  )
	    { 
	      for ( int k=0;k<=nr;k++ )
		{  B(ig,k)=0;	}
	    }
	  else
	    {
	   
	       for ( int r=1;r<=nr;r++ )
		 
		 //////////////////////////////////////
		 //{  B(ig,r)+=K.area*U2(r,G_K,Th)/3; }//marche donc pb avec la matrice
	       ////////////////////////////////////////
		     
		 { B(ig,r)+=K.area*fr(v,G_K,C[r-1])/3; }
      
	    }
	  
	}
    }
    
    //cout<<"B="<<B<<endl;
  
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
	{ 
	  // resultat pour gnuplot 
	
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



  //solution du pb direct
  KN<R> xxx(N);
  xxx=xx[0];

  if (converge[0]&&converge[1]&&converge[2]&&converge[3]&&converge[4]&&converge[5]&&converge[6])
    {
      for (int k=1;k<=nr;k++)
	{

	  for (int l=0;l<N;l++ )
	    { xxx[l]+=q1*xx[k][l]; }
	}
	  savesplot(Th,"solution0",xxx);
	  
	  ofstream fsol2("solution1");
	  
	  fsol2<<xxx; 
	  gnuplot_iso("solution2", Th, xxx , 20);
	  
  
    }
  else
    {
      cerr << "Erreur : non convergence du gradient conjuge " << endl;
      return 1; // pour que la commande unix retourne une erreur
    }
  
  

 //Résolution du systeme Ax=b avec A(6*6)
   
 KNM<R> A(nr,nr);
 A=0;
 KN<R>  X(nr);
 X=0;
 KN<R>  BB(nr);
 BB=0;
 // R Uopt=250.;
 MatriceIdentite<R> I2;
 R Uopt=temperaturecible;
 cout<<"Uopt="<<Uopt<<endl;

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
