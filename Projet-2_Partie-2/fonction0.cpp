#include "fonction0.hpp"

using namespace std;  //introduces namespace std


int main()
{	
  R *xxx= new R;
  R *yyy= new R;
  R *rrr= new R;
  cout << " ---" << endl;
  {
  Fonction0 X(xxx),Y(yyy),R(rrr);

  Fonction0 psi(R*R*compose(exp,-(X*X+Y*Y)/R));
  *xxx=1;
  *yyy=2;
  *rrr=3;
  cout  << psi() << endl;
  *xxx=1;
  *yyy=1;
  *rrr=0.1;
  cout  << psi() << endl;
  }
  cout << " ---" << endl;

  delete xxx;
  delete yyy;
  delete rrr;

  return 0;
}



