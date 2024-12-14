#include <iostream> 
#include <iomanip>
#include <string>
using namespace std;

#include <imatrix.hpp>
using namespace cxsc;

#include "sys/time.h"
double time() {
   struct timeval _tp;
   gettimeofday(&_tp,0);  
   return _tp.tv_sec + _tp.tv_usec / 1000000.0;
}

int n(100000);    //Number of summands/repetitions, GLOBAL variable
                 //and length of vectors.  
int repMax(100);  //Number of loop repetitions
double baseTime; //Time used for dotproduct of IEEE double vectors, GLOBAL
double diffTime; //Only used in print...

void print_time_penalty(const double startTime) {
  if (baseTime <1e-10) cout << "*** baseTime too small! " << endl; 
  diffTime= time()-startTime;
  cout.setf(ios::showpoint|ios::fixed, ios::floatfield);
  cout.precision(1);  
  cout << "Penalty " << setw(6) << diffTime/baseTime;
  cout.precision(10);
  cout << ",  time used  " << diffTime << endl;
} 

int main() {
  double startTime;
  cout << endl;
  cout << "*****************************************" << endl;
  cout << "Time measurements for different kinds of " << endl;
  cout << "dot product computations" << endl;
  cout << "Vector length: " << n << endl;
  cout << "repMax: " << repMax << endl;
  
  double ds(0), dx, dy;
  dx= 1.0/3;
  dy= 7.0/19;
 
  //warm up:
  for (int rep=1; rep<repMax; rep++)
  for (int i=1; i<10000000; i++) {
    ds+= dx*dy;
  }  
  double h(ds); //use ds to avoid compiler optimization
  ds= h-h;      //clear ds

  cout << "baseTime using IEEE double: ";  
  startTime= time();  
  for (int rep=1; rep<repMax; rep++)
  for (int i=1; i<n; i++) {
    ds+= dx*dy;
  }  
  baseTime= time()-startTime;
  cout << baseTime << endl;
  cout << "ds: " << ds << endl;

  //DOUBLE data, ordinary loop
  ds= 0;
  cout << "1)  Double ds+= dx*dy: " << endl;  
  startTime= time();  
  for (int rep=1; rep<repMax; rep++)
  for (int i=1; i<n; i++) {
    ds+= dx*dy;
  }  
  print_time_penalty(startTime); 
  cout << "ds: " << ds << endl;
  
  //DOUBLE array
  double dax[n], day[n];
  for (int i=0; i<n; i++) {
    dax[i]= 1.0/3;
    day[i]= 7.0/19;
  }  
  ds= 0; 
  
  cout << "2)  Double array ds+= dax[i]*day[i]: " << endl;  
  startTime= time();  
  for (int rep=1; rep<repMax; rep++)
  for (int i=0; i<n; i++) {
    ds+= dax[i]*day[i];
  }  
  print_time_penalty(startTime);
  cout << "ds: " << ds << endl;

  double *dnewx, *dnewy;
  dnewx= new double[n];
  dnewy= new double[n];
  for (int i=0; i<n; i++) {
    dnewx[i]= 1.0/3;
    dnewy[i]= 7.0/19;
  }  
  ds= 0; 

    
  cout << "2b)  Double arrays created with new: " << endl;  
  startTime= time();  
  for (int rep=1; rep<repMax; rep++)
  for (int i=0; i<n; i++) {
    ds+= dnewx[i]*dnewy[i];
  }  
  print_time_penalty(startTime);
  cout << "ds: " << ds << endl;

  //REAL data
  real rs(0), rx, ry;
  rx= 1.0/3;
  ry= 7.0/19;
  
  cout << "3)  Real rs+= rx*ry: " << endl;    
  startTime= time(); 
  for (int rep=1; rep<repMax; rep++)   
  for (int i=1; i<n; i++) {
    rs+= rx*ry;
  }  
  print_time_penalty(startTime);
  cout << "rs: " << rs << endl;
  
  //REAL data using rvector
  rvector rvx(n), rvy(n);
  rvx= real(1.0)/3;
  rvy= real(7.0)/19;
  rs= 0;
  
  cout << "4)  Real using rvector rs+= rvx[i]*rvy[i]: " << endl;    
  startTime= time(); 
  for (int rep=1; rep<repMax; rep++) 
  for (int i=1; i<n; i++) {
    rs+= rvx[i]*rvy[i];
  }  
  print_time_penalty(startTime);  
  cout << "rs: " << rs << endl;
  
  //INTERVAL data
  interval is(0), ix, iy;
  ix= interval(1.0)/3;
  iy= interval(7.0)/19;
  
  cout << "5)  Interval is+= ix*iy: " << endl;    
  startTime= time(); 
  for (int rep=1; rep<repMax; rep++) 
  for (int i=1; i<n; i++) {
    is+= ix*iy;
  }  
  print_time_penalty(startTime);
  cout << "is: " << is << endl;
 
  //DOTPRECISION data
  dotprecision rdots(0);
  rx= 1.0/3;
  ry= 7.0/19;
  
  cout << "6)  Dotprecision rdots+= rx*ry: " << endl;  
  startTime= time(); 
  for (int rep=1; rep<repMax; rep++)  
  for (int i=1; i<n; i++) {
    rdots+= rx*ry;
  }  
  print_time_penalty(startTime);
  cout << "rdots: " << rdots << endl;
  
  //IDOTPRECISION data
  idotprecision idots(0);
  
  cout << "7)  Idotprecision with idots+= ix*iy: " << endl;  
  startTime= time();  
  for (int rep=1; rep<repMax; rep++) 
  for (int i=1; i<n; i++) {
    idots+= ix*iy;
  }  
  print_time_penalty(startTime);
  cout << "idots: " << idots << endl;
       
  //DotK data with k equal to 2
  idots= 0;
 
 /*  idots.set_k( ) is only in effect when using vector
      arguments when calling accumulate()!
  cout << "8)  Dot2 and accumulate(idots,ix,iy): " << endl;  
  idots.set_k(2);
  startTime= time();  
  for (int i=1; i<n; i++) {
    accumulate(idots,ix,iy);
  }  
  print_time_penalty(startTime);
*/  
        
  //DOTK in real dot product with k equal to 1
  rdots= 0;
  //rvector rvx(n);
  rvx= 1.0/3;
  //rvector rvy(n);
  rvy= 7.0/19;
  
  opdotprec= 1;
  cout << "9)  Dot product using rvectors and Dot1: " << endl;  
  startTime= time(); 
  for (int rep=1; rep<repMax; rep++)  
    rdots+= rvx*rvy;
  print_time_penalty(startTime);
  cout << "rdots: " << rdots << endl;

  //DotK in real dot product with k equal to 2  
  rdots= 0;
  opdotprec= 2;
  rdots.set_k(2);
  
  cout << "10) Dot product using rvectors and Dot2: " << endl;  
  startTime= time();  
  for (int rep=1; rep<repMax; rep++)  
    rdots+= rvx*rvy;
  print_time_penalty(startTime);
  cout << "rdots: " << rdots << endl;
  
  //DotK in real dot product with k equal to 0 (accu)  
  rdots= 0;
  opdotprec= 0;
  rdots.set_k(0);
  
  cout << "11) Dot product using rvectors and Dot0 (accu): " << endl;  
  startTime= time(); 
  for (int rep=1; rep<repMax; rep++)  
    rdots+= rvx*rvy;
  print_time_penalty(startTime);
  cout << "rdots: " << rdots << endl;   
  
  //DotK in interval dot product with k equal to 1
  idots= 0;
  ivector ivx(n);
  ivx= interval(1.0)/3;
  ivector ivy(n);
  ivy= interval(7.0)/19;
  opdotprec= 1;
  idots.set_k(1); //???
  
  cout << "12) Interval dot product using ivectors and Dot1: " << endl;  
  startTime= time();  
  for (int rep=1; rep<repMax; rep++) 
    idots+= ivx*ivy;
  print_time_penalty(startTime);
  cout << "idots: " << idots << endl;
         
  //DotK in interval dot product with k equal to 2
  idots= 0;
  opdotprec= 2;
  idots.set_k(2);
  
  cout << "13) Interval dot product using ivectors and Dot2: " << endl;  
  startTime= time(); 
  for (int rep=1; rep<repMax; rep++)  
    idots+= ivx*ivy; 
  print_time_penalty(startTime);
  cout << "idots: " << idots << endl;
  
  //DotK in interval dot product with k equal to 0 (accu)
  idots= 0;
  opdotprec= 0;
  idots.set_k(0);
  
  cout << "14) Interval dot product using ivectors and Dot0 (accu): " << endl;  
  startTime= time(); 
  for (int rep=1; rep<repMax; rep++)  
    idots+= ivx*ivy; 
  print_time_penalty(startTime);
  cout << "idots: " << idots << endl;
}
