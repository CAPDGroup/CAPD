#include <iostream>
#include <iomanip>
#include <ctime> // Zeitmessung 
#include <cstdlib> // Zufallszahlen rand(), srand()
#include <cmath> 
#include "interval.hpp"
#include "ivector.hpp"
#include "idot.hpp"
#include "rmath.hpp"
#include "imath.hpp"

using namespace std;
using namespace cxsc;

void start_clock(clock_t& t1); // Startet den Timer 
void print_time_used(clock_t t1); 
double return_time_used(clock_t t1); 
   
 
int main(void)
{
   int kMax= 100000000; // Schleifendurchlaeufe zur Zeitmessung Grundoperationen
   int eMax= 5000000;   // Schleifendurchlaeufe zur Zeitmessung element. Fktn.
//   int kMax= 1; // Schleifendurchlaeufe zur Zeitmessung Grundoperationen
//   int eMax= 1;   // Schleifendurchlaeufe zur Zeitmessung element. Fktn.
   int vMax= 20000; // Dimension Vektor fuer Skalarprodukt
   int sMax= 1000;   // Schleifendurchlaeufe zur Zeitmessung Skalarprodukt
   clock_t t; // Datentyp in <ctime> definiert
   double ret_time1,ret_time2,ret_time3,ret_time4;
   
   string
   comp_name,comp_cpu,comp_ram,comp_os,comp_c,comp_copt,comp_cpp,comp_cppopt;
   cout << "Benchmark Program for C-XSC Library" << endl;
   cout << "Please insert some information about your system:" << endl;
   cout << "Name of your computer?  "; getline(cin,comp_name); //cin >> comp_name;
   cout << "Which CPU?              "; getline(cin,comp_cpu); //cin >> comp_cpu;
   cout << "RAM?                    "; getline(cin,comp_ram); //cin >> comp_ram;
   cout << "Operation System (OS)?  "; getline(cin,comp_os); //cin >> comp_os;
   cout << "C-Compiler (version)?   "; getline(cin,comp_c); //cin >> comp_c;
   cout << "Used Compiler options?  "; getline(cin,comp_copt); //cin >> comp_copt;
   cout << "C++-Compiler (version)? "; getline(cin,comp_cpp); //cin >> comp_cpp;
   cout << "Used Compiler options?  "; getline(cin,comp_cppopt); //cin >> comp_cppopt;
   cout << endl << endl;;                


   cout << "Benchmark Program for C-XSC Library - Protocol" << endl;
   cout << "==============================================" << endl << endl;
   cout << "Some information about the computer system: " << endl << endl;;
   cout << "  Computer name: " << comp_name << endl;
   cout << "  CPU:           " << comp_cpu << " with " << comp_ram << " RAM" << endl;
   cout << "  OS:            " << comp_os << endl;
   cout << "  C-Compiler:    " << comp_c << " with options: " << comp_copt << endl;
   cout << "  C++-Compiler:  " << comp_cpp << " with options: " << comp_cppopt << endl;
   cout << endl; 

   cout << endl;
   cout << "Basic Operations:" << endl;
   cout << "=================" << endl << endl;

   double ad, bd, cd, dummyd=0.0;
   ad = 1.23456789e4;
   bd = 4.6536727255666354e-1;
   real ar, br, cr, dummyr=0.0;
   ar = 1.23456789e4;
   br = 4.6536727255666354e-1;
   cout << "Number of operations: " << kMax << "  -  All time measurements in msec." << endl << endl;
 
   cout << "Column:    (1)  |    (2)     |    (3)     |    (4)" << endl;
   cout << "Datatyp: double | cxsc::real | cxsc::real | cxsc::real  factor |factor |factor" << endl;
   cout << "Rnd.:    nearest|   nearest  | downwards  |   upwards   (2)/(1)|(3)/(1)|(4)/(1)" << endl;
   cout << "----------------+------------+------------+-----------  -------+-------+-------" << endl;

   cout << fixed << setprecision(1);
   start_clock(t); 
   for (int k= 0; k< kMax; k++) { ad = ad + bd; }  
   ret_time1=return_time_used(t); 
   cout << "Addition" << setw(7) <<  ret_time1 << " |  ";
   dummyd += ad; ad = 1.23456789e4;

   start_clock(t); 
   for (int k= 0; k< kMax; k++) { ar = ar + br; }  
   ret_time2=return_time_used(t); 
   cout << setw(9) << ret_time2 << " | ";
   dummyr += ar; ar = 1.23456789e4;

   start_clock(t); 
   for (int k= 0; k< kMax; k++) { ar = addd(ar,br); }  
   ret_time3=return_time_used(t); 
   cout << setw(9) << ret_time3 << "  |";
   dummyr += ar; ar = 1.23456789e4;
 
   start_clock(t); 
   for (int k= 0; k< kMax; k++) { ar = addu(ar,br); }  
   ret_time4=return_time_used(t); 
   cout << setw(9) << ret_time4 << "   ";
   dummyr += ar; ar = 1.23456789e4;

   cout << "  " << setw(4)<< ret_time2/ret_time1;
   cout << "  | " << setw(4) << ret_time3/ret_time1;
   cout << "  | " << setw(4) << ret_time4/ret_time1 << endl;

// -----------------------------------------------------------------
 
   start_clock(t); 
   for (int k= 0; k< kMax; k++) { ad = ad - bd; }  
   ret_time1=return_time_used(t); 
   cout << "Subtr.  " << setw(7) <<  ret_time1 << " |  ";
   dummyd += ad; ad = 1.23456789e4;

   start_clock(t); 
   for (int k= 0; k< kMax; k++) { ar = ar - br; }  
   ret_time2=return_time_used(t); 
   cout << setw(9) << ret_time2 << " | ";
   dummyr += ar; ar = 1.23456789e4;

   start_clock(t); 
   for (int k= 0; k< kMax; k++) { ar = subd(ar,br); }  
   ret_time3=return_time_used(t); 
   cout << setw(9) << ret_time3 << "  |";
   dummyr += ar; ar = 1.23456789e4;
 
   start_clock(t); 
   for (int k= 0; k< kMax; k++) { ar = subu(ar,br); }  
   ret_time4=return_time_used(t); 
   cout << setw(9) << ret_time4 << "   ";
   dummyr += ar; ar = 1.23456789e4;

   cout << "  " << setw(4)<< ret_time2/ret_time1;
   cout << "  | " << setw(4) << ret_time3/ret_time1;
   cout << "  | " << setw(4) << ret_time4/ret_time1 << endl;

// -----------------------------------------------------------------


   start_clock(t); 
   for (int k= 0; k< kMax; k++) { ad = ad * bd; }  
   ret_time1=return_time_used(t); 
   cout << "Mult.   " << setw(7) <<  ret_time1 << " |  ";
   dummyd += ad; ad = 1.23456789e4;

   start_clock(t); 
   for (int k= 0; k< kMax; k++) { ar = ar * br; }  
   ret_time2=return_time_used(t); 
   cout << setw(9) << ret_time2 << " | ";
   dummyr += ar; ar = 1.23456789e4;

   start_clock(t); 
   for (int k= 0; k< kMax; k++) { ar = muld(ar,br); }  
   ret_time3=return_time_used(t); 
   cout << setw(9) << ret_time3 << "  |";
   dummyr += ar; ar = 1.23456789e4;
 
   start_clock(t); 
   for (int k= 0; k< kMax; k++) { ar = mulu(ar,br); }  
   ret_time4=return_time_used(t); 
   cout << setw(9) << ret_time4 << "   ";
   dummyr += ar; ar = 1.23456789e4;

   cout << "  " << setw(4)<< ret_time2/ret_time1;
   cout << "  | " << setw(4) << ret_time3/ret_time1;
   cout << "  | " << setw(4) << ret_time4/ret_time1 << endl;

// -----------------------------------------------------------------
   start_clock(t); 
   for (int k= 0; k< kMax; k++) { ad = ad / bd; }  
   ret_time1=return_time_used(t); 
   cout << "Division" << setw(7) <<  ret_time1 << " |  ";
   dummyd += ad; ad = 1.23456789e4;

   start_clock(t); 
   for (int k= 0; k< kMax; k++) { ar = ar / br; }  
   ret_time2=return_time_used(t); 
   cout << setw(9) << ret_time2 << " | ";
   dummyr += ar; ar = 1.23456789e4;

   start_clock(t); 
   for (int k= 0; k< kMax; k++) { ar = divd(ar,br); }  
   ret_time3=return_time_used(t); 
   cout << setw(9) << ret_time3 << "  |";
   dummyr += ar; ar = 1.23456789e4;
 
   start_clock(t); 
   for (int k= 0; k< kMax; k++) { ar = divu(ar,br); }  
   ret_time4=return_time_used(t); 
   cout << setw(9) << ret_time4 << "   ";
   dummyr += ar; ar = 1.23456789e4;

   cout << "  " << setw(4)<< ret_time2/ret_time1;
   cout << "  | " << setw(4) << ret_time3/ret_time1;
   cout << "  | " << setw(4) << ret_time4/ret_time1 << endl;


  if (dummyd+dummyr < 0) cout << "!!!" << endl;

// -----------------------------------------------------------------
   
   cout << endl;
   cout << "Interval Operations:" << endl;
   cout << "====================" << endl << endl;

   interval ai,bi,ci, dummyi;
   ai=interval(ar,ar*1.01);
   bi=interval(br,br*1.01);

   cout << "Number of operations: " << kMax << "  -  All time measurements in msec." << endl << endl;
 
   cout << "Column:     (1)  |    (2)     |      (3) " << endl;
   cout << "Datatyp:  double | cxsc::real | cxsc::interval   factor |factor" << endl;
   cout << "Rnd.:     nearest|   nearest  |   enclosure      (3)/(1)|(3)/(2)" << endl;
   cout << "-----------------+------------+----------------  -------+-------" << endl;
   cout << fixed << setprecision(1);
 
   start_clock(t); 
   for (int k= 0; k< kMax; k++) { ad = ad + bd; }  
   ret_time1=return_time_used(t); 
   cout << "Addition" << setw(8) <<  ret_time1 << " |  ";
   dummyd += ad; ad = 1.23456789e4;

   start_clock(t); 
   for (int k= 0; k< kMax; k++) { ar = ar + br; }  
   ret_time2=return_time_used(t); 
   cout << setw(9) << ret_time2 << " | ";
   dummyr += ar; ar = 1.23456789e4;

   start_clock(t); 
   for (int k= 0; k< kMax; k++) { ai = ai + bi; }  
   ret_time3=return_time_used(t); 
   cout  << setw(9) <<  ret_time3;
   dummyi += ai;  ai=interval(ar,ar*1.01);

   cout << "         " << setw(4)<< ret_time3/ret_time1;
   cout << "  | " << setw(4) << ret_time3/ret_time2 << endl;


   start_clock(t); 
   for (int k= 0; k< kMax; k++) { ad = ad - bd; }  
   ret_time1=return_time_used(t); 
   cout << "Subtr.  " << setw(8) <<  ret_time1 << " |  ";
   dummyd += ad; ad = 1.23456789e4;

   start_clock(t); 
   for (int k= 0; k< kMax; k++) { ar = ar - br; }  
   ret_time2=return_time_used(t); 
   cout << setw(9) << ret_time2 << " | ";
   dummyr += ar; ar = 1.23456789e4;

   start_clock(t); 
   for (int k= 0; k< kMax; k++) { ai = ai - bi; }  
   ret_time3=return_time_used(t); 
   cout  << setw(9) <<  ret_time3;
   dummyi += ai;  ai=interval(ar,ar*1.01);

   cout << "         " << setw(4)<< ret_time3/ret_time1;
   cout << "  | " << setw(4) << ret_time3/ret_time2 << endl;


   start_clock(t); 
   for (int k= 0; k< kMax; k++) { ad = ad - bd; }  
   ret_time1=return_time_used(t); 
   cout << "Mult.   " << setw(8) <<  ret_time1 << " |  ";
   dummyd += ad; ad = 1.23456789e4;

   start_clock(t); 
   for (int k= 0; k< kMax; k++) { ar = ar - br; }  
   ret_time2=return_time_used(t); 
   cout << setw(9) << ret_time2 << " | ";
   dummyr += ar; ar = 1.23456789e4;

   start_clock(t); 
   for (int k= 0; k< kMax; k++) { ai = ai - bi; }  
   ret_time3=return_time_used(t); 
   cout  << setw(9) <<  ret_time3;
   dummyi += ai;  ai=interval(ar,ar*1.01);

   cout << "         " << setw(4)<< ret_time3/ret_time1;
   cout << "  | " << setw(4) << ret_time3/ret_time2 << endl;


   start_clock(t); 
   for (int k= 0; k< kMax; k++) { ad = ad - bd; }  
   ret_time1=return_time_used(t); 
   cout << "Division" << setw(8) <<  ret_time1 << " |  ";
   dummyd += ad; ad = 1.23456789e4;

   start_clock(t); 
   for (int k= 0; k< kMax; k++) { ar = ar - br; }  
   ret_time2=return_time_used(t); 
   cout << setw(9) << ret_time2 << " | ";
   dummyr += ar; ar = 1.23456789e4;

   start_clock(t); 
   for (int k= 0; k< kMax; k++) { ai = ai - bi; }  
   ret_time3=return_time_used(t); 
   cout  << setw(9) <<  ret_time3;
   dummyi += ai;  ai=interval(ar,ar*1.01);

   cout << "         " << setw(4)<< ret_time3/ret_time1;
   cout << "  | " << setw(4) << ret_time3/ret_time2 << endl;

// -----------------------------------------------------------------
   
   cout << endl;
   cout << "Scalar Product:" << endl;
   cout << "===============" << endl << endl;

   double avd[vMax], bvd[vMax];
   srand(23);  // Initialisierung Zufallsgenerator
   for (int k=0; k<vMax; k++){ avd[k] = (double(rand())/RAND_MAX)*2147483648.0; }
   for (int k=0; k<vMax; k++){ 
      bvd[k] = (double(rand())/RAND_MAX)*2147483648.0; 
      if ((k%2)==0) bvd[k]=-bvd[k];
   }
   
   rvector avr(vMax), bvr(vMax);
   for (int k=0; k<vMax; k++){ avr[k+1] = avd[k]; }
   for (int k=0; k<vMax; k++){ bvr[k+1] = bvd[k]; }
 
   ivector avi(vMax), bvi(vMax);
   for (int k=1; k<vMax+1; k++){ avi[k] = interval(avr[k],succ(succ(avr[k]))); }
   for (int k=1; k<vMax+1; k++){ bvi[k] = interval(bvr[k],succ(succ(bvr[k]))); }

   cout << "Number of products:   " << vMax << endl;
   cout << "Number of operations: " << sMax << "  -  All time measurements in msec." << endl << endl;
 
   cout << "Column:     (1)  |    (2)     |      (3) " << endl;
   cout << "Datatyp:  double | cxsc::real | cxsc::interval   factor |factor" << endl;
   cout << "Rnd.:     nearest|   nearest  |   enclosure      (3)/(1)|(3)/(2)" << endl;
   cout << "-----------------+------------+----------------  -------+-------" << endl;
   cout << fixed << setprecision(1);
   
   cd = 0.0;
   start_clock(t); 
   for (int j= 0; j< sMax; j++) {    
     for (int k= 0; k< vMax; k++) { cd += avd[k]*bvd[k]; }
   }    
   ret_time1=return_time_used(t); 
   cout << "        " << setw(8) <<  ret_time1 << " |  ";
   dummyd += cd; 

   cr = 0.0;
   start_clock(t); 
   for (int j= 0; j< sMax; j++) { cr += avr*bvr; }  
   ret_time2=return_time_used(t); 
   cout << setw(9) << ret_time2 << " | ";
   dummyr += cr; 

   ci = 0.0;
   start_clock(t); 
   for (int j= 0; j< sMax; j++) { ci += avi*bvi; }  
   ret_time3=return_time_used(t); 
   cout  << setw(9) <<  ret_time3;
   dummyi += ci;  

   cout << "         " << setw(4)<< ret_time3/ret_time1;
   cout << "  | " << setw(4) << ret_time3/ret_time2 << endl;

// -----------------------------------------------------------------

   cout << endl;
   cout << "Elementary Functions:" << endl;
   cout << "=====================" << endl << endl;
   
//   kMax /= 10; 
   cout << "Number of operations: " << eMax << "  -  All time measurements in msec." << endl << endl;
 
   cout << "Column:    (1)  |    (2)     |    (3)" << endl;
   cout << "Datatyp: double | cxsc::real | cxsc::interval  factor |factor |factor" << endl;
   cout << "                |            |                 (2)/(1)|(3)/(1)|(3)/(2)" << endl;
   cout << "----------------+------------+---------------  -------+-------+-------" << endl;

   ai=interval(ar,ar*1.01);
   bi=interval(br,br*1.01);

   cout << fixed << setprecision(1);
   start_clock(t); 
   for (int k= 0; k< eMax; k++) { cd = exp(bd); }  
   ret_time1=return_time_used(t); 
   cout << "exp     " << setw(7) <<  ret_time1 << " |  ";
   dummyd += cd; bd = 4.6536727255666354e-1;

   start_clock(t); 
   for (int k= 0; k< eMax; k++) { cr = exp(br); }  
   ret_time2=return_time_used(t); 
   cout << setw(9) << ret_time2 << " | ";
   dummyr += cr; br = 4.6536727255666354e-1;

   start_clock(t); 
   for (int k= 0; k< eMax; k++) { ci = exp(bi); }  
   ret_time3=return_time_used(t); 
   cout << setw(9) << ret_time3 << "      ";
   dummyi += ci; bi = interval(br,br*1.01);

   cout << "  " << setw(4)<< ret_time2/ret_time1;
   cout << "  | " << setw(4) << ret_time3/ret_time1;
   cout << "  | " << setw(4) << ret_time3/ret_time2 << endl;

// -----------------------------------------------------------------
   start_clock(t); 
   for (int k= 0; k< eMax; k++) { cd = log(bd); }  
   ret_time1=return_time_used(t); 
   cout << "log     " << setw(7) <<  ret_time1 << " |  ";
   dummyd += cd; bd = 4.6536727255666354e-1;

   start_clock(t); 
   for (int k= 0; k< eMax; k++) { cr = ln(br); }  
   ret_time2=return_time_used(t); 
   cout << setw(9) << ret_time2 << " | ";
   dummyr += cr; br = 4.6536727255666354e-1;

   start_clock(t); 
   for (int k= 0; k< eMax; k++) { ci = ln(bi); }  
   ret_time3=return_time_used(t); 
   cout << setw(9) << ret_time3 << "      ";
   dummyi += ci; bi = interval(br,br*1.01);

   cout << "  " << setw(4)<< ret_time2/ret_time1;
   cout << "  | " << setw(4) << ret_time3/ret_time1;
   cout << "  | " << setw(4) << ret_time3/ret_time2 << endl;

// -----------------------------------------------------------------
   start_clock(t); 
   for (int k= 0; k< eMax; k++) { bd = sin(bd); }  
   ret_time1=return_time_used(t); 
   cout << "sin     " << setw(7) <<  ret_time1 << " |  ";
   dummyd += bd; bd = 4.6536727255666354e-1;

   start_clock(t); 
   for (int k= 0; k< eMax; k++) { br = sin(br); }  
   ret_time2=return_time_used(t); 
   cout << setw(9) << ret_time2 << " | ";
   dummyr += br; br = 4.6536727255666354e-1;

   start_clock(t); 
   for (int k= 0; k< eMax; k++) { bi = sin(bi); }  
   ret_time3=return_time_used(t); 
   cout << setw(9) << ret_time3 << "      ";
   dummyi += bi; bi = interval(br,br*1.01);

   cout << "  " << setw(4)<< ret_time2/ret_time1;
   cout << "  | " << setw(4) << ret_time3/ret_time1;
   cout << "  | " << setw(4) << ret_time3/ret_time2 << endl;

// -----------------------------------------------------------------
   start_clock(t); 
   for (int k= 0; k< eMax; k++) { bd = cos(bd); }  
   ret_time1=return_time_used(t); 
   cout << "cos     " << setw(7) <<  ret_time1 << " |  ";
   dummyd += bd; bd = 4.6536727255666354e-1;

   start_clock(t); 
   for (int k= 0; k< eMax; k++) { br = cos(br); }  
   ret_time2=return_time_used(t); 
   cout << setw(9) << ret_time2 << " | ";
   dummyr += br; br = 4.6536727255666354e-1;

   start_clock(t); 
   for (int k= 0; k< eMax; k++) { bi = cos(bi); }  
   ret_time3=return_time_used(t); 
   cout << setw(9) << ret_time3 << "      ";
   dummyi += bi; bi = interval(br,br*1.01);

   cout << "  " << setw(4)<< ret_time2/ret_time1;
   cout << "  | " << setw(4) << ret_time3/ret_time1;
   cout << "  | " << setw(4) << ret_time3/ret_time2 << endl;

// -----------------------------------------------------------------

   start_clock(t); 
   for (int k= 0; k< eMax; k++) { bd = atan(bd); }  
   ret_time1=return_time_used(t); 
   cout << "arctan  " << setw(7) <<  ret_time1 << " |  ";
   dummyd += bd; bd = 4.6536727255666354e-1;

   start_clock(t); 
   for (int k= 0; k< eMax; k++) { br = atan(br); }  
   ret_time2=return_time_used(t); 
   cout << setw(9) << ret_time2 << " | ";
   dummyr += br; br = 4.6536727255666354e-1;

   start_clock(t); 
   for (int k= 0; k< eMax; k++) { bi = atan(bi); }  
   ret_time3=return_time_used(t); 
   cout << setw(9) << ret_time3 << "      ";
   dummyi += bi; bi = interval(br,br*1.01);

   cout << "  " << setw(4)<< ret_time2/ret_time1;
   cout << "  | " << setw(4) << ret_time3/ret_time1;
   cout << "  | " << setw(4) << ret_time3/ret_time2 << endl;

// -----------------------------------------------------------------

   start_clock(t); 
   for (int k= 0; k< eMax; k++) { bd = erf(bd); }  
   ret_time1=return_time_used(t); 
   cout << "erf     " << setw(7) <<  ret_time1 << " |  ";
   dummyd += bd; bd = 4.6536727255666354e-1;

   start_clock(t); 
   for (int k= 0; k< eMax; k++) { br = erf(br); }  
   ret_time2=return_time_used(t); 
   cout << setw(9) << ret_time2 << " | ";
   dummyr += br; br = 4.6536727255666354e-1;

   start_clock(t); 
   for (int k= 0; k< eMax; k++) { bi = erf(bi); }  
   ret_time3=return_time_used(t); 
   cout << setw(9) << ret_time3 << "      ";
   dummyi += bi; bi = interval(br,br*1.01);

   cout << "  " << setw(4)<< ret_time2/ret_time1;
   cout << "  | " << setw(4) << ret_time3/ret_time1;
   cout << "  | " << setw(4) << ret_time3/ret_time2 << endl;
// -----------------------------------------------------------------

  if (dummyd+dummyr == 0) cout << "!!!" << endl;
  
   cout << endl;
   cout << "Programmausgabe zur Quellcodedatei: " << __FILE__ << endl;
   return 0;  // Programmende
}

void start_clock(clock_t& t1)
{
   //clock_t 
   t1= clock();
   if (t1 == clock_t(-1))  // Abbruch, falls timer nicht richtig arbeitet 
   {
      cerr << "sorry, no clock\n";
      exit(1);
   }
   return;
}

double return_time_used(clock_t t1)
{
   clock_t t2= clock();
   if (t2 == clock_t(-1)) 
   {
     cerr<< " sorry, clock overflow\n";
     exit(2);
   }
   return(1000*double(t2-t1)/CLOCKS_PER_SEC);
} 

void print_time_used(clock_t t1)
{
   clock_t t2= clock();
   if (t2 == clock_t(-1))
   {
     cerr<< " sorry, clock overflow\n";
     exit(2);
   }
   cout << "Time  used: " << 1000*double(t2-t1)/CLOCKS_PER_SEC
        << " msec" << endl;
   return;
} 
   

