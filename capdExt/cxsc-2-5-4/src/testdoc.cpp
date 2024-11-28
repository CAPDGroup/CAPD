/*
**  CXSC is a C++ library for eXtended Scientific Computing (V 2.5.4)
**
**  Copyright (C) 1990-2000 Institut fuer Angewandte Mathematik,
**                          Universitaet Karlsruhe, Germany
**            (C) 2000-2014 Wiss. Rechnen/Softwaretechnologie
**                          Universitaet Wuppertal, Germany   
**
**  This library is free software; you can redistribute it and/or
**  modify it under the terms of the GNU Library General Public
**  License as published by the Free Software Foundation; either
**  version 2 of the License, or (at your option) any later version.
**
**  This library is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
**  Library General Public License for more details.
**
**  You should have received a copy of the GNU Library General Public
**  License along with this library; if not, write to the Free
**  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/* CVS $Id: testdoc.cpp,v 1.23 2014/01/30 17:23:49 cxsc Exp $ */

#include <string>
#include <map>
#include <fstream>

namespace cxsc {

using namespace std;

string parsearg(istream &in)
{
   string arg;
   int depth=0;
   char c;
   in >> c;
   do {
      if(c=='{')
      {
         if(depth)
            arg+=c;
         depth++;
      } else if(c=='}')
      {
         depth--;
         if(depth)
            arg+=c;
      } else
         arg+=c;
         
      if(depth>0)
         in >> c;
   } while(depth>0 && !in.eof() && !in.fail());
   return arg;
}

int main(int argc,char *argv[])
{
   map<string,int> tested;
   ifstream in;
   if(argc>=2)
      in.open(argv[1]);

   if(argc<2 || in.fail())
   {
      cerr << argv[0] << ": <test-output> <fail-output> <test-tex infile> <test-tex outfile> [<test-tex infile> <test-tex outfile>[<test-tex...]]" << endl; 
      return 1;
   }
   while(!in.eof())
   {
      string instring;
      getline(in,instring);
      instring=" "+instring;
      if(int(instring.find(":"))>0)
         tested[instring]=1;
   }
   in.close();
   
   typedef map<string,int>::const_iterator CI;
   
   int i=0;
   
   CI p;
   
   for(p=tested.begin();p!=tested.end();++p) i++;

   cout << "read " << i << " lines..." << endl;

   ofstream out;
   
   for(i=3;i+1<argc;i+=2)
   {
//      cout << "in:  " << argv[i] << endl;  
      in.open(argv[i]);
//      cout << "out: " << argv[i+1] << endl;
      out.open(argv[i+1]);
      
      in.unsetf(ios::skipws);
      in.setf(ios::binary);
      
      char c;
      in >> c;
      while(!in.eof() && !in.fail())
      {
         if(c=='\\')
         {
            in >> c;
            if(c=='t')
            {
               in >> c;
               if(c=='s')
               {
                  in >> c;
                  if(c=='t')
                  {
                     string first(parsearg(in));
                     string second(parsearg(in));
                     string dummy(parsearg(in));
                     int missing=0;
                     int failed=0;
                     int ok=0;
                     
                     
                     while(second.length()>0)
                     {
                        int pos=second.find(";");
                        string tfunc(second);
                        if(pos>0)
                        {
                           tfunc=second.substr(0,pos);
                           second=second.substr(pos+1);
                        } else
                           second="";
                           
                        tfunc=" "+tfunc;
                        cout << argv[i] << ": Checking " << tfunc << "         \r";
                        cout.flush();
                             
                        int found=0;
                           
                        for(p=tested.begin();p!=tested.end();++p)
                           if(p->first!="" && int(p->first.find(tfunc))>0)
                           {
                              found=1;
                              break;
                           }
                        
                        
                        
                        if(found)
                        {
                           if(p->second)
                           {
                              if(int(p->first.find("failed"))>=0)
                                 failed=1;
                              else
                                 ok=1;
                              tested[string(p->first)]++;
                           } else
                           {
                              missing=1;
                           }
                        } else
                        {
                           missing=1;
                           tested[tfunc]=0;
                        }
                     }
                     string status="\\color{green}";
                     if(ok && (missing || failed))
                        status="$\\frac12$";
                     if(missing && failed)
                        status+="\\color{magenta}";
                     else if(missing)
                        status+="\\color{red}";
                     else if(failed)
                        status+="\\{cyan}"; 
                        
                     out << "\\tst{" << first << "}{}{" << status << "}";
                  } else
                     out << "\\ts" << c;
               } else
                  out << "\\t" << c;
            } else
               out << '\\' << c;
         } else
            out << c;
         in >> c;
      };
      out.close();      
      in.close();
   }
   out.open(argv[2]);

   if(!out.fail())
   {
      int none=1;
      out << "\\subsection{Missing tests}"<< endl;
      out << "\\begin{verbatim}" << endl;
      long count(0);
      long all(0);
      for(p=tested.begin();p!=tested.end();++p)
      {
         if(p->second==0)
         {
            none=0;
            out << p->first << endl;
            count++;
         }
         all++;
      }
      out << "\\end{verbatim}" << endl;
      if(none)
         out << "none" << endl;
      cout << endl;
      cout << count << " missing tests ("<<100*count/all<<"%)" << endl;
      
      out << "\\subsection{Unused tests}" << endl;
      out << "\\begin{verbatim}" << endl;
      none=1;
      count=0;
      for(p=tested.begin();p!=tested.end();++p)
      {
         if(p->second==1)
         {
            none=0;
            count++;
            out << p->first << endl;
         }
      }
      out << "\\end{verbatim}" << endl;
      if(none) 
         out << "none" << endl;
      out.close();
      cout << count << " unused tests (" <<100*count/all<<"%)" << endl;
      cout << all << " checked tests" << endl;
   }
   cout << endl;
   return 0;
}

} // namespace cxsc

