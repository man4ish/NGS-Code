 #include <iostream>
 #include <fstream>
 #include <vector>
 using namespace std;
 vector<int> v1 ;
 void Formatdata(string asd,int width,bool firstcolflag)
 {
      int len = asd.size();
      int pad = width -len;
      if(pad > 0)
      {
          for(int i  =0; i < pad; i++)
          {
              if(firstcolflag == 1)
                 cout<< '0';
              else 
                 cout<< ' ';
          }          
          for(int ii  =0; ii < asd.size(); ii++)
          {
              cout << asd[ii];
          }
      }
      else 
      {
          for(int ii  =0; ii < width; ii++)
          {
              cout<< asd[ii];
          } 
      }
 }

 void init_vector(int numfields)
 {
    for (int i =0; i < numfields; i++)
    v1.push_back(0);                  
 } 

 int main(int argc ,char* argv[])
 {
    ifstream in(argv[1]);
    char line[1000];
    bool fieldflag = 0;
    while(in)
    {
       in.getline(line, 1000);
       vector <const char*> vec;
       if(in)
       {
          char* pch;
          pch = strtok(line,"\t");
          while(pch != NULL)
          {
                 vec.push_back(pch);
                 pch = strtok(NULL,"\t");
          }
          if(fieldflag == 0)
          {
            init_vector(vec.size()); 
            fieldflag = 1;
          }
          vector <int> vecsize;
          for (int i = 0; i < vec.size(); i++ )
          {
               vecsize.push_back(strlen(vec[i]));
          }
          for (int s = 0; s < vecsize.size(); s++)
          {
               if(v1[s] <  vecsize[s])
                  v1[s] = vecsize[s];      
          }   
       }
    }
    in.close();     

   ifstream featurefile(argv[1]);
   char featureline[1000];
   while(featurefile)
   {
       featurefile.getline(featureline, 1000);
       vector <const char*> vecbest;
       if(featurefile)
       {
          char* ppch;
          ppch = strtok(featureline,"\t");
          while(ppch != NULL)
          {
                 vecbest.push_back(ppch);
                 ppch = strtok(NULL,"\t");
          }
          Formatdata(vecbest[0],v1[0],1);
          cout << "\t";
          for(int ss =1; ss < vecbest.size() -1; ss++)
          {
              Formatdata(vecbest[ss],v1[ss],0);
              cout << "\t";
          }
          Formatdata(vecbest[vecbest.size() -1],v1[vecbest.size() -1],0);
          cout << "\n";
       }
    }
    featurefile.close();
    return 0;
 }
