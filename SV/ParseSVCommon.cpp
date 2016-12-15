//GenomeID\tUnique\tChromosome\tOffset\tVariant Type\tNucleotide_Variant - single base\tSupporting Reads (forward strand)\tSupporting Reads (reverse strand)\tSupporting Reads (inverse strand)\tValue: Total_Read_Density or CN\tAvg_QScore\tdbSNP ID\tGene_Name\tGene ID\tmiRNA ID\tPromoter\tExon\tZygosity\tLocal_Copy_Number\tKnown_CNVRegion\tTransversions(Tv)_Transitions(Ts)\tSynonymous(S)_Non-Synonymous(NS)\tProtein_Variant\tMissense\tNonsense\tSNV Confident Score\tSV - Length or CN - end of offset\tSV Link ID\tSV Edge Support

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
using namespace std;
struct Key
{
    int gi;
    size_t offset;
    string allele;
    Key(int _gi, size_t  _offset, string _allele)
    {
      gi = _gi ;
      offset = _offset;
      allele = _allele;
    }
};

struct ltstr
{
    inline bool operator()(const Key* _lhs, const Key* _rhs) const
    {
        if(_lhs->gi != _rhs->gi)
         return ( _lhs->gi < _rhs->gi);
        else if(_lhs->offset != _rhs->offset)
           return  ( _lhs->offset < _rhs->offset );
          return (strcmp((_lhs->allele).c_str(),(_rhs->allele).c_str())< 0);
    }
};

typedef map<Key*,string,ltstr> m;

int main(int argc, char* argv[])
{
 cout << "GenomeID\tUnique\tChromosome\tOffset\tVariant Type\tNucleotide_Variant - single base\tSupporting Reads (forward strand)\tSupporting Reads (reverse strand)\tSupporting Reads (inverse strand)\tValue: Total_Read_Density or CN\tAvg_QScore\tdbSNP ID\tGene_Name\tmiRNA ID\tPromoter\tExon\tZygosity\tLocal_Copy_Number\tKnown_CNVRegion\tTransversions(Tv)_Transitions(Ts)\tSynonymous(S)_Non-Synonymous(NS)\tProtein_Variant\tMissense\tNonsense\tSNV Confident Score\tSV - Length or CN - end of offset\tSV Link ID\tSV Edge Support"<<endl; 
ifstream in(argv[1]);
char line[10000];
int count =0;
while(in)
{
  in.getline(line,10000);
  vector <const char*> vec; 
  if(in )
  {
    char* pch;
    pch = strtok(line,"\t");
    while(pch != NULL)
    {
       vec.push_back(pch);  
       pch = strtok(NULL,"\t");   
    }
    count++;
    if(vec[0][0] != '#')
    {
     string h = vec[16];
     int i = h.find(":");    
     string side1 = h.substr(0,i);
     string side2 = h.substr(i+1,h.size());
     string ID = argv[2];
     int i2 = ID.find("_");
     string side1genome = ID.substr(0,i2);
     string side2genome = ID.substr(i2+1,ID.size());
     cout<< side1genome<<"\t"<< argv[3]<<"\t"<< vec[0]<<"\t"<<vec[1]<<"\t"<<argv[4]<<"1"<<"\t"<<"-"<<"\t"<<vec[7]<<"\t"<< vec[8]<<"\t"<<atoi(vec[10])+atoi(vec[9])<<"\t"<<"-"<<"\t"<<side1<<"\t-\t"<<vec[12]<<"\t"<<vec[23]<<"\t"<<vec[25]<<"\t"<<vec[14]<<"\t"<<vec[19]<<"\t"<<"-"<<"\t"<<"-"<<"\t"<<"-"<<"\t"<<"-"<<"\t"<<"-"<<"\t"<<"-"<<"\t"<<"-"<<"\t"<<vec[20]<<"\t"<<atoi(vec[1])+atoi(vec[2])<<"\t"<<argv[2]<<count<<"_"<<side1genome<<"\t-\t"<<endl;
     cout<< side1genome<<"\t"<< argv[3]<<"\t"<< vec[4]<<"\t"<<vec[5]<<"\t"<<argv[4]<<"2"<<"\t"<<"-"<<"\t"<<vec[7]<<"\t"<< vec[8]<<"\t"<<atoi(vec[10])+atoi(vec[9])<<"\t"<<"-"<<"\t"<<side1<<"\t-\t"<<vec[13]<<"\t"<<vec[24]<<"\t"<<vec[26]<<"\t"<<vec[15]<<"\t"<<vec[19]<<"\t"<<"-"<<"\t"<<"-"<<"\t"<<"-"<<"\t"<<"-"<<"\t"<<"-"<<"\t"<<"-"<<"\t"<<"-"<<"\t"<<vec[20]<<"\t"<<atoi(vec[5])+atoi(vec[6])<<"\t"<<argv[2]<<count<<"_"<<side1genome<<"\t-\t"<<endl;

cout<< side2genome<<"\t"<< argv[3]<<"\t"<< vec[0]<<"\t"<<vec[1]<<"\t"<<argv[4]<<"1"<<"\t"<<"-"<<"\t"<<vec[7]<<"\t"<< vec[8]<<"\t"<<atoi(vec[10])+atoi(vec[9])<<"\t"<<"-"<<"\t"<<side2<<"\t-\t"<<vec[12]<<"\t"<<vec[23]<<"\t"<<vec[25]<<"\t"<<vec[14]<<"\t"<<vec[19]<<"\t"<<"-"<<"\t"<<"-"<<"\t"<<"-"<<"\t"<<"-"<<"\t"<<"-"<<"\t"<<"-"<<"\t"<<"-"<<"\t"<<vec[20]<<"\t"<<atoi(vec[1])+atoi(vec[2])<<"\t"<<argv[2]<<count<<"_"<<side2genome<<"\t-\t"<<endl;
     cout<< side2genome<<"\t"<< argv[3]<<"\t"<< vec[4]<<"\t"<<vec[5]<<"\t"<<argv[4]<<"2"<<"\t"<<"-"<<"\t"<<vec[7]<<"\t"<< vec[8]<<"\t"<<atoi(vec[10])+atoi(vec[9])<<"\t"<<"-"<<"\t"<<side2<<"\t-\t"<<vec[13]<<"\t"<<vec[24]<<"\t"<<vec[26]<<"\t"<<vec[15]<<"\t"<<vec[19]<<"\t"<<"-"<<"\t"<<"-"<<"\t"<<"-"<<"\t"<<"-"<<"\t"<<"-"<<"\t"<<"-"<<"\t"<<"-"<<"\t"<<vec[20]<<"\t"<<atoi(vec[5])+atoi(vec[6])<<"\t"<<argv[2]<<count<<"_"<<side2genome<<"\t-\t"<<endl;
}
}

}
  in.close();
  return 0;
}


