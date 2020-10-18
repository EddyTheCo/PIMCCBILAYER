#include<random>
#include<fstream>
#include<stack>
#include<iostream>
using namespace std;

#ifndef CONSTANTS_H
#define CONSTANTS_H 1


namespace Constants
{

    // forward declarations only
extern mt19937_64 mt;

inline double giveRanD(const double &A){
  uniform_real_distribution<double> dist(0,A);
 return dist(mt);
}
inline size_t giveRanI(const size_t &A){
 uniform_int_distribution<size_t> dist(0,A);
 return dist(mt);
}
inline size_t giveT(const size_t &A){
 uniform_int_distribution<size_t> dist(0,A);
 auto var=dist(mt);

 return (var==A)?A+1:var;
}
inline double giveRanDNormal(const double &A, const double &B){
 normal_distribution<double> dist(A,B);
    return dist(mt);

}

     double calculateStd(stack<double> vect, double ave);
     void saveRandom(void);
     void readRandom(void);
extern uniform_real_distribution<double> distd;
    extern const double pi;
    extern const double avogadro;
    extern const double my_gravity;
}

template <class myType>
 myType ReadFromInput(const int line, const string in="input")
{
    fstream  myfilein;
      myfilein.open (in);
      if(!myfilein.is_open())
      {
          cout<<"Error; The file "<<in<<" couldn't be opened"<<endl;
          throw std::exception();
      }
      string linestr;
      myType A;
     for(int i=1;i<line;i++)
     {
       getline(myfilein,linestr);
     }
     myfilein>>A>>linestr;
     cout<<linestr<<"="<<A<<endl;
     myfilein.close();
    return A;
}

#endif
