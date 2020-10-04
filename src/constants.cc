#include<random>
#include<iostream>
#include<fstream>
#include<string>
#include<stack>
#include<string>
#include<constants.hh>
#include"position.hh"
using namespace std;

   // actual global variables


 //random_device rd;
 mt19937_64 Constants::mt(47814181);



     const double Constants::pi(3.141592654);
     const double Constants::avogadro(6.0221413e23);
     void Constants::readRandom(void)
    {
        ifstream fin(".seed.dat");
        fin>>mt;
        fin.close();

    }
     void Constants::saveRandom(void)
    {

        ofstream fout(".seed.dat");
        fout<<mt;
        fout.close();
    }

    inline double Constants::calculateStd(stack<double> vect, double ave)
    {
        double sum=0;
        size_t k=vect.size();
        if(k==0)
            return 0;
        while(!vect.empty())
        {
            sum+=(vect.top()-ave)*(vect.top()-ave);
            vect.pop();
        }
        return sqrt(sum/k);

    }



extern template  int ReadFromInput<int>( const int,const string in="input");
extern template  double ReadFromInput<double>(const int,const string in="input");
extern template  position ReadFromInput<position>(const int,const string in="input");
extern template  size_t ReadFromInput<size_t>(const int,const string in="input");
extern template  string ReadFromInput<string>(const int,const string in="input");


