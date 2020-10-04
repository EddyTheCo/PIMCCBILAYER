#include<iostream>
#include<fstream>

#include<stack>
#include<vector>
#ifdef USEROOT
#include<TH2D.h>
#include<TH1.h>
#endif
#ifndef BLOCK_H
#define BLOCK_H


using namespace std;
#include"site.hh"
class block
{
    public:
    block(array<vector<Site>,10000>* particles, const size_t &NTimeSlices,  const size_t &NSweeps, const bool &realB
      #ifdef USEROOT
      ,  TH2D * const Greens
      #endif
          );

    ~block();

    double SumofDisplacement, NumberOfParticles,SumOfPotential;
    double Wormlenght;
    double SumofWinding;
    inline double getKineticEnergy(void)const{
        return (NumberOfParticles*d*0.5/tao-SumofDisplacement/(4*landa*tao*tao*NTimeSlices));
        }
    inline double getSuperfluidDensity(void)const{
        return SumofWinding/(d*2*landa*beta*NumberOfParticles);
        }


private:

};



#endif // BLOCK_H
