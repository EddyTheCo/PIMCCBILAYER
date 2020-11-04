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


extern const int Warmup;
class block
{
    public:
    block(array<vector<Site>,10000>* particles, const size_t &NTimeSlices,  const size_t &NSweeps
      #ifdef USEROOT
      ,  TH2D * const Greens
      #endif
          );

    ~block();

    double SumofDisplacement, NumberOfParticles,NumberOfParticlesUp,SumOfPotential;
    double Wormlenght;
    double SumofWindingUp,SumofWindingDown;
    inline double getKineticEnergy(void)const{
        return (NumberOfParticles*d*0.5/tao-SumofDisplacement/(4*landa*tao*tao*NTimeSlices));
        }
    inline double getSuperfluidDensityUp(void)const{
        return SumofWindingUp/(d*2*landa*beta*NumberOfParticlesUp);
        }
    inline double getSuperfluidDensityDown(void)const{
        return SumofWindingDown/(d*2*landa*beta*(NumberOfParticles-NumberOfParticlesUp));
        }



private:

};



#endif // BLOCK_H
