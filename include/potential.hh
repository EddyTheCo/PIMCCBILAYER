#include<fstream>
#include<iostream>

using namespace std;
#ifndef POTENTIAL_HH
#define POTENTIAL_HH

class Site;
#include"position.hh"


class potential
{
public:
    static ifstream PotConf;
    potential(){}
      void U(const Site* const bead,double& dU,position& graddU ) const;


     void ExternPot  (const Site* const ,double& dU,position& graddU )const;
     void PairInteraction  (const Site* const bead, double& dU, position& graddU )const;
     void LiebLini(double & dU,position& graddU,const Site* const bead,const Site* const ptr)const;
     void Dipolar(double & dU,position& graddU,const Site* const bead,const Site* const ptr)const;
    void Softcore(double & dU,position& graddU,const Site* const bead,const Site* const ptr)const;
    void BoninPot(double & dU,position& graddU,const Site* const bead,const Site* const ptr)const;

    static const bool harmonic,freeExt,infinteWell;
    static const bool LiebLin, Dipola,Softcor,Bonin,freeInt;
    static const double CLieb,CDip,V_0Softcore,SigmaBonin;
    static const position k;
};

#endif // POTENTIAL_HH
