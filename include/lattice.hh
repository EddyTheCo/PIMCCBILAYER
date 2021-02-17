#include <cstdlib>
#include<stack>
#include <fstream>
#include<iostream>
#ifdef USEROOT
#include<TFile.h>
#include<TH1.h>
#include<TH3.h>

#include<TVector.h>
#endif

#include <complex>
#ifndef LATTICE_HH
#define LATTICE_HH

#include"block.hh"


using namespace std;


extern  bool isGrandCanonical;
extern const bool restart;
extern size_t Warmup;

extern const size_t NPartiIni,SAMPLING;
class lattice
{
public:
    lattice();

    void setup()const;
    void move()const;
#ifdef WARMUP
    void Warm() const;
#endif
     void PrintConfiguration (
        #ifdef USEROOT
        const size_t step
        #else
        void
        #endif
             ) const;


void initialize_histos(void);
inline double getOpenRatio(void)const
{
    return Site::NOpen*1./Site::NOpenP;
}
inline double getCloseRatio(void)const
{
    return Site::NClose*1./Site::NCloseP;
}
inline double getMoveRatio(void)const
{
    return Site::NMove*1./Site::NMoveP;
}
inline double getSwapRatio(void)const
{
    return Site::NSwap*1./Site::NSwapP;
}
inline double getWiggleRatio(void)const
{
    return Site::NWiggle*1./Site::NWiggleP;
}
inline double getShiftRatio(void)const
{
    return Site::NShift*1./Site::NShiftP;
}
inline double getInsertRatio(void)const
{
    return Site::NInsert*1./Site::NInsertP;
}
inline double getRemoveRatio(void)const
{
    return Site::NRemo*1./Site::NRemoP;
}
  static array<vector<Site>,100000>*  const grid;

    static ofstream  thesweep,theratios;



    static const size_t NRep,NSweeps;
#ifdef USEROOT
    static TFile *RootFile;
    static TH1* hpos1,*hpos2;
    static TH1D * PCFUp,*PCFDown,*PCFMix;
    static TH2D *Greens;

static TVectorD * v;
#endif

};

#endif // LATTICE_HH
