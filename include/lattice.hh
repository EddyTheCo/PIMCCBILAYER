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


extern const bool isGrandCanonical,restart;

extern const size_t NPartiIni,SAMPLING;
class lattice
{
public:
    lattice();

    void setup()const;
    void move()const;

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
    return grid->at(0).at(0).NOpen*1./grid->at(0).at(0).NOpenP;
}
inline double getCloseRatio(void)const
{
    return grid->at(0).at(0).NClose*1./grid->at(0).at(0).NCloseP;
}
inline double getMoveRatio(void)const
{
    return grid->at(0).at(0).NMove*1./grid->at(0).at(0).NMoveP;
}
inline double getSwapRatio(void)const
{
    return grid->at(0).at(0).NSwap*1./grid->at(0).at(0).NSwapP;
}
inline double getInsertRatio(void)const
{
    return grid->at(0).at(0).NInsert*1./grid->at(0).at(0).NInsertP;
}
inline double getRemoveRatio(void)const
{
    return grid->at(0).at(0).NRemo*1./grid->at(0).at(0).NRemoP;
}
  static array<vector<Site>,10000>*  const grid;

    static ofstream  thesweep,theratios;



    static const size_t NRep,NSweeps;
#ifdef USEROOT
    static TFile *RootFile;
    static TH1* hpos;
    static TH1D * PCFUp,*PCFDown,*PCFMix;
    static TH2D *Greens;

static TVectorD * v;
#endif

};

#endif // LATTICE_HH
