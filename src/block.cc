#include"block.hh"
#include<iostream>
#include<fstream>
#include"lattice.hh"
#include<vector>
#include"constants.hh"
#ifdef USEROOT
#include<TH2D.h>
#endif
using namespace std;

#include"site.hh"




block::block(array<vector<Site>,10000>* particles, const size_t &NTimeSlices, const size_t &NSweeps, const bool &realB
             #ifdef USEROOT
             , TH2D * const Greens
             #endif
             )
{


double TSumOfdisplacement=0,TSumOfPotential=0,TNumberOfParticles=0;
size_t TWormlenght=0,step=0,measureCounter=0,measureCounter1=0;
double TWinding=0;
Site* const start=&(particles->at(0).at(0));



while(step<NSweeps)
{



//start->printLattice();
//cout<<start->NParti_<<" "<<step<<endl;


         if(start->ThereIsAWorm)
        {

             if(realB)
             {
                measureCounter1++;
                TWormlenght+=start->NInactiveLinks();
#ifdef USEROOT
                Greens->Fill(sqrt((start->Rbead->pos-start->Lbead->pos).norm()),abs(1.*start->Rbead->TimeSliceOnBead-1.*start->Lbead->TimeSliceOnBead));
#endif
             }
            switch ((!isGrandCanonical)?giveRanI(2):giveRanI(3)) {
            case 0:
            {
               //cout<<"closing worm"<<endl;
                    start->NCloseP++;

                    if(start->Lbead->CloseWorm(0))
                    {
                        step++;
                    }
                break;
            }
            case 1:
            {
               //cout<<"MoveWorm"<<endl;
                    start->MoveWorm();
                     break;
            }
            case 2:
            {
              // cout<<"swap"<<endl;
                start->NSwapP++;
               if(start->NParti_>1)
               start->PrepareSwap();

               break;
            }
            case 3:
            {

                //cout<<"removeWorm"<<endl;
                start->removeWorm();
                break;
            }



            }

        }
        else
        {
            if(realB)
            {
                TSumOfdisplacement+=start->TEnergy;
                TSumOfPotential+=start->TPotential;
                TNumberOfParticles+=start->NParti_;
                TWinding+=start->TWinding.norm();
                measureCounter++;
            }

             switch ((isGrandCanonical)?giveRanI(2):giveRanI(1)) {
             case 0:
             {

               if(start->NParti_)
               {
               //    cout<<"OpenWorm"<<endl;
                   const size_t posiTimes=giveRanI(NTimeSlices-1) ;
                   const size_t posiParti=giveRanI(particles->at(posiTimes).size()-1);
                   const size_t var2=  giveRanI(MBar-2);
                   Site* const Ranbead=&(particles->at(posiTimes).at(posiParti));
                   start->ThereIsAWorm= Ranbead->OpenWorm(var2,var2+1,0,Ranbead->pos);
               }

                break;
             }
             case 1:
            {
                 if(start->NParti_)
                 {

             //      cout<<"wiggle"<<endl;

                     const size_t posiTimes=giveRanI(NTimeSlices-1) ; //Choose a random time slice
                     const size_t posiParti=giveRanI(particles->at(posiTimes).size()-1); //Choose the particle

                        start->Lbead=&(particles->at(posiTimes).at(posiParti)); //LBEAD is proposed (but dosent mean theres is a worm)
                        const size_t var2= giveRanI(MBar-3)+1;

                        start->Rbead=start->Lbead->searchBead(true,var2);
                        start->Lbead->oldpos=start->Lbead->pos;

                        start->Lbead->Wiggle(0);

                        start->Lbead=nullptr;
                        start->Rbead=nullptr;

                }
                break;
            }
              case 2:
             {
                 //cout<<"insertworminclose "<<Constants::giveRanD(1)<<" "<<Constants::giveRanDNormal(0,1)<<endl;
                 start->insertWorm();
                 break;
             }
             }


        }

    }


if(realB)
{
        SumofDisplacement=TSumOfdisplacement/measureCounter;
        SumOfPotential=TSumOfPotential/measureCounter;
        NumberOfParticles=TNumberOfParticles/measureCounter;
        Wormlenght=1.*TWormlenght/measureCounter1;
        SumofWinding=TWinding/measureCounter;
}

}




block::~block()
{

}
