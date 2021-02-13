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




block::block(array<vector<Site>,100000>* particles, const size_t &NTimeSlices, const size_t &NSweeps
             #ifdef USEROOT
             , TH2D * const Greens
             #endif
             )
{


double TSumOfdisplacement=0,TSumOfPotential=0,TNumberOfParticles=0,TNumberOfParticlesUp=0;
size_t TWormlenght=0,step=0,measureCounter=0,measureCounter1=0;
double TWindingUp=0,TWindingDown=0;
Site* const start=&(particles->at(0).at(0));

string svar="";


while(step<NSweeps)
{

//start->printLattice(svar);

         if(start->ThereIsAWorm)
        {

#ifdef WARMUP
             if(!Warmup)
#endif
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
              // svar="closing worm";

                    start->NCloseP++;


                    if(!(start->cantClose(MBar)))
                    {

                            if(start->Lbead->CloseWorm(0))
                            {
                                step++;

                            }



                    }

                break;
            }
            case 1:
            {
              //cout<<"MoveWorm"<<endl;
              //svar="MoveWorm";
                    start->NMoveP++;
                    start->MoveWorm();
                     break;
            }
            case 2:
            {
              //cout<<"swap"<<endl;
             // svar="swap";

                start->NSwapP++;
               if(start->NParti_>1)
               start->PrepareSwap();

               break;
            }
            case 3:
            {

                //cout<<"removeWorm"<<endl;
               // svar="removeWorm";
                start->NRemoP++;
                if(start->removeWorm())
                {
                    step++;

                }
                break;
            }



            }

        }
        else
        {
            #ifdef WARMUP
            if(!Warmup)
            #endif
            {
                TSumOfdisplacement+=start->TEnergy;
                TSumOfPotential+=start->TPotential;
                TNumberOfParticles+=start->NParti_;
                TNumberOfParticlesUp+=start->Nparti_UpxNT/NTimeSlices;
                TWindingUp+=start->TWindingUp.normxy();

                TWindingDown+=start->TWindingDown.normxy();
                measureCounter++;
            }

             switch ((isGrandCanonical)?giveRanI(2):giveRanI(1)) {
             case 0:
             {

               if(start->NParti_)
               {
                 //cout<<"OpenWorm"<<endl;
                   // svar="OpenWorm";
                   const size_t posiTimes=giveRanI(NTimeSlices-1) ;
                   const size_t posiParti=giveRanI(particles->at(posiTimes).size()-1);
                   const size_t var2=  giveRanI(MBar-2);
                   Site* const Ranbead=&(particles->at(posiTimes).at(posiParti));
                   start->NOpenP++;

                           start->ThereIsAWorm= Ranbead->OpenWorm(var2,var2+1,0,Ranbead->pos);

               }

                break;
             }
             case 1:
            {
                 if(start->NParti_)
                 {

                  //cout<<"wiggle"<<endl;
        //svar="wiggle";

                     const size_t posiTimes=giveRanI(NTimeSlices-1) ; //Choose a random time slice
                     const size_t posiParti=giveRanI(start->NParti_-1); //Choose the particle

                        start->Lbead=&(particles->at(posiTimes).at(posiParti)); //LBEAD is proposed (but dosent mean theres is a worm)
                        const size_t var2= giveRanI(MBar-3)+1;

                        start->Rbead=start->Lbead->searchBead(true,var2);
                        start->Lbead->oldpos=start->Lbead->pos;
                        start->NWiggleP++;
                        start->Lbead->Wiggle(0);

                        start->Lbead=nullptr;
                        start->Rbead=nullptr;

                }
                break;
            }

              case 3:
             {
                 //cout<<"insertworminclose "<<endl;
            //     svar="insertworminclose";
                 start->NInsertP++;
                 start->insertWorm();
                 break;
             }
             case 5:
            {
                 start->NShiftP++;
                 if(start->NParti_)
                 {

                 //cout<<"ShiftParticle"<<endl;
                // svar="ShiftParticle";



                     const size_t posiParti=giveRanI(start->NParti_-1); //Choose the particle

                        start->Lbead=&(particles->at(0).at(posiParti)); //LBEAD is proposed (but dosent mean theres is a worm)


                        vector<double> varVec;

                        varVec.push_back(-position::L.TheX()/1000+giveRanD(position::L.TheX()/500));
                        varVec.push_back(-position::L.TheY()/1000+giveRanD(position::L.TheY()/500));
                        varVec.push_back(0.);

                        const position p=position(varVec);


                        start->Lbead->shiftParticle(0,p);

                        start->Lbead=nullptr;


                }
                break;
            }

             }


        }


    }

#ifdef WARMUP
if(!Warmup)
#endif
{
        SumofDisplacement=TSumOfdisplacement/measureCounter;
        SumOfPotential=TSumOfPotential/measureCounter;
        NumberOfParticles=TNumberOfParticles/measureCounter;
        NumberOfParticlesUp=TNumberOfParticlesUp/measureCounter;
        Wormlenght=1.*TWormlenght/measureCounter1;
        SumofWindingUp=TWindingUp/measureCounter;
        SumofWindingDown=TWindingDown/measureCounter;

}

}




block::~block()
{

}
