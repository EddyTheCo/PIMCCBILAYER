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




block::block(const size_t &NSweeps
             #ifdef USEROOT
             , TH2D * const Greens
             #endif
             )
{


double TSumOfdisplacement=0,TSumOfPotential=0,TNumberOfParticles=0,TNumberOfParticlesUp=0;
size_t TWormlenght=0,step=0,measureCounter=0,measureCounter1=0;
double TWindingUp=0,TWindingDown=0;



string svar;

while(step<NSweeps)
{



         if(Site::ThereIsAWorm)
        {

#ifdef WARMUP
             if(!Warmup)
#endif
             {
                measureCounter1++;
                TWormlenght+=Site::NInactiveLinks();

#ifdef USEROOT
                Greens->Fill(sqrt((Site::Rbead->pos-Site::Lbead->pos).norm()),abs(1.*Site::Rbead->TimeSliceOnBead-1.*Site::Lbead->TimeSliceOnBead));
#endif
             }
            switch ((!isGrandCanonical)?giveRanI(2):giveRanI(3)) {
            case 0:
            {
               //cout<<"closing worm"<<endl;
               //svar="closing worm";

                    Site::NCloseP++;


                    if(!(Site::cantClose(MBar)))
                    {

                            if(Site::Lbead->CloseWorm(0))
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
                    Site::NMoveP++;
                    Site::MoveWorm();
                     break;
            }
            case 2:
            {
              //cout<<"swap"<<endl;
              //svar="swap";

                Site::NSwapP++;
               if(Site::getNParti()>1)
               {
                    Site::PrepareSwap();
               }

               break;
            }
            case 3:
            {

                //cout<<"removeWorm"<<endl;
                //svar="removeWorm";
                Site::NRemoP++;

                if(Site::removeWorm())
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
                TSumOfdisplacement+=Site::TEnergy;
                TSumOfPotential+=Site::TPotential;
                TNumberOfParticles+=Site::getNParti();
                TNumberOfParticlesUp+=Site::getNParti_UP();
                TWindingUp+=Site::TWindingUp.normxy();
                TWindingDown+=Site::TWindingDown.normxy();
                measureCounter++;
            }

             switch ((isGrandCanonical)?giveRanI(2):giveRanI(1)) {
             case 0:
             {

               if(Site::getNParti())
               {
                 //cout<<"OpenWorm"<<endl;
                   // svar="OpenWorm";
                   const size_t posiTimes=giveRanI(NTimeSlices-1) ;
                   const size_t posiParti=giveRanI(Site::getNParti()-1);
                   const size_t var2=  giveRanI(MBar-2);

                   Site* const Ranbead=&(Site::theParticles->at(posiTimes).at(posiParti));
                   Site::NOpenP++;

                           Site::ThereIsAWorm= Ranbead->OpenWorm(var2,var2+1,0,Ranbead->pos);

               }

                break;
             }
             case 1:
            {
                 if(Site::getNParti())
                 {

                  //cout<<"wiggle"<<endl;
        //svar="wiggle";

                     const size_t posiTimes=giveRanI(NTimeSlices-1) ; //Choose a random time slice
                     const size_t posiParti=giveRanI(Site::getNParti()-1); //Choose the particle

                        Site::Lbead=&(Site::theParticles->at(posiTimes).at(posiParti)); //LBEAD is proposed (but dosent mean theres is a worm)
                        const size_t var2= giveRanI(MBar-3)+1;

                        Site::Rbead=Site::Lbead->searchBead(true,var2);
                        Site::Lbead->oldpos=Site::Lbead->pos;
                        Site::NWiggleP++;
                        Site::Lbead->Wiggle(0);

                        Site::Lbead=nullptr;
                        Site::Rbead=nullptr;

                }
                break;
            }

              case 2:
             {
          //       cout<<"insertworminclose "<<endl;
            //    svar="insertworminclose";
                 Site::NInsertP++;
                Site::insertWorm();
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
