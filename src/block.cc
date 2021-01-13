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


size_t h=0;
size_t War=Warmup;
static size_t CorrectNpart=0;


while(step<NSweeps)
{


    if(!(h%1000)&&War&&isGrandCanonical)
    {
        if(start->NParti_<War)
        {
            Site::mu+=1;
        }
        if(start->NParti_>War)
        {
            Site::mu-=1;
        }
        if(start->NParti_==War)
        {
            CorrectNpart++;
            if(CorrectNpart>=100&&start->Nparti_UpxNT/NTimeSlices==Warmup/2&&start->NParti_==Warmup)
            {
                War=0;
                isGrandCanonical=0;

            }
        }

        cout<<"mu="<<Site::mu<<" eta="<<Site::eta<<" Npar="<<start->NParti_<<" NPartUp="<<start->Nparti_UpxNT/NTimeSlices<<" C="<<CorrectNpart<<endl;


    }
    h++;


         if(start->ThereIsAWorm)
        {

             if(!Warmup)
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
//               cout<<"closing worm"<<endl;
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
//               cout<<"MoveWorm"<<endl;
                    start->MoveWorm();
                     break;
            }
            case 2:
            {
//              cout<<"swap"<<endl;
                start->NSwapP++;
               if(start->NParti_>1)
               start->PrepareSwap();

               break;
            }
            case 3:
            {

        //        cout<<"removeWorm"<<endl;
                start->removeWorm();
                break;
            }



            }

        }
        else
        {
            if(!Warmup)
            {
                TSumOfdisplacement+=start->TEnergy;
                TSumOfPotential+=start->TPotential;
                TNumberOfParticles+=start->NParti_;
                TNumberOfParticlesUp+=start->Nparti_UpxNT/NTimeSlices;
                (d>1)?TWindingUp+=start->TWindingUp.normxy():TWindingUp+=start->TWindingUp.norm();
                TWindingDown+=start->TWindingDown.normxy();
                measureCounter++;
            }

             switch ((isGrandCanonical)?giveRanI(3):giveRanI(2)) {
             case 0:
             {

               if(start->NParti_)
               {
//                 cout<<"OpenWorm"<<endl;
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

//                  cout<<"wiggle"<<endl;

                     const size_t posiTimes=giveRanI(NTimeSlices-1) ; //Choose a random time slice
                     const size_t posiParti=giveRanI(start->NParti_-1); //Choose the particle

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

              case 3:
             {
              //   cout<<"insertworminclose "<<endl;
                 start->insertWorm();
                 break;
             }
             case 2:
            {
                 if(start->NParti_)
                 {

//                 cout<<"ShiftParticle"<<endl;


                     const size_t posiParti=giveRanI(start->NParti_-1); //Choose the particle

                        start->Lbead=&(particles->at(0).at(posiParti)); //LBEAD is proposed (but dosent mean theres is a worm)

//cout<<"Lbead="<<start->Lbead->ParticleOnBead<<" "<<start->Lbead->TimeSliceOnBead<<endl;
                        vector<double> varVec;

                        varVec.push_back(giveRanDNormal(0,position::L.TheX()));
                        varVec.push_back(giveRanDNormal(0,position::L.TheY()));
                        varVec.push_back(0.);

                        const position p=position(varVec);
//cout<<"vectShift="<<p<<endl;

                        start->Lbead->shiftParticle(0,p);

                        start->Lbead=nullptr;


                }
                break;
            }

             }


        }

    }


if(!Warmup)
{
        SumofDisplacement=TSumOfdisplacement/measureCounter;
        SumOfPotential=TSumOfPotential/measureCounter;
        NumberOfParticles=TNumberOfParticles/measureCounter;
        NumberOfParticlesUp=TNumberOfParticlesUp/measureCounter;
        Wormlenght=1.*TWormlenght/measureCounter1;
        SumofWindingUp=TWindingUp/measureCounter;
        SumofWindingDown=TWindingDown/measureCounter;
}
if(!War&&Warmup)
{
    Warmup=0;
}


}




block::~block()
{

}
