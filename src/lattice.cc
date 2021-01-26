#include<iostream>
#include <iomanip>
#include "site.hh"
#include"constants.hh"
#include"lattice.hh"
#include"block.hh"
#include<stack>
#include <cstring>
#include <cstdio>


const bool restart=ReadFromInput<string>(10)=="restart";
#ifdef USEROOT
#include<TFile.h>
#include<TH1.h>
#include<TH2D.h>
#include<TH3D.h>


TFile *lattice::RootFile = (restart)?new TFile("RootFile.root","UPDATE"):new TFile("RootFile.root","RECREATE");
TH1* lattice::hpos = nullptr;
TH1D* lattice::PCFUp = nullptr;
TH1D* lattice::PCFDown = nullptr;
TH1D* lattice::PCFMix = nullptr;
TH2D* lattice::Greens = nullptr;
#endif

using namespace std;
using namespace Constants;


const size_t SAMPLING=ReadFromInput<size_t>(20);
const size_t NPartiIni=(restart)?ReadFromInput<size_t>(1,".restart.conf"):ReadFromInput<size_t>(1);
ofstream  lattice::thesweep("sweep",(restart)?std::ofstream::out | std::ofstream::app:std::ofstream::out),lattice::theratios("ratios",(restart)?std::ofstream::out | std::ofstream::app:std::ofstream::out);

bool isGrandCanonical=ReadFromInput<string>(11)=="GrandCanonical";
size_t Warmup=ReadFromInput<int>(22);
const size_t lattice::NRep=ReadFromInput<size_t>(5);

const size_t lattice::NSweeps=ReadFromInput <size_t> (8);









array<vector<Site>,100000>*  const lattice::grid=new array<vector<Site>,100000>;
lattice::lattice()
{
cout<<"init Lattice"<<endl;

    for(size_t i=0;i<NTimeSlices;i++)
    {
            grid->at(i).reserve(10000);

            for(size_t j = 0; j<NPartiIni; j++)
            {
                 grid->at(i).push_back(Site(j,i));

            }


    }

#ifdef USEROOT
 if (RootFile->IsZombie()) {
      cout << "Error opening file" << endl;
  }
    if(!RootFile->IsOpen())
    {
        RootFile=new TFile("RootFile.root","RECREATE");
    }
    initialize_histos();
#endif



cout<<"end Lattice"<<endl;
}
#ifdef USEROOT
TVectorD *lattice::v=nullptr;
 void lattice::initialize_histos(void)
{

    if(restart)
    {
        Greens = (TH2D*)gDirectory->Get("Greens");
        PCFUp = (TH1D*)gDirectory->Get("PCFUp");
        PCFDown = (TH1D*)gDirectory->Get("PCFDown");
        PCFMix = (TH1D*)gDirectory->Get("PCFMix");


        if(!Greens)Greens=new TH2D("Greens","",5000,0,sqrt((position::L).norm()),NTimeSlices,-0.5,NTimeSlices-0.5);
        if(!PCFUp)PCFUp=new TH1D("PCFUp","",5000,0,sqrt((position::L).norm()-(position::L).TheZ()*(position::L).TheZ()));
        if(!PCFDown)PCFDown=new TH1D("PCFDown","",5000,0,sqrt((position::L).norm()-(position::L).TheZ()*(position::L).TheZ()));
        if(!PCFMix)PCFMix=new TH1D("PCFMix","",5000,0,sqrt((position::L).norm()-(position::L).TheZ()*(position::L).TheZ()));
        v = (TVectorD*)gDirectory->Get("v");
        if(!v)
        {
            TVectorD v(1);
            v[0]=0.;
            v.Write("v");
        }

    }
    else
    {
        Greens=new TH2D("Greens","",5000,0,sqrt((position::L).norm()),NTimeSlices,-0.5,NTimeSlices-0.5);
        PCFUp=new TH1D("PCFUp","",5000,0,sqrt((position::L).norm()-(position::L).TheZ()*(position::L).TheZ()));
        PCFDown=new TH1D("PCFDown","",5000,0,sqrt((position::L).norm()-(position::L).TheZ()*(position::L).TheZ()));
        PCFMix=new TH1D("PCFMix","",5000,0,sqrt((position::L).norm()-(position::L).TheZ()*(position::L).TheZ()));

        TVectorD v(1);
        v[0]=0.;
        v.Write("v");
    }


    v = (TVectorD*)gDirectory->Get("v");

    if(d==3)
    hpos=new TH3D("pos","",100,-position::L.x.at(0)/2,position::L.x.at(0)/2,100,-position::L.x.at(1)/2,position::L.x.at(1)/2,100,-position::L.x.at(2)/2,position::L.x.at(2)/2);
    if(d==2)
    hpos=new TH2D("pos","",1000,-position::L.x.at(0)/2,position::L.x.at(0)/2,1000,-position::L.x.at(1)/2,position::L.x.at(1)/2);
    if(d==1)
    hpos=new TH1D("pos","",10000,-position::L.x.at(0)/2,position::L.x.at(0)/2);

}
#endif




void lattice::setup()const
{
cout<<"init Lattice Setup"<<endl;
    grid->at(0).at(0).setLattice(grid);


    for(size_t i=0;i<NTimeSlices;i++)
    {
            for(size_t j = 0; j<NPartiIni; j++)
            {
                 grid->at(i).at(j).setPtr();
            }
    }

    for(size_t i=0;i<NTimeSlices;i++)
    {

            for(size_t j = 0; j<NPartiIni; j++)
            {

                 grid->at(i).at(j).totalEnergy();
            }


    }
cout<<"end Lattice Setup"<<endl;
}


#ifdef WARMUP
void lattice::Warm() const
{

    double Np=0;
    while(Np<Warmup&&isGrandCanonical)
    {
        Site::mu++;

        const auto myBlock=block(grid,NTimeSlices,100
#ifdef USEROOT
                                 ,nullptr
#endif
                                     );



        Np=grid->at(0).at(0).NParti_;
        cout<<"Nparticles="<<Np<<" mu="<<Site::mu<<" eta="<<Site::eta<<endl;
            system(("sed -i 's/^.*\\#mu\\b.*$/" + to_string(Site::mu)   + "              \\#mu/' input").c_str());
    }

    while((grid->at(0).at(0).NParti_!=Warmup||grid->at(0).at(0).Nparti_UpxNT/NTimeSlices!=Warmup/2)&&isGrandCanonical)
    {
        const auto myBlock=block(grid,NTimeSlices,1
#ifdef USEROOT
                                 ,nullptr
#endif
                                     );

        cout<<"Nparticles="<<grid->at(0).at(0).NParti_<<" NPUP="<<grid->at(0).at(0).Nparti_UpxNT/NTimeSlices<<endl;

    }



system(("sed -i 's/^.*\\#eta\\b.*$/" + to_string(Site::eta)   + "              \\#eta factor increase for more oppening/' input").c_str());

    system(("sed -i 's/^.*\\#Warmup\\b.*$/" +to_string(0) +     "              \\#Warmup/' input").c_str());
    system("sed -i '/Canonical/s/.*/Canonical       #GrandCanonical or Canonical/' input");
    system("sed -i '/start or restart/s/.*/restart           #start or restart/' input");


    Warmup=0;
    isGrandCanonical=false;
    ofstream RestartConf(".restartVAR.conf");
    ofstream RestartPtrConf(".restartPtrVAR.conf");

    RestartConf<<grid->at(0).at(0).NParti_<<" #Particles"<<endl;
    RestartConf.precision(12);


    for(size_t i=0;i<NTimeSlices;i++)
    {

            for(size_t j = 0; j<grid->at(0).at(0).NParti_; j++)
            {

                const auto var=grid->at(i).at(j);
                RestartConf<<var.pos;

                RestartPtrConf<<var.left->ParticleOnBead<<" "<<var.left->TimeSliceOnBead<<" "<<var.right->ParticleOnBead<<" "<<var.right->TimeSliceOnBead<<" ";
            }
    }
    RestartConf.close();
    RestartPtrConf.close();
    rename(".restartVAR.conf", ".restart.conf");
    rename(".restartPtrVAR.conf", ".restartPtr.conf");





}
#endif
void lattice::move()const
{

cout<<"init Lattice Move"<<endl;
if(!restart)thesweep<< left << setw(16) <<"KEnergy"<< left << setw(16) <<"PEnergy"<< left << setw(16) <<"TEnergy"<< left << setw(16) <<"WoLen"<< left << setw(16) <<"SFDUp"<< left << setw(16) <<"SFDDown"<< left << setw(16) <<"NpartiUp"<< right << setw(16) <<"NParti"<</* right << setw(12) <<"S(k,w=0)"<<*/endl;

theratios<< left << setw(12)<<"Ropen"<<left << setw(12)<<"RClose"<<left << setw(12)<<"Rmove"<<left << setw(12)<<"Rswap"<<left << setw(12)<<"RWiggle"<<left << setw(12)<<"Rshift"<<left << setw(12)<<"RInsert"<<left << setw(12)<<"RRmove"<<endl;
    for(size_t step=0;step<NRep;step++)
    {
        //grid->at(0).at(0).printLattice("");
        const auto myBlock=block(grid,NTimeSlices,NSweeps
#ifdef USEROOT
                                 ,Greens
#endif
                                 );


cout<<"finished block "<<step<<endl;
            PrintConfiguration(
#ifdef USEROOT
                        (*v)[0]+1
#endif
  );


thesweep << left << setw(16) << myBlock.getKineticEnergy()<<" "<< left << setw(16) <<myBlock.SumOfPotential/(NTimeSlices)
         <<" "<< left << setw(16) << myBlock.getKineticEnergy()+myBlock.SumOfPotential/(NTimeSlices)
          <<" "<<left << setw(16)    << myBlock.Wormlenght/(NTimeSlices-1)<<" "<<left << setw(16)    << myBlock.getSuperfluidDensityUp()
         <<" "<<left << setw(16)    << myBlock.getSuperfluidDensityDown()<<" "<<left << setw(16)    << myBlock.NumberOfParticlesUp<<" "<<right << setw(16)    <<myBlock.NumberOfParticles<</*" "<<right << setw(12)    <<StructureFactor<<*/endl;

if(grid->at(0).size())
{
    theratios<< left << setw(12) <<getOpenRatio()<< left << setw(12) <<getCloseRatio()<< left << setw(12) <<getMoveRatio()<< left << setw(12)
                <<getSwapRatio()<< left << setw(12) <<getWiggleRatio()<< left << setw(12) <<getShiftRatio()<< left << setw(12) <<getInsertRatio()<< left << setw(12) <<getRemoveRatio()<<endl;
}



    }

#ifdef USEROOT
    RootFile->Close();
#endif
    thesweep.close();
    theratios.close();
    cout<<"end Lattice Move"<<endl;
}

void lattice::PrintConfiguration (
        #ifdef USEROOT
        const size_t step=0
        #else
        void
        #endif
        ) const
{

   ofstream RestartConf(".restartVAR.conf");
   ofstream RestartPtrConf(".restartPtrVAR.conf");
#ifdef SAVECONF
   static size_t step=0;
   std::ofstream Data;
     Data.open ("conf.dat", std::ofstream::out | std::ofstream::app);
   static int var=1;
   if(var)
   {
       Data<<"#Beta="<<beta<<" Nparticles="<<grid->at(0).at(0).NParti_<<" NtimesLices="<<NTimeSlices<<endl;
       var=0;
   }
#endif
   size_t NP=0;
   if(grid->at(0).size())
   {
       NP=grid->at(0).at(0).NParti_;

   }
   RestartConf<<NP<<" #Particles"<<endl;
   RestartConf.precision(12);

    if(NP)for(size_t i=0;i<NTimeSlices;i++)
   {

           for(size_t j = 0; j<NP; j++)
           {

               const auto var=grid->at(i).at(j);
               RestartConf<<var.pos;
#ifdef SAVECONF
               if(!(step%SAMPLING)&&step!=0)Data<<var.pos;
#endif

               RestartPtrConf<<var.left->ParticleOnBead<<" "<<var.left->TimeSliceOnBead<<" "<<var.right->ParticleOnBead<<" "<<var.right->TimeSliceOnBead<<" ";
#ifdef USEROOT

                if(d==3)
                ((TH3D*)hpos)->Fill(grid->at(i).at(j).pos.perio(0),grid->at(i).at(j).pos.perio(1),grid->at(i).at(j).pos.perio(2));
                if(d==2)
                ((TH2D*)hpos)->Fill(grid->at(i).at(j).pos.perio(0),grid->at(i).at(j).pos.perio(1));
                if(d==1)
                ((TH1D*)hpos)->Fill(grid->at(i).at(j).pos.perio(0));
                if(Site::NParti_>4)
                for(size_t k = 0; k<grid->at(0).at(0).NParti_; k++)
                {
                    if(j!=k)
                    {


                    if(((grid->at(i).at(j).pos.TheZ()>0.)&&(grid->at(i).at(k).pos.TheZ()<0.))||((grid->at(i).at(j).pos.TheZ()<0.)&&(grid->at(i).at(k).pos.TheZ()>0.)))
                    {
                        const double dis=sqrt((grid->at(i).at(j).pos-grid->at(i).at(k).pos).norm()-dplanes*dplanes);
                        PCFMix->Fill(dis);
                    }
                    else
                    {
                        const double dis=sqrt((grid->at(i).at(j).pos-grid->at(i).at(k).pos).norm());
                        if(grid->at(i).at(j).pos.TheZ()>0)
                        {
                            PCFUp->Fill(dis);
                        }
                        else
                        {
                            PCFDown->Fill(dis);
                        }
                    }
                    }
                }

#endif
            }
#ifdef SAVECONF
           if(!(step%SAMPLING)&&step!=0)Data<<endl;
#endif
   }
#ifdef SAVECONF
   if(!(step%SAMPLING)&&step!=0)Data<<endl;
   Data.close();
   step++;
#endif
#ifdef USEROOT
    (*v)[0] = step;



       (*v).Write("v",TObject::kOverwrite);

       hpos->Write(("pos" + to_string(step%1000)).c_str(),TObject::kOverwrite);
        hpos->Reset("ICESM");

   gDirectory->Write("", TObject::kOverwrite);
#endif


Constants::saveRandom();
   RestartConf.close();
   RestartPtrConf.close();
   rename(".restartVAR.conf", ".restart.conf");
   rename(".restartPtrVAR.conf", ".restartPtr.conf");




}

