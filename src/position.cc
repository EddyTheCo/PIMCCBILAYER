#include <cstdlib>
#include <fstream>
#include<iostream>
#include"position.hh"
#include"lattice.hh"

using namespace std;

const size_t d=ReadFromInput<size_t>(6);
const bool periodic=ReadFromInput<string>(14)=="periodic";

const position position::L=ReadFromInput<position>(7);
const double position::volumen=position::getVolumen();

  ifstream position::RestartConf(".restart.conf");
 ifstream position::PosiConf("Ini.conf");

ifstream* position::inFile=(restart)?&position::RestartConf:&position::PosiConf;

double position::foovar=position::foo();
position::position(const double &a)
{

    for(size_t i=0;i<d;i++)
    {
        x.push_back(a);
    }
}

position::position(const string str)
{
    if(str=="random")
        for(size_t i=0;i<d;i++)
        {
            x.push_back(Constants::giveRanD(L.x.at(i))-L.x.at(i)/2);
        }
    if(str=="ini")
        for(size_t i=0;i<d;i++)
        {
            if(restart)
            {
                double var;
                (*inFile)>>var;
                x.push_back(var);
            }
            else
            {
                //double var;

                //     (*inFile)>>var;
                 //     x.push_back(var);
                //x.push_back(0.);
               if(i==2)
	       {
		       x.push_back(0.0);
	       }
	       else
	       {
		       x.push_back(Constants::giveRanD(L.x.at(i))-L.x.at(i)/2);
	       }
            }

        }
}

position::position(const bool & isRight, const Site * const bead)
{
    const Site* edge;
    double ab;
    if(isRight)
    {
         edge=bead->Rbead;
         ab=edge->TimeSliceOnBead-static_cast<double>(bead->TimeSliceOnBead);
    }
    else
    {
         edge=bead->Lbead;
         ab=bead->TimeSliceOnBead-static_cast<double>(edge->TimeSliceOnBead);
    }

    if(ab<0)
        ab+=NTimeSlices;
    ab--;


     const auto aveVect=(bead->pos-edge->pos);
     const double mult=ab/(1.+ab);
     const auto var=variance*mult;

    for(size_t i=0;i<d;i++)
    {
        x.push_back(edge->pos.x.at(i)+mult*(aveVect.x.at(i))+Constants::giveRanDNormal(0,sqrt(var)));
    }



}
