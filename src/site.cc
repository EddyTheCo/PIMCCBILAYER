#include<iostream>
#include "site.hh"
#include"constants.hh"
#include"position.hh"
#include"potential.hh"
#include"lattice.hh"
#include<stack>
#include<map>
using namespace std;
using namespace Constants;


const bool Site::fourthOrder=ReadFromInput<string>(19)=="true";

int Site::Nparti_UpxNT=0;
int Site::NpartixNT=0;

const double landa=ReadFromInput<double>(4);
const double tao=ReadFromInput<double>(2);
const double beta=ReadFromInput<double>(3);
const double dplanes=ReadFromInput<double>(21);
double Site::mu=ReadFromInput<double>(9);
double Site::eta=ReadFromInput<double>(13);


double Site::TEnergyVar=0,Site::TPotentialVar=0;

const double variance=2*landa*tao;
ifstream Site::RestartPtrConf(".restartPtr.conf");
Site* Site::theZeta=nullptr;
const size_t NTimeSlices=beta/tao;
const size_t MBar=(ReadFromInput<size_t>(15)<1||ReadFromInput<size_t>(15)>NTimeSlices)?NTimeSlices:ReadFromInput<size_t>(15);

double Site::TEnergy=0,Site::TPotential=0;

const potential Site::ThePotential=potential();




bool Site::ThereIsAWorm=false;
array<vector<Site>,100000>* Site::theParticles=nullptr;

Site* Site::Lbead=nullptr;
Site* Site::Rbead=nullptr;
position Site::TWindingUp=position(0.);
position Site::TWindingDown=position(0.);
position Site::TWindingVarUp=position(0.);
position Site::TWindingVarDown=position(0.);
size_t Site::NClose=1,Site::NWiggle=1,Site::NWiggleP=1,Site::NShift=1,Site::NShiftP=1,Site::NOpen=1,Site::NMove=1,Site::NSwap=1,Site::NCloseP=1,Site::NOpenP=1,Site::NMoveP=1,Site::NSwapP=1,Site::NInsertP=1,Site::NInsert=1,Site::NRemoP=1,Site::NRemo=1;

Site::Site(const size_t i, const size_t j, const bool UP, const position var, const bool act):active(act),UPplane(UP),pos(var),oldpos(pos),TimeSliceOnBead(j),ParticleOnBead(i),left(nullptr),right(nullptr),up(nullptr),down(nullptr)
{

if(UP)Nparti_UpxNT++;
NpartixNT++;
}



Site::~Site(){
    if(UPplane)Nparti_UpxNT--;
    NpartixNT--;
}

bool Site::OpenWorm(const size_t step, const size_t ab, double dU, const position& start)
{

    dU+=-mu*tao;

    if(step>0)
    {
        double U=0;
        right->ChangeInU(true,dU,U);

        if(right->OpenWorm(step-1,ab,dU,start))
        {            
            Lbead=this;
            right->active=false;
            const auto Dist=pos-right->pos;
            TEnergy-=Dist.normxy();

            (UPplane)?(TWindingUp=TWindingUp+Dist):(TWindingDown=TWindingDown+Dist);
            TPotential-=U;
            return true;
        }
        return false;

    }
    else
    {


        if(eta*getNParti()/(propagator(start,right->pos,ab,-dU)*position::volumen)>giveRanD(1.))
        {

            const auto Dist=pos-right->pos;
            TEnergy-=Dist.normxy();
            (UPplane)?(TWindingUp=TWindingUp+Dist):(TWindingDown=TWindingDown+Dist);
            Rbead=this->right;
            ThereIsAWorm=true;
            Lbead=this;
            NOpen++;
            return true;
        }
        return false;
    }

}
bool Site::CloseWorm(double dU)
{

    dU+=mu*tao;

    if(this->TimeSliceOnBead!=Rbead->left->TimeSliceOnBead)
    {

        right->pos=position(true,this);
        double U=0;
        right->ChangeInU(false,dU,U);


        if(right->CloseWorm(dU))
        {
            right->active=true;
            const auto Dist=right->pos-pos;
            TEnergy+=Dist.normxy();

            (UPplane)?(TWindingUp=TWindingUp+Dist):(TWindingDown=TWindingDown+Dist);

            TPotential+=U;

            return true;

        }
        return false;

    }
    else
    {

        if(propagator(Lbead->pos,Rbead->pos,NInactiveLinks(),dU)*position::getVolumen()/(eta*getNParti())>giveRanD(1.))
        {

            const auto Dist=right->pos-pos;
            TEnergy+=Dist.normxy();

           (UPplane)?(TWindingUp=TWindingUp+Dist):(TWindingDown=TWindingDown+Dist);

            Rbead=nullptr;
            Lbead=nullptr;
            ThereIsAWorm=false;
            NClose++;
            return true;
        }
        return false;
    }


}
bool Site::Wiggle(double dU)
{

right->oldpos=right->pos;

    if(this->TimeSliceOnBead!=Rbead->left->TimeSliceOnBead)
    {
        const auto Dist=oldpos-right->pos;
        TEnergyVar-=Dist.normxy();

        (UPplane)?TWindingVarUp=TWindingVarUp+Dist:TWindingVarDown=TWindingVarDown+Dist;

        double U=0;
        right->ChangeInU(true,dU,U);

        TPotentialVar+=U;

        right->pos=position(true,this);


        right->ChangeInU(false,dU,U);

        if(right->Wiggle(dU))
        {
            const auto Dist=right->pos-pos;
            TEnergy+=Dist.normxy();

            (UPplane)?(TWindingUp=TWindingUp+Dist):(TWindingDown=TWindingDown+Dist);

            TPotential+=U;

            return true;
        }
        right->pos=right->oldpos;
        return false;
    }
    else
    {
        const auto Dist=oldpos-right->pos;
        TEnergyVar-=Dist.normxy();
        (UPplane)?TWindingVarUp=TWindingVarUp+Dist:TWindingVarDown=TWindingVarDown+Dist;

        if(exp(dU)>giveRanD(1.))
        {

            TPotential-=TPotentialVar;
            const auto Dist=right->pos-pos;
            TEnergy+=TEnergyVar+Dist.normxy();

            (UPplane)?(TWindingUp=TWindingUp+Dist+TWindingVarUp):(TWindingDown=TWindingDown+Dist+TWindingVarDown);

            Lbead=nullptr;
            Rbead=nullptr;
            TPotentialVar=0;
            TEnergyVar=0;
            NWiggle++;
            (UPplane)?TWindingVarUp=position(0.):TWindingVarDown=position(0.);
            return true;
        }
        TPotentialVar=0;
        TEnergyVar=0;
        (UPplane)?TWindingVarUp=position(0.):TWindingVarDown=position(0.);

        return false;
    }
}



bool Site::deleteToRight(const size_t step, double dU)
{

dU+=-mu*tao;
double U=0;
ChangeInU(true,dU,U);


        if(step>0)
        {
            if(right->deleteToRight(step-1,dU))
            {
                active=false;
                const auto Dist=pos-right->pos;
                TEnergy-=Dist.normxy();

                (UPplane)?(TWindingUp=TWindingUp+Dist):(TWindingDown=TWindingDown+Dist);

                TPotential-=U;
                return true;
            }
            return false;

        }
        else
        {

            if(this==Lbead)dU+=mu*tao;



            if(exp(dU)>giveRanD(1.))
            {

                if(this!=Lbead)
                {
                    const auto Dist=pos-right->pos;
                    TEnergy-=Dist.normxy();

                    (UPplane)?(TWindingUp=TWindingUp+Dist):(TWindingDown=TWindingDown+Dist);

                }
                active=false;
                TPotential-=U;
                Rbead=this->right;
                NMove++;
                return true;
            }
            return false;
        }


}
bool Site::deleteToLeft(const size_t step,double dU)
{

dU+=-mu*tao;
double U=0;
ChangeInU(true,dU,U);
            if(step>0)
            {
                if(left->deleteToLeft(step-1,dU))
                {
                    active=false;
                    const auto Dist=left->pos-pos;
                    TEnergy-=Dist.normxy();

                    (UPplane)?(TWindingUp=TWindingUp+Dist):(TWindingDown=TWindingDown+Dist);

                    TPotential-=U;
                    return true;
                }
                return false;

            }
            else
            {
                NMoveP++;

                if(exp(dU)>giveRanD(1.))
                {

                    active=false;
                    TPotential-=U;
                    const auto Dist=left->pos-pos;
                    TEnergy-=Dist.normxy();

                    (UPplane)?(TWindingUp=TWindingUp+Dist):(TWindingDown=TWindingDown+Dist);

                    Lbead=this->left;
                    NMove++;
                    return true;
                }
                return false;
            }


}
bool Site::insertToRight(const size_t step,double dU)
{
    static bool aParticleisInserted=false;
    if (right->active)
    {
        if(!isGrandCanonical)
            return false;
        else {
            insertParticle(UPplane);
            aParticleisInserted=true;
            right=theParticles->at(TimeSliceOnBead).back().right;
            auto var=right->left;
            right->left=this;
            Rbead->left=var;
            var->right=Rbead;
        }
    }

dU+=mu*tao;
right->pos=position(pos,variance);
double U=0;
right->ChangeInU(false,dU,U);


            if(step>0)
            {
                if(right->insertToRight(step-1,dU))
                {
                    right->active=true;
                    const auto Dist=right->pos-pos;
                    TEnergy+=Dist.normxy();

                    (UPplane)?(TWindingUp=TWindingUp+Dist):(TWindingDown=TWindingDown+Dist);

                    TPotential+=U;
                    return true;
                }
                return false;

            }
            else
            {

                NMoveP++;

                if(exp(dU)>giveRanD(1.))
                {

                    right->active=true;
                    const auto Dist=right->pos-pos;
                    TEnergy+=Dist.normxy();

                    (UPplane)?(TWindingUp=TWindingUp+Dist):(TWindingDown=TWindingDown+Dist);

                    TPotential+=U;
                    Lbead=this->right;
                    if(aParticleisInserted)
                    {
                        aParticleisInserted=false;
                    }
                    NMove++;
                    return true;
                }

                if(aParticleisInserted)//but the move is not accepted
                {
                    auto var=&(theParticles->at(Rbead->TimeSliceOnBead).back());
                    auto varLeft=var->left;
                    var->left=Rbead->left;
                    var->left->right=var;
                    Rbead->left=varLeft;
                    varLeft->right=Rbead;
                    removeLastParticle();
                    aParticleisInserted=false;
                }
                return false;

            }

}
bool Site::insertToLeft(const size_t step,double dU)
{


    static bool aParticleisInserted=false;

    if (left->active)
    {
        if(!isGrandCanonical)
            return false;
        else {

            insertParticle(UPplane);

            aParticleisInserted=true;
            left=theParticles->at(TimeSliceOnBead).back().left;
            auto var=left->right;
            left->right=this;
            Lbead->right=var;
            var->left=Lbead;
        }
    }

dU+=mu*tao;

    left->pos=position(pos,variance);


    double U=0;
    left->ChangeInU(false,dU,U);
            if(step>0)
            {

                if(left->insertToLeft(step-1,dU))
                {
                    left->active=true;
                    const auto Dist=pos-left->pos;
                    TEnergy+=Dist.normxy();


                    (UPplane)?(TWindingUp=TWindingUp+Dist):(TWindingDown=TWindingDown+Dist);


                    TPotential+=U;
                    return true;
                }
                return false;

            }
            else
            {

                NMoveP++;

                if(exp(dU)>giveRanD(1.))
                {

                    left->active=true;
                    const auto Dist=pos-left->pos;
                    TEnergy+=Dist.normxy();


                    (UPplane)?(TWindingUp=TWindingUp+Dist):(TWindingDown=TWindingDown+Dist);


                    TPotential+=U;
                    Rbead=this->left;
                    if(aParticleisInserted)
                    {

                        aParticleisInserted=false;
                    }

                    NMove++;
                    return true;
                }
                if(aParticleisInserted)//but the move is not accepted
                {
                    auto var=&(theParticles->at(Lbead->TimeSliceOnBead).back());
                    auto varRight=var->right;
                    var->right=Lbead->right;
                    var->right->left=var;
                    Lbead->right=varRight;
                    varRight->left=Lbead;
                    removeLastParticle();
                    aParticleisInserted=false;
                }
                return false;
            }
}






Site* Site::chooseTheBead(double &SumI, const size_t& vae, const Site * const startBead)const
{
     multimap <double, Site*> prop;

    auto ptr=this->up;
    double var;
    while(ptr!=this)
    {


        if(ptr->active&&ptr!=Rbead&&ptr!=Lbead&&ptr->UPplane==startBead->UPplane)
        {

            var=propagator(startBead->pos,ptr->pos,vae);
            if(var>0){
                SumI+=var;
                prop.insert(pair <double, Site*> (var, ptr));
            }

        }
        ptr=ptr->up;
    }
    if(prop.size()==0)
        return nullptr;
    auto ran=giveRanD(1);

    double pr=0;
    map<double, Site*>::reverse_iterator rit;
      for (rit=prop.rbegin(); rit!=prop.rend(); ++rit)
        {
            pr+=rit->first/SumI;
            if(pr>ran)
            {
                return rit->second;
            }

        }
}

bool Site::swap(Site* const zeta, const double& SumI, const double& SumZ, double dU, const bool &isRight)
{



if(isRight)
{
    static bool aParticleisInserted=false;
    static Site* ri;
    if(this->TimeSliceOnBead!=Rbead->left->TimeSliceOnBead)
    {
        //cout<<"pos="<<pos<<endl;
        if (right->active)
        {
           if(!isGrandCanonical)
           {
               theZeta=nullptr;
                return false;
           }
            else {
                insertParticle(UPplane);
                aParticleisInserted=true;
                //cout<<"seinserto aparticle"<<endl;
                ri=right;
                right=theParticles->at(TimeSliceOnBead).back().right;
                auto var=right->left;
                right->left=this;
                ri->left=var;
                var->right=ri;
            }
        }

        right->pos=position(true,this);
        double Ualpha=0,Uzeta=0;
        theZeta=zeta->right;
        right->ChangeInU(false,dU,Ualpha);
        zeta->right->ChangeInU(true,dU,Uzeta);

        if(right->swap(zeta->right,SumI,SumZ,dU,true))
        {
            this->right->active=true;
            zeta->right->active=false;
            const auto Dist1=right->pos-pos;
            const auto Dist2=zeta->pos-zeta->right->pos;
            TEnergy+=Dist1.normxy()-Dist2.normxy();

            (UPplane)?(TWindingUp=TWindingUp+Dist1+Dist2):(TWindingDown=TWindingDown+Dist1+Dist2);

            TPotential+=Ualpha;
            TPotential-=Uzeta;

            return true;
        }

        return false;

    }
    else
    {

         if(exp(dU)*SumI/SumZ>giveRanD(1.))
         {
             Site* const prev=this->right;
             this->right=zeta->right;
             const auto Dist1=right->pos-pos;
             const auto Dist2=zeta->pos-zeta->right->pos;
             TEnergy+=Dist1.normxy()-Dist2.normxy();

             (UPplane)?(TWindingUp=TWindingUp+Dist1+Dist2):(TWindingDown=TWindingDown+Dist1+Dist2);

            zeta->right=prev;
            zeta->right->left=zeta;
            this->right->left=this;
            NSwap++;
            if(aParticleisInserted)
            {
                aParticleisInserted=false;
            }
            theZeta=nullptr;
            return true;
         }
         if(aParticleisInserted)//but the move is not accepted
         {

             const auto var=&(theParticles->at(ri->TimeSliceOnBead).back());
             auto varLeft=var->left;
             var->left=ri->left;
             var->left->right=var;
             ri->left=varLeft;
             varLeft->right=ri;
             removeLastParticle();
             aParticleisInserted=false;

         }
         theZeta=nullptr;
         return false;
    }

}
else {
    static bool aParticleisInserted=false;
    static Site* le;
    if(this->TimeSliceOnBead!=Lbead->right->TimeSliceOnBead)
    {

        if (left->active)
        {
            if(!isGrandCanonical)
            {
                theZeta=nullptr;
                return false;

            }
            else {
                insertParticle(UPplane);
                aParticleisInserted=true;

                le=left;
                left=theParticles->at(TimeSliceOnBead).back().left;
                auto var=left->right;
               left->right=this;
                le->right=var;
                var->left=le;

            }
        }

        left->pos=position(false,this);
        double Ualpha=0,Uzeta=0;
        theZeta=zeta->left;
        left->ChangeInU(false,dU,Ualpha);
        zeta->left->ChangeInU(true,dU,Uzeta);

        if(left->swap(zeta->left,SumI,SumZ,dU,false))
        {

            this->left->active=true;
            zeta->left->active=false;
            const auto Dist1=pos-left->pos;
            const auto Dist2=zeta->left->pos-zeta->pos;
            TEnergy+=Dist1.normxy()-Dist2.normxy();

            (UPplane)?(TWindingUp=TWindingUp+Dist1+Dist2):(TWindingDown=TWindingDown+Dist1+Dist2);

            TPotential+=Ualpha;
            TPotential-=Uzeta;

            return true;
        }

        return false;
    }
    else
    {


        if(exp(dU)*SumI/SumZ>giveRanD(1.))
        {

            Site * const prev=this->left;
            this->left=zeta->left;
            const auto Dist1=pos-left->pos;
            const auto Dist2=zeta->left->pos-zeta->pos;
            TEnergy+=Dist1.normxy()-Dist2.normxy();

            (UPplane)?(TWindingUp=TWindingUp+Dist1+Dist2):(TWindingDown=TWindingDown+Dist1+Dist2);

            zeta->left=prev;
            zeta->left->right=zeta;
            this->left->right=this;
            NSwap++;
            if(aParticleisInserted)
            {


                aParticleisInserted=false;
            }
            theZeta=nullptr;
            return true;
        }
        if(aParticleisInserted)//but the move is not accepted
        {

            const auto var=&(theParticles->at(le->TimeSliceOnBead).back());
            auto varRight=var->right;
            var->right=le->right;
            var->right->left=var;
            le->right=varRight;
            varRight->left=le;
            removeLastParticle();
            aParticleisInserted=false;

        }
        theZeta=nullptr;
        return false;
    }

}





}


bool Site::shiftParticle(double dU, const position& shift)const
{
    right->oldpos=right->pos;

    double U=0;
    right->ChangeInU(true,dU,U);
    TPotentialVar+=U;

    //cout<<"right="<<right->ParticleOnBead<<" "<<right->TimeSliceOnBead<<" "<<right->pos<<endl;
    right->pos=right->pos+shift;
    //cout<<"right="<<right->ParticleOnBead<<" "<<right->TimeSliceOnBead<<" "<<right->pos<<endl;
    right->ChangeInU(false,dU,U);

        if(this!=Lbead->left)
        {

            if(right->shiftParticle(dU,shift))
            {
                TPotential+=U;
                //cout<<oldpos-right->oldpos<<" "<<pos-right->pos<<endl;
                return true;
            }
            right->pos=right->oldpos;
            return false;
        }
        else
        {


            if(exp(dU)>giveRanD(1.))
            {
                NShift++;
                TPotential+=U;
                TPotential-=TPotentialVar;
                TPotentialVar=0.;
                //cout<<oldpos-right->oldpos<<" "<<pos-right->pos<<endl;
               return true;
            }
            right->pos=right->oldpos;
            TPotentialVar=0.;
            return false;
        }


}
