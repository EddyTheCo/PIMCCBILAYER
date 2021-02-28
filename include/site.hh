#include <cstdlib>
#include<stack>
#include <fstream>
#include<iostream>
#include<array>

using namespace std;


#ifndef SITE_H
#define SITE_H 1

#include"potential.hh"

const extern double beta,tao,landa,variance,dplanes;
const extern size_t NTimeSlices,d, MBar;
extern const bool restart;


class Site
{
	public:

    ~Site();
    Site(const size_t i,const size_t j,const bool UP,const position var,const bool act );

    static Site* Lbead;
    static Site* Rbead;
    friend ostream & operator << (ostream &out, const Site & obj)
    {
                if(obj.active)
                out<<obj.pos;
        return out;
    }
    void setLattice(array<vector<Site>,100000>* const particles)const
    {
        theParticles=particles;
    }
    inline static size_t  getNParti_UP(void)
    {
        return Nparti_UpxNT/NTimeSlices;
    }
    inline  static size_t getNParti(void)
    {
        return NpartixNT/NTimeSlices;
    }
    void static PrepareSwap(void)
    {

        const auto oldRbead=Rbead;
        double SumI=0,SumZ=0;
        const auto oldLbead=Lbead;


            const auto vae=giveRanI(MBar-2);
            Site* alpha, *zeta;

            if ( giveRanI(1) )
            {


                const auto NewRbead=oldLbead->searchBeadForced(true,vae);
                alpha=NewRbead->chooseTheBead(SumI,vae+1,oldLbead);

                if(alpha==nullptr)return;

                zeta=alpha->searchBead(false,vae);
                if (zeta==nullptr||zeta==Rbead||zeta==Lbead)return;


                    NewRbead->chooseTheBead(SumZ,vae+1,zeta);
                    Rbead=alpha;
                    if(Lbead->swap(zeta,SumI,SumZ,0,true))
                    {
                        Lbead=zeta;
                    }
                    Rbead=oldRbead;


            }
            else
            {
                const auto NewLbead=oldRbead->searchBeadForced(false,vae);
                if(NewLbead==nullptr)return;

                alpha=NewLbead->chooseTheBead(SumI,vae+1,oldRbead);

                if(alpha==nullptr)return;


                zeta=alpha->searchBead(true,vae);
                if (zeta==nullptr||zeta==Rbead||zeta==Lbead)return;

                NewLbead->chooseTheBead(SumZ,vae+1,zeta);
                Lbead=alpha;
                if(Rbead->swap(zeta,SumI,SumZ,0,false))
                {
                    Rbead=zeta;
                }
                Lbead=oldLbead;

            }



    }
    static bool removeWorm(void)
    {

        const auto wormLenght=CalculateWormLenght();
        if(wormLenght>=MBar)
            return false;




        if(Rbead->deleteToRight(wormLenght,-log(eta)))
        {
            NRemo++;
            Site* ptr=Rbead;
            const Site* markPtr;
            size_t step=0;
            do
            {
                markPtr=ptr->right;
                 Site* newPtr=ptr->left;
                const auto timesl=ptr->TimeSliceOnBead;
                const auto topParticle=&(theParticles->at(timesl).back());
                if(newPtr==topParticle)newPtr=topParticle->left;
                    topParticle->up->down=topParticle->down;
                    topParticle->down->up=topParticle->up;

                    if(topParticle!=ptr&&topParticle->active)
                    {
                        ptr->active=true;
                        if(Rbead==topParticle)Rbead=ptr;
                        if(Lbead==topParticle)Lbead=ptr;
                        topParticle->left->right=ptr;
                        topParticle->right->left=ptr;
                        Site* const oldPtrLeft=ptr->left;
                        Site* const oldPtrRight=ptr->right;
                        ptr->left=topParticle->left;
                        ptr->right=topParticle->right;
                        ptr->pos=topParticle->pos;
                        ptr->oldpos=topParticle->oldpos;
                        topParticle->left=oldPtrLeft;
                        topParticle->right=oldPtrRight;



                    }
                            if(!topParticle->right->active)
                            {
                                topParticle->right->left=topParticle->left;
                                topParticle->left->right=topParticle->right;

                            }

                step++;

                theParticles->at(timesl).pop_back();
                if(ptr==markPtr)break;

                ptr=newPtr;


            }while(1);




            Lbead=nullptr;
            Rbead=nullptr;
            ThereIsAWorm=false;
            return true;


        }
        return false;

    }
    static void insertParticle(void)
    {

        const bool var=giveRanI(1);
        vector<double> x;
        for(size_t k=0;k<d;k++)
        {
            x.push_back(Constants::giveRanD(position::L.x.at(k))-position::L.x.at(k)/2);
        }
        for(size_t i=0;i<NTimeSlices;i++)
        {

            theParticles->at(i).push_back(Site(getNParti(),i,(var)?true:false,position(x),false));


            theParticles->at(i).back().up=&(theParticles->at(i).front());

            if(theParticles->at(i).size()>1)
                theParticles->at(i).back().down=&(theParticles->at(i).back())-1;
            else {
                theParticles->at(i).back().down=&(theParticles->at(i).back());
            }
            theParticles->at(i).back().down->up=&(theParticles->at(i).back());
            theParticles->at(i).back().up->down=&(theParticles->at(i).back());

            if(i)
            {
                theParticles->at(i).back().left=&(theParticles->at(i-1).back());
                theParticles->at(i).back().left->right=&(theParticles->at(i).back());
            }

        }

        theParticles->at(NTimeSlices-1).back().right=&(theParticles->at(0).back());
        theParticles->at(0).back().left=&(theParticles->at(NTimeSlices-1).back());
    }
    static void insertWorm(void)
    {

    insertParticle();



        const size_t posiTimes=giveRanI(NTimeSlices-1) ;

        Rbead=&(theParticles->at(posiTimes).back());
        Lbead=Rbead;
            const auto var2=giveRanI(MBar-2);



            Rbead->active=true;



        double U=0,dU=0;
        Rbead->ChangeInU(false,dU,U);

            if(Rbead->insertToLeft(var2,dU+log(eta)))
            {
                NInsert++;

                TPotential+=U;
                ThereIsAWorm= true;
            }
            else
            {

                removeLastParticle();

            }

    }
    Site& operator=(const Site&& v) noexcept {

        pos = std::move(v.pos);
        oldpos = std::move(v.oldpos);
        ParticleOnBead=std::move(v.ParticleOnBead);
        TimeSliceOnBead=std::move(v.TimeSliceOnBead);
        active=std::move(v.active);
        left= std::move(v.left);
        right=std::move(v.right);
        up=std::move(v.up);
        down=std::move(v.down);
        if(v.UPplane)
        {
            Nparti_UpxNT++;
        }
        NpartixNT++;
        return *this;
      }
    Site(const Site& v):UPplane(v.UPplane),pos(v.pos){
        oldpos = v.oldpos;
        ParticleOnBead=v.ParticleOnBead;
        TimeSliceOnBead=v.TimeSliceOnBead;
        active=v.active;
        left= v.left;
        right=v.right;
        up=v.up;
        down=v.down;
        if(v.UPplane)
        {
            Nparti_UpxNT++;
        }
        NpartixNT++;

    }
    Site& operator=(const Site& v)
    {
        auto tmp = v;
        (*this) = std::move(tmp);
        if(v.UPplane)
        {
            Nparti_UpxNT++;
        }
        NpartixNT++;
        return *this;
    }
 inline static size_t  CalculateWormLenght(void)
 {
     size_t wormL=0;

     const Site* ptr=Rbead;
     while(ptr!=Lbead)
     {
         wormL++;
         ptr=ptr->right;
     }
     return wormL;

 }
 inline static bool cantClose(size_t del)
 {

const Site*  ptr=Lbead->right;
while(del>1)
{

    if(ptr==Rbead)
        return false;
    ptr=ptr->right;
    del--;
}
return true;
 }

 inline static bool NfindHead(size_t del,const bool &isRight)
 {
    if(isRight)
    {
        const Site* ptr=Rbead->right;
        while(del>0)
        {
            if(ptr==Lbead)
                return false;
            ptr=ptr->right;
            del--;
        }
        if(ptr==Lbead)
            return false;
        else {
            return true;
        }

    }
    else
    {
        const Site* ptr=Lbead->left;
        while(del>0)
        {
            if(ptr==Rbead)
                return false;
            ptr=ptr->left;
            del--;
        }
        if(ptr==Rbead)
            return false;
        else {
            return true;
        }


    }
 }
 inline size_t  CalculateLenght(void)const
 {
     size_t wormL=1;

     const Site* ptr=right;
     while(ptr!=this)
     {
         wormL++;
         ptr=ptr->right;
     }
     return wormL;
 }
inline static size_t  CalculateNoWormLenght(void)
 {
     size_t wormL=0;

     const Site* ptr=Lbead;
     while(ptr!=Rbead)
     {
         wormL++;
         ptr=ptr->right;
     }
     return wormL;
 }
     bool OpenWorm(const size_t, const size_t, double dU, const position & start);
     bool CloseWorm(double dU);
     bool Wiggle(double dU);
     bool shiftParticle(double dU, const position& shift)const;

     bool swap(Site * const , const double& SumI, const double& SumZ, double dU, const bool &isRight);



     bool deleteToRight(const size_t step, double dU);
     bool deleteToLeft(const size_t step, double dU);
     static inline double propagator(const position& a, const position& b,const double& ab)
     {
         return exp(-fact1(ab)-(a-b).norm()*fact2(ab));

     }
     static  inline double propagator(const position& a, const position& b,const double& ab, const double& dU )
     {

         return exp(-fact1(ab)-(a-b).norm()*fact2(ab)+dU);

     }
     bool insertToRight(const size_t step, double dU);


     bool insertToLeft(const size_t step, double dU);



    size_t TimeSliceOnBead,ParticleOnBead;
    static double TEnergy,mu,TPotential,eta;
    static position TWindingUp,TWindingDown;
    static double TEnergyVar,TPotentialVar;
    static position TWindingVarUp,TWindingVarDown;


    Site* chooseTheBead(double &SumI, const size_t &vae, const Site * const startBead) const;
    static bool ThereIsAWorm;

    static ifstream RestartPtrConf;
    inline void readPtr(size_t &LeftPa,size_t &LeftTi,size_t & RightPa,size_t & RightTi)const
   {
        if(!RestartPtrConf.is_open())
        {
            cout<<"Error; The file .restartPtr.Conf couldn't be opened"<<endl;
            throw std::exception();
        }
       RestartPtrConf>>LeftPa>>LeftTi>>RightPa>>RightTi;

   }
    void setPtr(void)
    {
        if(restart)
        {
            size_t LeftPa,LeftTi,RightPa,RightTi;
            readPtr(LeftPa,LeftTi,RightPa,RightTi);
            left=&(theParticles->at(LeftTi).at(LeftPa));
            right=&(theParticles->at(RightTi).at(RightPa));

        }
        else
        {
            left=&(theParticles->at((TimeSliceOnBead-1+NTimeSlices)%NTimeSlices).at(ParticleOnBead));
            right=&(theParticles->at((TimeSliceOnBead+1)%NTimeSlices).at(ParticleOnBead));
        }

        up=&(theParticles->at(TimeSliceOnBead).at((ParticleOnBead+1)%getNParti()));
        down=&(theParticles->at(TimeSliceOnBead).at((ParticleOnBead-1+getNParti())%getNParti()));



    }
    inline void totalEnergy(void)const
    {

        TEnergy+=(right->pos-pos).normxy();



        (UPplane)?TWindingUp=TWindingUp+(right->pos-pos):TWindingDown=TWindingDown+(right->pos-pos);

        double U=0;
        position gU=position(0.);
        ThePotential.PairInteraction(this,U,gU);
        U/=2;
        ThePotential.ExternPot(this,U,gU);
        TPotential+=U;


    }

    inline static void MoveWorm(void)
    {        
        const auto del=giveRanI(MBar-1);
        switch ( giveRanI(4) )
        {
        case 0:
        {
            if(NfindHead(del,true))
                Rbead->deleteToRight(del,0);
            break;
        }
        case 1:
        {

            Rbead->insertToLeft(del,0);
            break;
        }
        case 2:
        {

            if(NfindHead(del,false))
                Lbead->deleteToLeft(del,0);
            break;
        }
        case 3:
        {
            Lbead->insertToRight(del,0);
            break;
        }



        }


 }
    static void insertParticle(const bool &UP)
    {


        for(size_t i=0;i<NTimeSlices;i++)
        {

            theParticles->at(i).push_back(Site(getNParti(),i,UP,position(0.),false));


            theParticles->at(i).back().up=&(theParticles->at(i).front());

            if(theParticles->at(i).size()>1)
                theParticles->at(i).back().down=&(theParticles->at(i).back())-1;
            else {
                theParticles->at(i).back().down=&(theParticles->at(i).back());
            }
            theParticles->at(i).back().down->up=&(theParticles->at(i).back());
            theParticles->at(i).back().up->down=&(theParticles->at(i).back());

            if(i)
            {
                theParticles->at(i).back().left=&(theParticles->at(i-1).back());
                theParticles->at(i).back().left->right=&(theParticles->at(i).back());
            }

        }

        theParticles->at(NTimeSlices-1).back().right=&(theParticles->at(0).back());
        theParticles->at(0).back().left=&(theParticles->at(NTimeSlices-1).back());
    }
    const static bool fourthOrder;
inline void ChangeInU(const bool & isRemove, double& dU ,double & U )const
{


    const bool isEven=this->TimeSliceOnBead%2==0;
    position gU=position(0.);
    U=0;
    ThePotential.U(this,U,gU);
    double var;
    if(fourthOrder)
    {
        var=2*tao*U/3+((isEven)?0:2*tao*U/3+tao*tao*tao*landa*2/9*(gU).norm());
    }
    else {
        var=tao*U;
    }



    if(isRemove)
        {
            dU+=var;

             return ;
        }
        if(!isRemove)
        {            
            dU-=var;

             return ;
        }

}


    position pos,oldpos;
    Site* left,* right,* up,* down;
    const static potential ThePotential;
    static size_t NClose,NWiggle,NShift,NOpen,NMove,NSwap,NInsert,NInsertP,NRemo,NRemoP,NCloseP,NWiggleP,NShiftP,NOpenP,NMoveP,NSwapP;
    static int Nparti_UpxNT,NpartixNT;

 static inline void restartRatios(void){
     NClose=1;
     NOpen=1;NMove=1;NSwap=1;NInsert=1;NInsertP=1;NRemo=1;NRemoP=1;NCloseP=1;NOpenP=1;NMoveP=1;NSwapP=1;
 }
    static inline double fact1(const double N=1)
    {

        return log(4*landa*tao*pi*N)*d/2;
    }
    static inline double fact2(const double N=1)
    {
        return 1/(4*landa*tao*N);
    }
    inline static size_t NInactiveLinks(void) {
        if(Rbead!=nullptr&&Lbead!=nullptr)
        {
            const int dis=Rbead->TimeSliceOnBead-Lbead->TimeSliceOnBead;
            if(dis<0)
                return NTimeSlices+ dis;
            if(dis>0)
                return dis;
            if(dis==0)
                return NTimeSlices;
        }
        else
            cout<<"#errorSite.hh"<<Rbead->TimeSliceOnBead<<" "<<Lbead->TimeSliceOnBead<<endl;

    }
    inline static size_t NActiveLinks(void){
        if(Rbead!=nullptr&&Lbead!=nullptr)
        {
            const int dis=Rbead->TimeSliceOnBead-Lbead->TimeSliceOnBead;
            if(dis<0)
                return  -1*dis;
            if(dis>0)
                return NTimeSlices-dis;
            if(dis==0)
                return NTimeSlices;

        }
        else
            cout<<"#errorSite.hh"<<Rbead->TimeSliceOnBead<<" "<<Lbead->TimeSliceOnBead<<endl;

    }


    inline static void removeLastParticle(void)
    {
        size_t step=0;
        Site* ptr;

        while(step<NTimeSlices)
        {
            ptr=&theParticles->at(step).back();
            ptr->up->down=ptr->down;
            ptr->down->up=ptr->up;
            theParticles->at(ptr->TimeSliceOnBead).pop_back();
            step++;
        }

    }
    inline Site* searchBead(const bool & isRight, size_t step)const
    {
        Site* var=nullptr;
        if(isRight)
         var=this->right;
        if(!isRight)
            var=this->left;
        while (step>0) {
            if(var==Rbead||var==Lbead)
                return nullptr;
            if(isRight)
             var=var->right;
            if(!isRight)
                var=var->left;
            step--;
        }

        return var;
    }
    inline Site* searchBeadForced(const bool & isRight, size_t step)const
    {
        Site* var=nullptr;
        if(isRight)
         var=this->right;
        if(!isRight)
            var=this->left;
        while (step>0) {
            if(isRight)
             var=var->right;
            if(!isRight)
                var=var->left;
            step--;
        }

        return var;
    }
    static void printLattice(const string svar)
    {

    static int step=1;

    ofstream PrintConf(("Print" +to_string(step) +     ".conf").c_str());
        Site *var;
        for(size_t j=0;j<getNParti();j++)
         {
             for (size_t i=0;i<NTimeSlices;i++)
             {
                 var=&(theParticles->at(i).at(j));
                 if(var->active)
                     PrintConf<<"1 "<<var->ParticleOnBead<<" "<<var->TimeSliceOnBead<<" "<<var->right->ParticleOnBead<<" "<<var->right->TimeSliceOnBead<<" "<<var->pos<<" "<<var->UPplane<<" "<<svar<<endl;
                 else {
                     PrintConf<<"0 "<<var->ParticleOnBead<<" "<<var->TimeSliceOnBead<<" "<<var->right->ParticleOnBead<<" "<<var->right->TimeSliceOnBead<<" "<<var->pos<<" "<<var->UPplane<<" "<<svar<<endl;
                 }



             }

         }
        step++;
        PrintConf.close();
    }

    bool active;
    static Site* theZeta;
    const bool UPplane;
    static array<vector<Site>,100000>* theParticles;
	private:





};


#endif
