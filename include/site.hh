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
extern size_t NpartiUp;

class Site
{
	public:
    Site();
    ~Site(){}
    Site(const size_t i,const size_t j);
    Site(const size_t i, const size_t j, const string var);
    static Site* Lbead;
    static Site* Rbead;
    friend ostream & operator << (ostream &out, const Site & obj)
    {
                if(obj.active)
                out<<obj.pos;
        return out;
    }
    void setLattice(array<vector<Site>,10000>* const particles)const
    {
        theParticles=particles;
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
        return *this;
      }
    Site(const Site& v){
        pos = v.pos;
        oldpos = v.oldpos;
        ParticleOnBead=v.ParticleOnBead;
        TimeSliceOnBead=v.TimeSliceOnBead;
        active=v.active;
        left= v.left;
        right=v.right;
        up=v.up;
        down=v.down;
    }
    Site& operator=(const Site& v)
    {
        auto tmp = v;
        (*this) = std::move(tmp);
        return *this;
    }
 inline size_t  CalculateWormLenght(void)const
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
 inline bool cantClose(size_t del)const
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

 inline bool NfindHead(size_t del,const bool &isRight)const
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
inline size_t  CalculateNoWormLenght(void)const
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
     bool shiftParticle(double dU, const position& shift, const Site * const &str)const;
     void PrepareSwap(void)const;
     bool swap(Site * const , const double& SumI, const double& SumZ, double dU, const bool &isRight);



     bool deleteToRight(const size_t step, double dU);
     bool deleteToLeft(const size_t step, double dU);
     inline double propagator(const position& a, const position& b,const double& ab)const
     {
         return exp(-fact1(ab)-(a-b).norm()*fact2(ab));

     }
     inline double propagator(const position& a, const position& b,const double& ab, const double& dU )const
     {
         return exp(-fact1(ab)-(a-b).norm()*fact2(ab)+dU);

     }
     bool insertToRight(const size_t step, double dU);
     void insertWorm(void);
     void removeWorm(void);
     bool insertToLeft(const size_t step, double dU);



    size_t TimeSliceOnBead,ParticleOnBead;
    static double TEnergy,mu,TPotential,eta;
    static position TWindingUp,TWindingDown;
    static double TEnergyVar,TPotentialVar;
    static position TWindingVar;


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

        up=&(theParticles->at(TimeSliceOnBead).at((ParticleOnBead+1)%NParti_));
        down=&(theParticles->at(TimeSliceOnBead).at((ParticleOnBead-1+NParti_)%NParti_));
    }
    inline void totalEnergy(void)const
    {

        TEnergy+=(right->pos-pos).norm();
        TWinding=TWinding+(right->pos-pos);
        double U=0;
        position gU=position(0.);
        ThePotential.PairInteraction(this,U,gU);
        U/=2;
        ThePotential.ExternPot(this,U,gU);
        TPotential+=U;


    }

    inline void MoveWorm(void)const
    {        
        const auto del=giveRanI(MBar-1);
        switch ( giveRanI(4) ) {
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
    static size_t NParti_,NClose,NOpen,NMove,NSwap,NInsert,NInsertP,NRemo,NRemoP,NCloseP,NOpenP,NMoveP,NSwapP;


 inline void restartRatios(void)const{
     NClose=1;
     NOpen=1;NMove=1;NSwap=1;NInsert=1;NInsertP=1;NRemo=1;NRemoP=1;NCloseP=1;NOpenP=1;NMoveP=1;NSwapP=1;
 }
    inline double fact1(const double N=1)const
    {
        return log(4*landa*tao*pi*N)*d/2;
    }
    inline double fact2(const double N=1) const
    {
        return 1/(4*landa*tao*N);
    }
    inline size_t NInactiveLinks(void) const{
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
    inline size_t NActiveLinks(void)const{
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
    void insertParticle(void)const;
    inline void removeLastParticle(void)const
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
    void printLattice(void)const
    {


        Site var;
        for(size_t j=0;j<NParti_;j++)
         {
             for (size_t i=0;i<NTimeSlices;i++)
             {
                 var=theParticles->at(i).at(j);
                 if(var.active)
                     cout<<var.pos;
                 else {
                     cout<<"#"<<var.pos<<" ";
                 }
                 cout<<" ";
                 //cout<<var.right<<"  ";

             }
             cout<<endl;

         }
    }

    bool active;
    static Site* theZeta;
	private:

    static array<vector<Site>,10000>* theParticles;



};


#endif
