#include <cstdlib>
#include <fstream>
#include<iostream>

using namespace std;

#ifndef Position_h
#define Position_h 1


#include"constants.hh"
using namespace Constants;

class Site;

const extern   size_t d;
const extern bool periodic;
class position
{
	public:

    position(const string str, const Site *const bead);
static  ifstream PosiConf;
static  ifstream RestartConf;
static  ifstream *inFile;


static double foovar;
inline static double foo(void)
{
    size_t N;
    string str;
    *inFile>>N>>str;
    return 0.;
}
position()
{

}
static const double volumen;
inline static double getVolumen(void)
{
    auto volu=L.x.at(0);
    for(size_t i=1;i<d;i++)
        volu*=L.x.at(i);
    return volu;
}
position& operator=(const position&& v) noexcept {

    x = std::move(v.x);

    return *this;
  }
position(const position& v){
    x=v.x;
}
position& operator=(const position& v) {

  const auto tmp = v;
  (*this) = std::move(tmp);
  return *this;


}

position(const vector<double> &a)
{
    size_t var=a.size();
    for(size_t i=0;i<d;i++)
    {
        if(var>0)
            x.push_back(a.at(i));
        else
            x.push_back(a.at(0));
        var--;
    }
}
position(const double &a);
position(const position& prevPos, const double &var)
{
    for(size_t i=0;i<d;i++)
    {
        x.push_back(prevPos.x.at(i)+Constants::giveRanDNormal(0,sqrt(var)));
    }
}
position(const string str);

double TheX  (void) const {return x.at(0);}
double TheY  (void) const {return (d==2)?x.at(1):0;}
double TheZ  (void) const {return (d==3)?x.at(2):0;}

friend std::ostream & operator << (std::ostream &out, const position & obj)
{
    for(size_t i=0;i<d;i++)
        out<<obj.x.at(i)<<" ";
    return out;
}
friend std::istream & operator >> (std::istream &in,  position & obj)
{
    for(size_t i=0;i<d;i++)
    {
        double var;
        in>>var;
        obj.x.push_back(var);
    }
    return in;
}

bool operator<(const position& rhs)const
{
    return this->norm()<rhs.norm();
}
bool operator>(const position& rhs)const
{
    return rhs<(*this);
}

position operator*(const double &a)const
{
    position var=position();
    for(size_t i=0;i<d;i++)
    {
        var.x.push_back(x.at(i)*a);
    }
    return var;
}
inline double perio(const size_t &i)const
{
    if(periodic)
    {
            auto pro=x.at(i);
            pro/=L.x.at(i);
            pro-=floor(pro+0.5);
           return pro*L.x.at(i);
    }
    return x.at(i);
}
position operator/(const double a)const
{
    return (*this)*(1/a);
}
inline double norm(void) const
{
    double var=0;
    for(size_t i=0;i<d;i++)
        var+=x.at(i)*x.at(i);
    return var;
}
inline double normxy(void) const
{
    double var=0;
    for(size_t i=0;i<2;i++)
        var+=x.at(i)*x.at(i);
    return var;
}
inline position getL(void) const {return L;}
double operator*(const position & rhs) const
{
    double pro=0;
    for(size_t i=0;i<d;i++)
    {
        pro+=x.at(i) * rhs.x.at(i);
    }
    return pro;
}
const position operator+(const position & rhs) const
{
    vector<double> var;
    for(size_t i=0;i<d;i++)
    {
        double pro=x.at(i) + rhs.x.at(i);

            var.push_back(pro);

    }

    return position(var);
}

const position operator-(const position & rhs) const
{
    vector<double> var;
    for(size_t i=0;i<d;i++)
    {
        double pro=x.at(i) - rhs.x.at(i);

       if(periodic)
        {

                pro/=L.x.at(i);
                pro-=floor(pro+0.5);
                var.push_back(pro*L.x.at(i));

        }
        else
       {
           var.push_back(pro);
        }
    }
    return position(var);
}


    vector<double> x;

    static  const position L;
	private:


};



#endif
