#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <fstream>
#include<iostream>
#include"lattice.hh"



using namespace std;

bool CheckRestartFiles(void)
{
    const size_t D=ReadFromInput<size_t>(6);
    const double Tao=ReadFromInput<double>(2);
    const double Beta=ReadFromInput<double>(3);
    const size_t NTImeSlices=Beta/Tao;
    size_t NPArti_=ReadFromInput<size_t>(1,".restart.conf");

    system("tail -n1 .restart.conf| tr -cd [:space:]| wc|awk '{print $3 \"   NPArti_*NTImeSlices*D\"}'>.GoodTogo");
    system("tail -n1 .restartPtr.conf| tr -cd [:space:]| wc|awk '{print $3  \"   NPArti_*NTImeSlices*4+1\"}'>>.GoodTogo");

    if(ReadFromInput<size_t>(1,".GoodTogo")!=NPArti_*NTImeSlices*D)
        return false;
    if(ReadFromInput<size_t>(2,".GoodTogo")!=NPArti_*NTImeSlices*4)
        return false;

    return true;
}


int main()
{
auto start = chrono::high_resolution_clock::now();
bool GoodToGo=true;
if(ReadFromInput<string>(10)=="restart")
{
Constants::readRandom();
GoodToGo=CheckRestartFiles();
}
else
{
int va=system("cp input .input.start"); //Makes a copy of the input file
}
int va=system("cp input .input.ini"); //Makes a copy of the input file


if(GoodToGo)
{
   auto theLattice=lattice();

    theLattice.setup();


if(size_t Warmup=ReadFromInput<int>(22))
{
    cout<<"Starting The WarmingUp"<<endl;
    theLattice.Warm();
    cout<<"Finishing The WarmingUp"<<endl;
}
    theLattice.move();
}
else
{
    cout<<"ERROR #4<< The restart files configurations are corrupted"<<endl;
}
auto finish = std::chrono::high_resolution_clock::now();
chrono::duration<double> elapsed = finish - start;
cout<<" Elapsed time: " << elapsed.count()<<"s"<<endl;

}


