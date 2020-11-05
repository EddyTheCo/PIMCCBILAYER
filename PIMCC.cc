#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <chrono>
#include <fstream>
#include<iostream>
#include"lattice.hh"



using namespace std;

void printinput(void);
int main()
{
auto start = chrono::high_resolution_clock::now();

if(ReadFromInput<string>(10)=="restart")    
{
    Constants::readRandom();
    int va=system("cp input .input.start");
}

int va=system("cp input .input.ini"); //Makes a copy of the input file

cout<<"Version 1.0"<<endl;

   auto theLattice=lattice();

    theLattice.setup();



    theLattice.move();

auto finish = std::chrono::high_resolution_clock::now();
chrono::duration<double> elapsed = finish - start;
cout<<" Elapsed time: " << elapsed.count()<<"s"<<endl;

}


