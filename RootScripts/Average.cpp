#include <iostream>
#include <fstream>
#include <vector>
TFile *MyFile = new TFile("RootFile.root","UPDATE");
TVectorD *v=nullptr;

void Average ()
{
gROOT->SetBatch(kTRUE);
        TH1 * Ave =nullptr;


    size_t stp=0;


    v = (TVectorD*)gDirectory->Get("v");
    TH1* hist = nullptr;
    bool var=1;

    TCanvas*c1 = new TCanvas("c1", "c1", 1200,1000);
    gStyle->SetOptStat(0);


    for(int i=1;i<=(((*v)[0]<1000)?((*v)[0]):(999));i++)
        {

            cout<<("pos" + to_string(i)).c_str()<<endl;
            hist=(TH1* )gDirectory->Get(("pos" + to_string(i)).c_str());
            if(hist!=nullptr)
            {

                if(var)
                {
                    Ave=(TH1*)hist->Clone();

                    var=0;

                }
                else
                {

                   Ave->Add(hist);


                }

            }
            delete hist;
            stp=i;

        }

    hist=(TH1* )gDirectory->Get(("pos" + to_string(stp)).c_str());
    hist->GetXaxis()->SetTitle("X");
    hist->GetYaxis()->SetTitle("Y");
    hist->GetZaxis()->SetTitle("Z");
    hist->GetYaxis()->CenterTitle(true);
    hist->GetZaxis()->CenterTitle(true);
    hist->GetXaxis()->CenterTitle(true);
    hist->Draw();
    c1->Print("AverageLast.png");


    Ave->GetXaxis()->SetTitle("X");
    Ave->GetYaxis()->SetTitle("Y");
    Ave ->GetZaxis()->SetTitle("Z");
    Ave->GetYaxis()->CenterTitle(true);
    Ave->GetZaxis()->CenterTitle(true);
    Ave->GetXaxis()->CenterTitle(true);
    Ave->Draw();
    c1->Print("Average.png");
    Ave->Scale(1.0/stp);
    Ave->Write("Average",TObject::kOverwrite);
    gDirectory->Write("", TObject::kOverwrite);





}
