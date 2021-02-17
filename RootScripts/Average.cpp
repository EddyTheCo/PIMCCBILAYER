#include <iostream>
#include <fstream>
#include <vector>
TFile *MyFile = new TFile("RootFile.root","UPDATE");
TVectorD *v=nullptr;


void XZProj (TH1 * hist,string Title)
{
        gROOT->SetBatch(kTRUE);
         TH1 * xzProj =nullptr;
        xzProj=((TH3D*)hist)->Project3D("xz");

    TCanvas*c1 = new TCanvas();
    gStyle->SetOptStat(0);
    xzProj->Draw("CONT4Z");
    xzProj->GetXaxis()->SetTitle("z");
    xzProj->GetYaxis()->CenterTitle(true);
    xzProj->GetYaxis()->SetTitle("x");
   xzProj->GetXaxis()->CenterTitle(true);
    xzProj->SetTitle("");
    c1->Print(("xzProj"+Title+".png").c_str());
    return xzProj;
}
void XYProj (TH1 * hist, string Title)
{
        gROOT->SetBatch(kTRUE);
    TH1 * xyProj =nullptr;
    if(hist->GetZaxis()->GetNbins()>1)
    {
        xyProj=((TH3D*)hist)->Project3D("xy");
    }
    else
    {
        if(hist->GetYaxis()->GetNbins()>1)
        {
            xyProj=(TH2D*)hist->Clone();
        }
        else
        {
            return;

        }
    }


    TCanvas*c1 = new TCanvas();
    gStyle->SetOptStat(0);
    xyProj->Draw("CONT4Z");
    xyProj->GetYaxis()->SetTitle("X");
    xyProj->GetYaxis()->CenterTitle(true);
    xyProj->GetXaxis()->SetTitle("Y");
   xyProj->GetXaxis()->CenterTitle(true);
    xyProj->SetTitle("");
    c1->Print(("xyProj"+Title+".png").c_str());
    return xyProj;
}

void Average (string pp)
{
gROOT->SetBatch(kTRUE);
        TH1 * Ave =nullptr;


    size_t stp=0;


    v = (TVectorD*)gDirectory->Get("v");
    cout<<"v="<<(*v)[0]<<endl;
    TH1* hist = nullptr;
    bool var=1;

    TCanvas* c1 = new TCanvas("c1", "c1", 1200,1000);
    gStyle->SetOptStat(0);


    for(int i=1;i<=(((*v)[0]<1000)?((*v)[0]):(999));i++)
        {


            hist=(TH1* )gDirectory->Get((pp + to_string(i)).c_str());
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

    hist=(TH1* )gDirectory->Get((pp + to_string(stp)).c_str());
    hist->GetXaxis()->SetTitle("X");
    hist->GetYaxis()->SetTitle("Y");
    hist->GetZaxis()->SetTitle("Z");
    hist->GetYaxis()->CenterTitle(true);
    hist->GetZaxis()->CenterTitle(true);
    hist->GetXaxis()->CenterTitle(true);
    hist->Draw();
    c1->Print(("AverageLast" + pp + ".png").c_str());

    XYProj(hist,"Last"+pp);
    if(hist->GetZaxis()->GetNbins()>1)
    {
        XZProj(hist,"Last"+pp);
    }

    TCanvas* c2 = new TCanvas("c2", "c2", 1200,1000);
    cout<<"Last NEntries="<<hist->GetEntries()<<endl;
    Ave->GetXaxis()->SetTitle("X");
    Ave->GetYaxis()->SetTitle("Y");
    Ave ->GetZaxis()->SetTitle("Z");
    Ave->GetYaxis()->CenterTitle(true);
    Ave->GetZaxis()->CenterTitle(true);
    Ave->GetXaxis()->CenterTitle(true);
    Ave->Draw();
    cout<<"Full NEntries="<<Ave->GetEntries()<<endl;

    c2->Print(("Average"+pp+".png").c_str());

    XYProj(Ave,"Average"+pp);
    if(hist->GetZaxis()->GetNbins()>1)
    {
            XZProj(Ave,"Average"+pp);
    }





}
