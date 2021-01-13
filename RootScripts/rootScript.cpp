#include <iostream>
#include <fstream>
#include <vector>
using namespace std;


#define D3
#define TYP TH3D

TVectorD *v=nullptr;
TFile *MyFile = new TFile("RootFile.root","UPDATE");



void StrucFact(TYP * hist,size_t step)
{

#ifdef D3
    TH1 * xyProj=hist->Project3D("xy");
#else
    TH1 * xyProj=(TH1 *)hist->Clone();
#endif

    xyProj->Scale(1./xyProj->GetEntries());
    TVirtualFFT::SetTransform(0);
    TH1 *hm =0,*hp=0;

    hm = xyProj->FFT(hm, "MAG");
    hp = xyProj->FFT(hp, "PH");

    hm->Write(("hm" + to_string(step)).c_str(),TObject::kOverwrite);
    hp->Write(("hp" + to_string(step)).c_str(),TObject::kOverwrite);

}

void CalcStrucFact(void)
{
    v = (TVectorD*)gDirectory->Get("v");
    cout<<"v="<<(*v)[0]<<endl;
TH1* hmi = nullptr,*hmj=nullptr,*hpi=nullptr,*hpj=nullptr;
TH1* sum=nullptr;
    for(int i=1;i<(((*v)[0]<1000)?((*v)[0]):(999));i++)
        {

            hmi=(TH1* )gDirectory->Get(("hm" + to_string(i)).c_str());
            if(!hmi)continue;
            hpi=(TH1* )gDirectory->Get(("hp" + to_string(i)).c_str());
            for(int j=i+1;j<=(((*v)[0]<1000)?((*v)[0]):(999));j++)
                {
                    cout<<("i=" + to_string(i)).c_str()<<(" j=" + to_string(j)).c_str()<<endl;
                     hmj=(TH1* )gDirectory->Get(("hm" + to_string(j)).c_str());
                     if(!hmj)continue;
                     hpj=(TH1* )gDirectory->Get(("hp" + to_string(j)).c_str());
                    hpi->Scale(-1.0);
                    hpj->Add(hpi);
                    for(size_t a=1;a<=hpj->GetXaxis()->GetNbins();a++)
                    {
                        for(size_t b=1;b<=hpj->GetYaxis()->GetNbins();b++)
                        {
                            hpj->SetBinContent( a,b,cos(hpj->GetBinContent(a,b)));
                        }

                    }

                    hpj->Multiply(hmi);
                    hpj->Multiply(hmj);
                    hpj->Scale(-2.0);
                     if(sum)
                    {
                        sum->Add(hpj);
                    }
                    else {
                        sum=(TH1*)hpj->Clone();
                    }


                     delete hmj;
                     delete hpj;
                }

            delete hmi;
            delete hpi;
        }
    TCanvas*c1 = new TCanvas("c1", "c1", 1200,1000);
    sum->Draw("SURF");
    TH1D * Sq=((TH2D*)sum)->ProjectionX();
    Sq->Scale(1./hpj->GetYaxis()->GetNbins());


    gStyle->SetOptStat(0);
    Sq->GetXaxis()->SetTitle("q");
    Sq->GetYaxis()->SetTitle("S(q)");

    Sq->GetYaxis()->CenterTitle(true);

    Sq->GetXaxis()->CenterTitle(true);

    Sq->SetMarkerStyle(kFullCircle);
    Sq->SetMarkerStyle(23);
    Sq->SetMarkerSize(1);
    Sq->SetMarkerColor(2);
    Sq->SetTitle("");
    //Sq->Draw("PLC PMC");
    c1->Print("Sq.png");



}


#ifdef D3
TH1 * XZProj (TYP * hist,int center, int num,int NT)
{
        gROOT->SetBatch(kTRUE);
         TH1 * xzProj =nullptr;


    int nbi=hist->GetNbinsY();
    hist->GetYaxis()->SetRange(nbi/2-num+center,nbi/2+num+center);
    xzProj=hist->Project3D("xz");
xzProj->Scale(1.0/NT);
TCanvas*c1 = new TCanvas();
    gStyle->SetOptStat(0);
    xzProj->Draw("CONT4Z");
xzProj->GetXaxis()->SetTitle("z");
    xzProj->GetYaxis()->CenterTitle(true);
    xzProj->GetYaxis()->SetTitle("x");
   xzProj->GetXaxis()->CenterTitle(true);
    xzProj->SetTitle("");
    c1->Print("xzProj.png");
    return xzProj;
}
TH1 * XYProj (TYP * hist,int center, int num, int NT)
{
        gROOT->SetBatch(kTRUE);
    TH1 * xyProj =nullptr;
    int nbi=hist->GetNbinsZ();
    hist->GetZaxis()->SetRange(nbi/2-num+center,nbi/2+num+center);
    xyProj=hist->Project3D("xy");
 xyProj->Scale(1.0/NT);
    TCanvas*c1 = new TCanvas();
    gStyle->SetOptStat(0);
    xyProj->Draw("CONT4Z");
    xyProj->GetYaxis()->SetTitle("X");
    xyProj->GetYaxis()->CenterTitle(true);
    xyProj->GetXaxis()->SetTitle("Y");
   xyProj->GetXaxis()->CenterTitle(true);
    xyProj->SetTitle("");
    c1->Print("xyProj.png");
    return xyProj;
}
#endif
void rootScript (int NT, int starte)
{
        gStyle->SetPalette(1);

}

