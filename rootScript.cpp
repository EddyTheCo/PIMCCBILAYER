#include <iostream>
#include <fstream>

using namespace std;


#define D1
#define TYP TH1D

TVectorD *v=nullptr;
TFile *MyFile = new TFile("RootFile.root","UPDATE");

TH1 * StrucFact(TH1 * hist,double &val, int & binMax)
{
    hist->Scale(1./hist->GetEntries());
    TVirtualFFT::SetTransform(0);
    TH1 *hm =0;

    hm = hist->FFT(hm, "MAG R2C M");
    hm->Multiply(hm);
    hm->SetBinContent(1,1,0);
    binMax=hm->GetMaximumBin();
    val=hm->GetBinContent(binMax);
    return hm;

}

TH1 * CalStruFact(TYP * hist)
{
#ifdef D3
    TH1 * xyProj=hist->Project3D("xy");
#else
    TH1 * xyProj=(TH1 *)hist->Clone();
#endif

    ofstream  theMaxXY("theMaxXY",std::ofstream::out | std::ofstream::app);
    double val;
    int binMax;
    TH1 * h1=StrucFact(xyProj,val,binMax);
    theMaxXY<<val<<" "<<binMax<<endl;
    return h1;

}

TYP * Average (int &ini)
{
    gROOT->SetBatch(kTRUE);
        TYP * Ave =nullptr;
        TH1 * AveStru =nullptr;

    size_t stp=0;
    int bin;
    int binSum=0;

    v = (TVectorD*)gDirectory->Get("v");

    bool var=1;

    TCanvas*c1 = new TCanvas("c1", "c1", 1200,1000);
    gStyle->SetOptStat(0);


    for(int i=ini;i<=(*v)[0];i++)
        {
            TYP* hist = nullptr;
            hist=(TYP* )gDirectory->Get(("pos" + to_string(i)).c_str());
            if(hist!=nullptr)
            {   

                if(var)
                {   
                    Ave=(TYP*)hist->Clone();
                    AveStru=(TH1*)(CalStruFact(hist))->Clone();
                    var=0;
			stp++;	
                }
                else
                {  
                    AveStru->Add(CalStruFact(hist));
                   Ave->Add(hist);
			stp++;
                }
                delete gROOT->FindObject("out_MAG R2C M");
            }
            delete hist;
        }
	ini=stp;
    Ave->GetXaxis()->SetTitle("X");
    Ave->GetYaxis()->SetTitle("Y");
    Ave ->GetZaxis()->SetTitle("Z");
    Ave->GetYaxis()->CenterTitle(true);
    Ave->GetZaxis()->CenterTitle(true);
    Ave->GetXaxis()->CenterTitle(true);
    Ave->Draw();
    c1->Print("Avera.png");
    Ave->Write("Ave");
    gDirectory->Write("", TObject::kOverwrite);

    TCanvas* c2 = new TCanvas("c2", "c2", 1300,1000);
    gStyle->SetOptStat(0);

    AveStru->Scale(1.0/stp);

    AveStru->SetTitle("Structure Factor");
    AveStru->GetXaxis()->SetTitle("Kx");
    AveStru->GetYaxis()->SetTitle("Ky");

    AveStru->GetYaxis()->CenterTitle(true);

    AveStru->GetXaxis()->CenterTitle(true);
    AveStru->Draw("CONT4Z");
    c1->Print("AverStruc.png");

    return Ave;


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
        TVectorD * v = (TVectorD*)gDirectory->Get("v");
       int start=starte;
	if(start==0)start=((*v)[0]);

         TYP * hist=Average (start);
#ifdef D3
        XYProj(hist,0,100,NT*start);
        XZProj(hist,0,100,NT*start);
#endif
}

