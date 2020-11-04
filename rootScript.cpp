#include <iostream>
#include <fstream>
#include <vector>
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
void BilayerPCF(size_t Npart, size_t NTimeSlices,double Rangetop,double Lx,double Ly)
{
    //gROOT->SetBatch(kTRUE);
    gStyle->SetOptStat(0);
    TH1* PCFUp = (TH1D*)gDirectory->Get("PCFUp");
    TH1* PCFDown = (TH1D*)gDirectory->Get("PCFDown");
    TH1* PCFMix = (TH1D*)gDirectory->Get("PCFMix");

    double pi=3.141592654;
    v = (TVectorD*)gDirectory->Get("v");
    size_t Nsamples=(*v)[0];
    PCFUp->Add(PCFDown);
    PCFUp->Scale(Lx*Ly/Nsamples/4/pi/(PCFUp->GetBinWidth(1)*(Npart/2-1)*NTimeSlices*(Npart/2)));
    PCFMix->Scale(Lx*Ly/Nsamples/4/pi/(PCFMix->GetBinWidth(1)*(Npart/2)*NTimeSlices*(Npart/2)));




    for(size_t j=1;j<=PCFUp->GetXaxis()->GetNbins();j++)
    {

            PCFUp->SetBinContent(j,PCFUp->GetBinContent(j)/PCFUp->GetXaxis()->GetBinCenter(j));
            PCFMix->SetBinContent(j,PCFMix->GetBinContent(j)/PCFMix->GetXaxis()->GetBinCenter(j));
    }



    PCFUp->GetXaxis()->SetTitle("r");
    PCFUp->GetYaxis()->SetTitle("g(r)");
    PCFUp->GetXaxis()->SetRangeUser(0.,Rangetop);


    PCFUp->GetYaxis()->CenterTitle(true);
    PCFUp->GetXaxis()->CenterTitle(true);
     PCFUp->SetMarkerStyle(kFullCircle);
     PCFUp->SetMarkerStyle(21);
     PCFUp->SetMarkerSize(0.8);
     PCFUp->SetMarkerColor(2);

     PCFMix->SetMarkerStyle(kFullTriangleUp);
     PCFMix->SetMarkerStyle(18);
     PCFMix->SetMarkerSize(0.8);
     PCFMix->SetMarkerColor(3);

     TLegend *leg=new TLegend(0.8,0.8,0.9,0.9);
        leg->SetFillColor(0);

       leg->AddEntry(PCFUp,"g_#alpha#alpha","P");
       leg->AddEntry(PCFMix,"g_#alpha#beta","P");
        leg->SetTextFont(132);



     TCanvas*c1 = new TCanvas("c1", "c1", 1200,1000);
     PCFUp->Draw("PLC PMC");
     PCFMix->Draw("PLC PMC SAME");
     leg->Draw();

     c1->Print("PCFUp.png");
    /*



    //PCFUp->Add(PCFDown);

    //PCFDown->Scale(Lx*Ly/Nsamples/2.0/pi/(PCFUp->GetBinWidth(1)*(Npart/2-1)*NTimeSlices*(Npart/2)))

*/
}


TYP * Average (int &ini,size_t NT)
{

        TYP * Ave =nullptr;
        TH1 * M4 =nullptr;
        TH1 * M2 =nullptr;
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

                    var=0;
			stp++;	
                }
                else
                {  

                   Ave->Add(hist);
			stp++;
                }

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

        // TYP * hist=Average (start,NT);
    //CalStruFact(start,NT);
#ifdef D3
        XYProj(hist,0,100,NT*start);
        XZProj(hist,0,100,NT*start);
#endif
}

