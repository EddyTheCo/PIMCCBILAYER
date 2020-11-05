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
void BilayerPCF(size_t Npart, size_t NTimeSlices,double Rangetop,double Lx,double Ly)
{
    gROOT->SetBatch(kTRUE);
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
     PCFUp->SetMarkerStyle(23);
     PCFUp->SetMarkerSize(1);
     PCFUp->SetMarkerColor(2);

     PCFMix->SetMarkerStyle(kFullTriangleUp);
     PCFMix->SetMarkerStyle(23);
     PCFMix->SetMarkerSize(1);
     PCFMix->SetMarkerColor(3);

     TLegend *leg=new TLegend(0.7,0.7,0.89,0.89);
        leg->SetFillColor(0);

       leg->AddEntry(PCFUp,"g_{#alpha#alpha}","P");
       leg->AddEntry(PCFMix,"g_{#alpha#beta}","P");
        leg->SetTextFont(42);

        leg->SetBorderSize(0);


        leg->SetTextSize(0.04);




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

void Average (int ini,bool SaveFourier)
{
gROOT->SetBatch(kTRUE);
        TYP * Ave =nullptr;


    size_t stp=0;


    v = (TVectorD*)gDirectory->Get("v");

    bool var=1;

    TCanvas*c1 = new TCanvas("c1", "c1", 1200,1000);
    gStyle->SetOptStat(0);


    for(int i=ini;i<=(((*v)[0]<1000)?((*v)[0]):(999));i++)
        {
            TYP* hist = nullptr;
            cout<<("pos" + to_string(i)).c_str()<<endl;
            hist=(TYP* )gDirectory->Get(("pos" + to_string(i)).c_str());
            if(hist!=nullptr)
            {   

                if(var)
                {   
                    Ave=(TYP*)hist->Clone();
                    if(SaveFourier)
                    {
                        StrucFact(hist,i);
                    }
                    var=0;
			stp++;	
                }
                else
                {  

                   Ave->Add(hist);
                   if(SaveFourier)
                   {
                       StrucFact(hist,i);
                   }
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
    Ave->Scale(1.0/stp);
    Ave->Write("Ave",TObject::kOverwrite);
    gDirectory->Write("", TObject::kOverwrite);





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

