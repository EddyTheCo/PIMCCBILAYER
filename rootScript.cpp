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

void CalStruFact(size_t ini,size_t NT)
{
    gROOT->SetBatch(kTRUE);

        TH1 * M4 =nullptr;
        TH1 * M2 =nullptr;
         TH1 * AveStru =nullptr;

    size_t stp=0;


    v = (TVectorD*)gDirectory->Get("v");

    bool var=1;

    TCanvas*c1 = new TCanvas("c1", "c1", 1200,1000);
    gStyle->SetOptStat(0);


    for(int i=ini;i<=(*v)[0];i++)
        {
        cout<<i<<endl;
            TYP* hist = nullptr;
            hist=(TYP* )gDirectory->Get(("pos" + to_string(i)).c_str());
            if(hist!=nullptr)
            {



#ifdef D3
    TH1 * xyProj=hist->Project3D("xy");
#else
    TH1 * xyProj=(TH1 *)hist->Clone();
#endif
    xyProj->Scale(1./NT);

    static vector<TH1 *> VecMagPlus;
    static vector<TH1 *> VecPhasePlus;
    static vector<TH1 *> VecMagMinus;
    static vector<TH1 *> VecPhaseMinus;
    static TH1* M4Hist=nullptr;
    static TH1* M2Hist=nullptr;
cout<<"xyproj"<<endl;
    for(size_t j=1;j<=xyProj->GetXaxis()->GetNbins();j++)
    {
            cout<<xyProj->GetBinContent(xyProj->GetBin(j))<<" ";

    }
    cout<<endl;
    TVirtualFFT::SetTransform(0);
    TH1 *hmP =0;
    TH1 *hmM =0;

    hmP = xyProj->FFT(hmP, "MAG C2CFORWARD");
    hmP->SetName(("MAG C2CFORWARD"  + to_string(i)).c_str());
    VecMagPlus.push_back(hmP);

    cout<<"hmp"<<endl;
        for(size_t j=1;j<=xyProj->GetXaxis()->GetNbins();j++)
        {
                cout<<hmP->GetBinContent(xyProj->GetBin(j))<<" ";

        }
    cout<<endl;


    hmM = xyProj->FFT(hmM, "MAG C2CBACKWARD");
    hmM->SetName(("MAG C2CBACKWARD"  + to_string(i)).c_str());
    VecMagMinus.push_back(hmM);

    cout<<"hmM"<<endl;
        for(size_t j=1;j<=xyProj->GetXaxis()->GetNbins();j++)
        {
                cout<<hmM->GetBinContent(xyProj->GetBin(j))<<" ";

        }
    cout<<endl;



    TH1 *hpP = 0;
    TH1 *hpM = 0;
    hpP = xyProj->FFT(hpP, "PH C2CFORWARD");
    hpP->SetName(("PH C2CFORWARD"  + to_string(i)).c_str());
    VecPhasePlus.push_back(hpP);
    cout<<"hpP"<<endl;
        for(size_t j=1;j<=xyProj->GetXaxis()->GetNbins();j++)
        {
                cout<<hpP->GetBinContent(xyProj->GetBin(j))<<" ";

        }
    cout<<endl;


    hpM = xyProj->FFT(hpM, "PH C2CBACKWARD");
    hpM->SetName(("PH C2CBACKWARD"  + to_string(i)).c_str());
    VecPhaseMinus.push_back(hpM);
    cout<<"hpM"<<endl;
        for(size_t j=1;j<=xyProj->GetXaxis()->GetNbins();j++)
        {
                cout<<hpM->GetBinContent(xyProj->GetBin(j))<<" ";

        }
    cout<<endl;
delete hist;
    delete xyProj;

    if(!M4Hist)
    {
        M4Hist=(TH1*)hmM->Clone();
        M4Hist->Multiply(hmP);
        M4Hist->Multiply(hmM);
        M4Hist->Multiply(hmP);
        M2Hist=(TH1*)hmP->Clone();
        M2Hist->Multiply(hmP);
    }
    else
    {
        TH1* varM4=(TH1*)hmM->Clone();
        varM4->Multiply(hmP);
        varM4->Multiply(hmM);
        varM4->Multiply(hmP);
        M4Hist->Add(varM4);
        TH1* varM2=(TH1*)hmP->Clone();
        varM2->Multiply(hmP);
        M2Hist->Add(varM2);
        delete varM2;
        delete varM4;
    }
    cout<<"M4Hist"<<endl;
        for(size_t j=1;j<=hpM->GetXaxis()->GetNbins();j++)
        {
                cout<<M4Hist->GetBinContent(hpM->GetBin(j))<<" ";

        }
    cout<<endl;
    cout<<"M2hist"<<endl;
        for(size_t j=1;j<=hpM->GetXaxis()->GetNbins();j++)
        {
                cout<<M2Hist->GetBinContent(hpM->GetBin(j))<<" ";

        }
    cout<<endl;


    if(!M4)
    {
        M4=(TH1*) M4Hist->Clone();
        M2=(TH1*) M2Hist->Clone();
    }
    else
    {
        M4->Add(M4Hist);
        M2->Add(M2Hist);
    }
    static TH1 * M2cos=nullptr;
    static TH1 * M4cos=nullptr;

bool first=true;
    if(VecMagPlus.size()>1)
    {
        if(!M2cos)
        {
            M2cos=(TH1*)VecMagPlus.back()->Clone();

            M4cos=(TH1*)VecMagPlus.back()->Clone();
            M4cos->Multiply(VecMagMinus.back());
        }
            TH1* PhaseSum=(TH1*)VecPhaseMinus.back()->Clone();
            PhaseSum->Add(VecPhasePlus.back());
            TH1* PhaseK=(TH1*)VecPhasePlus.back()->Clone();
            PhaseK->Scale(-1.);
            PhaseSum->Scale(-1.);
                    for(size_t j=0;j<VecMagPlus.size()-1;j++)
                    {


                        TH1* PhaseSum2=(TH1* )VecPhaseMinus.at(j)->Clone();
                                PhaseSum2->Add(VecPhasePlus.at(j));

                        PhaseSum2->Add(PhaseSum);

                        TH1* DiffPase4=PhaseSum2;
                        PhaseK->Add(VecPhasePlus.at(j));
                        TH1* DiffPase2=PhaseK;

                        for(size_t j=1;j<=DiffPase4->GetXaxis()->GetNbins();j++)
                        {
                            for(size_t k=1;k<=DiffPase4->GetYaxis()->GetNbins();k++)
                            {
                                DiffPase2->SetBinContent( DiffPase2->GetBin(j,k),cos(DiffPase2->GetBinContent(DiffPase2->GetBin(j,k))));
                                DiffPase4->SetBinContent( DiffPase4->GetBin(j,k),cos(DiffPase4->GetBinContent(DiffPase4->GetBin(j,k))));
                            }
                        }
                        if(first)
                        {
                            M2cos->Multiply(VecMagPlus.at(j));
                            M2cos->Multiply(DiffPase2);
                            M2cos->Scale(2.0);
                            M4cos->Multiply(VecMagMinus.at(j));
                            M4cos->Multiply(VecMagPlus.at(j));
                            M4cos->Multiply(DiffPase4);
                            M4cos->Scale(2.0);
                            first=false;
                        }
                        else
                        {
                            TH1* M2cosvar=(TH1*)VecMagPlus.back()->Clone();

                            TH1* M4cosvar=(TH1*)VecMagPlus.back()->Clone();
                            M4cosvar->Multiply(VecMagMinus.back());

                            M2cosvar->Multiply(VecMagPlus.at(j));
                            M2cosvar->Multiply(DiffPase2);
                            M2cosvar->Scale(2.0);
                            M4cosvar->Multiply(VecMagMinus.at(j));
                            M4cosvar->Multiply(VecMagPlus.at(j));
                            M4cosvar->Multiply(DiffPase4);
                            M4cosvar->Scale(2.0);
                              M4cos->Add(M4cosvar)      ;
                              M2cos->Add(M2cosvar)      ;
                              delete M2cosvar;
                              delete M4cosvar;
                        }

                    }

                    cout<<"M4cos"<<endl;
                        for(size_t j=1;j<=hpM->GetXaxis()->GetNbins();j++)
                        {
                                cout<<M4cos->GetBinContent(hpM->GetBin(j))<<" ";

                        }
                    cout<<endl;

                    cout<<"M2cos"<<endl;
                        for(size_t j=1;j<=hpM->GetXaxis()->GetNbins();j++)
                        {
                                cout<<M2cos->GetBinContent(hpM->GetBin(j))<<" ";

                        }
                    cout<<endl;

                    M4->Add(M4cos);
                    M2->Add(M2cos);



        }

                stp++;

            }

        }


    M2->Scale(-1./stp);
            for(size_t j=1;j<=M4->GetXaxis()->GetNbins();j++)
            {
                for(size_t k=1;k<=M4->GetYaxis()->GetNbins();k++)
                {
                    M4->SetBinContent( M4->GetBin(j,k),sqrt(M4->GetBinContent(M4->GetBin(j,k))));
                }
            }
    M4->Scale(1./stp);
    M4->Add(M2);
    AveStru=M4;

    AveStru->SetTitle("Static Structure Factor");
    AveStru->GetXaxis()->SetTitle("Kx");
#ifndef D1
    AveStru->GetYaxis()->SetTitle("Ky");
#endif
    AveStru->GetYaxis()->CenterTitle(true);

    AveStru->GetXaxis()->CenterTitle(true);
#ifndef D1
    AveStru->Draw("CONT4Z");
#else
      AveStru->Draw("Hist");
#endif
    c1->Print("AverStruc.png");












}

TYP * Average (int &ini,size_t NT)
{
    gROOT->SetBatch(kTRUE);
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
    CalStruFact(start,NT);
#ifdef D3
        XYProj(hist,0,100,NT*start);
        XZProj(hist,0,100,NT*start);
#endif
}

