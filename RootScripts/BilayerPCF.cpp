#include <iostream>
#include <fstream>
#include <vector>
TFile *MyFile = new TFile("RootFile.root","READ");
TVectorD *v=nullptr;

void BilayerPCF(size_t Npart, size_t NTimeSlices,double Rangextop=2,double Rangeytop=2,double Lx=1,double Ly=1,const char *text="")
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


    ofstream PCFOUT("PCF.dat");
    PCFOUT<<"r      PCFUP     PCFUP_ERR       PCFMIX      PCFMIX_ERR"<<endl;

    for(size_t j=1;j<=PCFUp->GetXaxis()->GetNbins();j++)
    {
        PCFOUT<<PCFUp->GetXaxis()->GetBinCenter(j)<<" "<<PCFUp->GetBinContent(j)/PCFUp->GetXaxis()->GetBinCenter(j)
             <<" "<<PCFUp->GetBinContent(j)/PCFUp->GetXaxis()->GetBinCenter(j)*sqrt(PCFUp->GetBinError(j)*PCFUp->GetBinError(j)/
                                            PCFUp->GetBinContent(j)/PCFUp->GetBinContent(j) +
                                            PCFUp->GetXaxis()->GetBinWidth(j)*PCFUp->GetXaxis()->GetBinWidth(j)/
                                            PCFUp->GetXaxis()->GetBinCenter(j)/PCFUp->GetXaxis()->GetBinCenter(j))
            <<" "<<PCFMix->GetBinContent(j)/PCFMix->GetXaxis()->GetBinCenter(j)
                             <<" "<<PCFMix->GetBinContent(j)/PCFMix->GetXaxis()->GetBinCenter(j)*
                               sqrt(PCFMix->GetBinError(j)*PCFMix->GetBinError(j)/
                                                            PCFMix->GetBinContent(j)/PCFMix->GetBinContent(j) +
                                                            PCFMix->GetXaxis()->GetBinWidth(j)*PCFMix->GetXaxis()->GetBinWidth(j)/
                                                            PCFMix->GetXaxis()->GetBinCenter(j)/PCFMix->GetXaxis()->GetBinCenter(j))<<endl;

            PCFUp->SetBinContent(j,PCFUp->GetBinContent(j)/PCFUp->GetXaxis()->GetBinCenter(j));
            PCFMix->SetBinContent(j,PCFMix->GetBinContent(j)/PCFMix->GetXaxis()->GetBinCenter(j));



    }
PCFOUT.close();


    PCFUp->GetXaxis()->SetTitle("r");
    PCFUp->GetYaxis()->SetTitle("g(r)");
    PCFUp->GetXaxis()->SetRangeUser(0.,Rangextop);
    PCFUp->GetYaxis()->SetRangeUser(0.,Rangeytop);


    PCFUp->GetYaxis()->CenterTitle(true);
    PCFUp->GetXaxis()->CenterTitle(true);
     PCFUp->SetMarkerStyle(kFullCircle);
     PCFUp->SetMarkerSize(1);
     PCFUp->SetMarkerColor(kBlack);

     PCFMix->SetMarkerStyle(kFullCircle);
     PCFMix->SetMarkerSize(1);
     PCFMix->SetMarkerColor(kRed);

     TLegend *leg=new TLegend(0.7,0.1,0.89,0.29);
        leg->SetFillColor(0);

       leg->AddEntry(PCFUp,"g_{#alpha#alpha}","P");
       leg->AddEntry(PCFMix,"g_{#alpha#beta}","P");
        leg->SetTextFont(42);

        leg->SetBorderSize(0);


        leg->SetTextSize(0.04);




     TCanvas*c1 = new TCanvas("c1", "c1", 1500,800);
     PCFUp->Draw("PLC ");
     PCFMix->Draw("PLC SAME");
     leg->Draw();



     TLatex latex;
     latex.SetTextSize(0.025);
     latex.SetTextAlign(13);  //align at top
     latex.DrawLatex(0.1,1.85,text);



     c1->Print("PCFUp.png");
    /*



    //PCFUp->Add(PCFDown);

    //PCFDown->Scale(Lx*Ly/Nsamples/2.0/pi/(PCFUp->GetBinWidth(1)*(Npart/2-1)*NTimeSlices*(Npart/2)))

*/
}
