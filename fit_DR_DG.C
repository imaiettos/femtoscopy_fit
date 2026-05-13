#include "TFile.h"
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include <iostream>
#include <fstream>
#include <cmath>

// ====================================================================
// CMS + Heavy Ion Label
// ====================================================================

void DrawCMSHeavyIonLabel()
{
    TLatex latex;
    latex.SetNDC();

    latex.SetTextFont(62);
    latex.SetTextSize(0.05);
    latex.DrawLatex(0.155, 0.955, "CMS");

    latex.SetTextFont(52);
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.275, 0.955, "Internal");

    latex.SetTextFont(42);
    latex.SetTextSize(0.04);
    latex.DrawLatex(0.5, 0.963, "2023 PbPb, #sqrt{s_{NN}} = 5.36 TeV");
}

static const double kHbarC = 0.1973; // GeV·fm

// ------------------ Background ------------------
double background(double* x, double* p){
  double q  = x[0];

  double N  = p[0];
  double a1 = p[1];
  double R1 = p[2];
  double a2 = p[3];
  double R2 = p[4];

  return N * (1.0 + a1 * exp(-q*q*R1*R1)) * (1.0 - a2 * exp(-q*q*R2*R2));
}

double expSource(double* x, double* p)
{
    double q = x[0];

    // p[0] = N (normalisation)
    // p[1] = R (source radius, fm)
    // p[2] = lambda (correlation strength)
    // p[3] = d (baseline)
    double N      = p[0];
    double R      = p[1];
    double lambda = p[2];
    double d      = p[3];

    double BKG = 1.0 + d*q;
    
    double QS  = 1.0 + lambda * exp(-(R * q) / kHbarC);

    return N * BKG * QS;
}

// ---------------------- Multiplicity & kT bins ----------------------

const int kNch = 6;
double mBin_low[kNch] = {0.0,  5.0, 10.0, 15.0, 25.0, 35.0};
double mBin_hig[kNch] = {100.0,10.0, 15.0, 25.0, 35.0,100.0};

const int kNkT = 6;
double kT_low[kNkT]  = {0.0, 0.1, 0.3, 0.5, 0.7, 0.9};
double kT_high[kNkT] = {1.5, 0.3, 0.5, 0.7, 0.9, 1.5};
const char* kT_labels[kNkT] = {
    "all kT",
    "0.1 < k_{T} < 0.3 GeV",
    "0.3 < k_{T} < 0.5 GeV",
    "0.5 < k_{T} < 0.7 GeV",
    "0.7 < k_{T} < 0.9 GeV",
    "0.9 < k_{T} < 1.5 GeV"
};

void fit_DR_DG()
{
    gStyle->SetOptStat(0);

    double qmin     = 0.02;
    double qmax     = 0.5;

    double qmin_bkg = 0.8;
    double qmax_bkg = 2.5;

    gSystem->Exec(Form(
        "mkdir -p results_DR_doublegaus/DR_sig%.2fto%.2f_bkg%.2fto%.2f",
        qmin, qmax, qmin_bkg, qmax_bkg));

    // CSV output
    std::ofstream fout(Form(
        "results_DR_doublegaus/DR_sig%.2fto%.2f_bkg%.2fto%.2f/fit_results_threefits.csv",
        qmin, qmax, qmin_bkg, qmax_bkg));

// Header needs:
fout << "iMult,iKt,Ntrk_low,Ntrk_high,kT_low,kT_high,"
     << "Exp_N,Exp_N_err,Exp_R,Exp_R_err,Exp_lambda,Exp_lambda_err,"
     << "Exp_d,Exp_d_err,"
     << "Exp_chi2,Exp_ndf\n";

    TFile* fhist = TFile::Open("/Users/isabelamaiettosilverio/Documents/root/Femtoscopy_UPC_analysis/files/Cq_pT0p2.root", "READ");
    if (!fhist || fhist->IsZombie()) {
        std::cerr << "Cannot open Cq_pT0p2.root\n";
        return;
    }

    // Multiplicity + kT loops
    for (int iM = 0; iM < kNch; iM++) {
        double nlow = mBin_low[iM];
        double nhig = mBin_hig[iM];

        for (int iKt = 0; iKt < kNkT; iKt++) {
            double ktlow = kT_low[iKt];
            double kthig = kT_high[iKt];

            TString path = Form("Cq_avg_m%d_kT%d", iM, iKt);

            TH1F* Cq_avg = (TH1F*)fhist->Get(path);
            if (!Cq_avg) {
                std::cerr << "Cannot find " << path << "\n";
                continue;
            }
            Cq_avg->SetDirectory(nullptr); // detach from file


TF1 fBG("fBG", background, qmin_bkg,  qmax_bkg, 5);
fBG.SetParameters(1.0, 0.05, 1.0, 0.05, 1.0);
fBG.SetParLimits(0, 0.5, 1.5);
fBG.SetParLimits(1, 0.0, 1.0);  // a1
fBG.SetParLimits(2, 0.0, 5.0);  // R1
fBG.SetParLimits(3, 0.0, 1.0);  // a2
fBG.SetParLimits(4, 0.0, 5.0);  // R2

TFitResultPtr rBG = Cq_avg->Fit(&fBG, "RSQ");

            // Read chi2 and ndf from the stored fit result (consistent source)
            double chi2_bkg = rBG->Chi2();
            double ndf_bkg  = rBG->Ndf();

fBG.SetRange(0.0, qmax_bkg);  // extend down to 0 after fitting

            TH1F* DR = (TH1F*)Cq_avg->Clone(Form("DR_%d_%d_%d", iM,(int)nlow, iKt));
            for(int ib = 1; ib <= DR->GetNbinsX(); ib++){
                double q      = DR->GetBinCenter(ib);
                double bgval  = fBG.Eval(q);
                double datval = DR->GetBinContent(ib);
                double daterr = DR->GetBinError(ib);
                if(bgval > 0){
                DR->SetBinContent(ib, datval / bgval);
                DR->SetBinError  (ib, daterr / bgval); 
                DR->SetBinContent(ib, 0.0);
                DR->SetBinError  (ib, 0.0);
                }
            }

            // ---- Fit fExp to DR ----
            TF1 fExp("fExp", expSource, qmin, qmax, 4);

            fExp.SetParameters(1.0, 2.0, 0.8, 0.0);
            fExp.SetParLimits(0, 0.5, 1.5);
            fExp.SetParLimits(1, 0.5, 8.0);
            fExp.SetParLimits(2, 0.25, 2.5);
            fExp.SetParLimits(3, -0.05, 0.05);
            fExp.FixParameter(0, 1.0);
            fExp.SetRange(qmin, 0.5);
            DR->Fit(&fExp, "R0Q");

            fExp.ReleaseParameter(0);
            fExp.SetRange(qmin, qmax);
            TFitResultPtr rExp = DR->Fit(&fExp, "RS");

            double chi2 = rExp->Chi2();
            double ndf  = rExp->Ndf();

auto MakePullHist = [](TH1F* hData, TF1* fit, const char* name) -> TH1F* {
    TH1F* hPull = (TH1F*)hData->Clone(name);
    hPull->Reset();
    hPull->SetTitle("");
    for (int i = 1; i <= hData->GetNbinsX(); i++) {
        double data  = hData->GetBinContent(i);
        double err   = hData->GetBinError(i);
        double fval  = fit->Eval(hData->GetBinCenter(i));
        if (err > 0)
            hPull->SetBinContent(i, (data - fval) / err);
        else
            hPull->SetBinContent(i, 0.0);
        hPull->SetBinError(i, 1.0); // pull error is always 1 by definition
    }
    return hPull;
};

TCanvas* c = new TCanvas(Form("c_%d_%d", iM, iKt), "", 1200, 700);

// ---- Left side: two pads stacked ----
TPad* pad1_main  = new TPad("pad1_main",  "", 0.00, 0.25, 0.50, 1.00);
TPad* pad1_ratio = new TPad("pad1_ratio", "", 0.00, 0.00, 0.50, 0.25);

// ---- Right side: two pads stacked ----
TPad* pad2_main  = new TPad("pad2_main",  "", 0.50, 0.25, 1.00, 1.00);
TPad* pad2_ratio = new TPad("pad2_ratio", "", 0.50, 0.00, 1.00, 0.25);

for (TPad* p : {pad1_main, pad2_main}) {
    p->SetBottomMargin(0.02);   // no gap between main and ratio
    p->SetTopMargin(0.06);
    p->SetLeftMargin(0.15);
    p->SetRightMargin(0.05);
    p->Draw();
}
for (TPad* p : {pad1_ratio, pad2_ratio}) {
    p->SetTopMargin(0.02);
    p->SetBottomMargin(0.35);   // room for x-axis labels
    p->SetLeftMargin(0.15);
    p->SetRightMargin(0.05);
    p->SetGridy();
    p->Draw();
}

// ============================================================
// PAD 1 MAIN: C(q) data + BKG fit  (your existing code)
// ============================================================
pad1_main->cd();

            // --- Canvas & drawing ---
            gPad->SetMargin(0.15, 0.05, 0.15, 0.05);

            Cq_avg->SetTitle("");
            Cq_avg->SetMarkerStyle(24);
            Cq_avg->SetMarkerSize(1.2);
            Cq_avg->SetLineColor(kBlack);
            Cq_avg->SetMarkerColor(kBlack);
            Cq_avg->GetXaxis()->SetRangeUser(0.0, 1.5);
            Cq_avg->GetYaxis()->SetRangeUser(0.7501, 2.501);
            Cq_avg->GetXaxis()->SetTitle("q_{inv} [GeV]");
            Cq_avg->GetYaxis()->SetTitle("C(q)");
            Cq_avg->GetXaxis()->CenterTitle();
            Cq_avg->GetYaxis()->CenterTitle();
            Cq_avg->GetXaxis()->SetTitleSize(0.05);
            Cq_avg->GetXaxis()->SetLabelSize(0.045);
            Cq_avg->GetYaxis()->SetTitleSize(0.05);
            Cq_avg->GetYaxis()->SetLabelSize(0.045);

            Cq_avg->Draw("E1");

            fBG.SetLineColor(kRed);
            fBG.SetLineWidth(1);
            fBG.Draw("SAME");

            TLatex t;
            t.SetNDC();
            t.SetTextFont(42);
            t.SetTextSize(0.04);
            t.DrawLatex(0.65, 0.88, Form("%d #leq N_{trk} < %d", (int)nlow, (int)nhig));
            t.DrawLatex(0.65, 0.82, kT_labels[iKt]);
            t.SetTextSize(0.035);
            t.DrawLatex(0.565, 0.5, Form("#chi^{2}/ndf = %.2f", chi2_bkg / ndf_bkg));
            t.DrawLatex(0.21, 0.88, "DG");
            t.DrawLatex(0.21, 0.82, "Fit range:");
            t.DrawLatex(0.21, 0.78, "0.5 #leq q #leq 2.5 GeV");

            TLine lRef(0.0, 1.0, 1.5, 1.0);
            lRef.SetLineColor(kGray + 2);
            lRef.SetLineStyle(7);
            lRef.Draw("SAME");

            TLegend leg(0.5, 0.55, 0.7, 0.7);
            leg.SetBorderSize(0);
            leg.SetTextSize(0.035);
            leg.AddEntry(Cq_avg, "C(q)", "lep");
            leg.AddEntry(&fBG, "Double gauss fit", "l");
            leg.Draw();

            DrawCMSHeavyIonLabel(); 

// ============================================================
// PAD 1 RATIO: (Cq_avg - fBG_ext) / stat error
// ============================================================
pad1_ratio->cd();

TH1F* hPull1 = MakePullHist(Cq_avg, &fBG, Form("pull1_%d_%d", iM, iKt));

hPull1->GetXaxis()->SetRangeUser(0.0, 1.5);
hPull1->GetYaxis()->SetRangeUser(-3.5, 3.5);
hPull1->GetXaxis()->SetTitle("q_{inv} [GeV]");
hPull1->GetYaxis()->SetTitle("Pull");
hPull1->GetXaxis()->CenterTitle();
hPull1->GetYaxis()->CenterTitle();
// ---- large sizes because this pad is short ----
hPull1->GetXaxis()->SetTitleSize(0.14);
hPull1->GetXaxis()->SetLabelSize(0.12);
hPull1->GetYaxis()->SetTitleSize(0.13);
hPull1->GetYaxis()->SetLabelSize(0.11);
hPull1->GetYaxis()->SetTitleOffset(0.4);
hPull1->GetYaxis()->SetNdivisions(504);
hPull1->SetMarkerStyle(20);
hPull1->SetMarkerSize(0.7);
hPull1->SetLineColor(kGray+1);
hPull1->SetMarkerColor(kGray+1);
hPull1->Draw("E1");

TLine zline1(0.0, 0.0, 1.5, 0.0);
zline1.SetLineColor(kBlack); zline1.SetLineStyle(7); zline1.SetLineWidth(1);
zline1.Draw("SAME");

// ±1σ and ±2σ bands for reference
TLine sig1p(0.0,  1.0, 1.5,  1.0); sig1p.SetLineColor(kGreen+2); sig1p.SetLineStyle(3); sig1p.Draw("SAME");
TLine sig1m(0.0, -1.0, 1.5, -1.0); sig1m.SetLineColor(kGreen+2); sig1m.SetLineStyle(3); sig1m.Draw("SAME");
TLine sig2p(0.0,  2.0, 1.5,  2.0); sig2p.SetLineColor(kOrange-3); sig2p.SetLineStyle(3); sig2p.Draw("SAME");
TLine sig2m(0.0, -2.0, 1.5, -2.0); sig2m.SetLineColor(kOrange-3); sig2m.SetLineStyle(3); sig2m.Draw("SAME");

// ============================================================
// PAD 2 MAIN: Double Ratio + fits  (your existing code)
// ============================================================
pad2_main->cd();

DR->GetListOfFunctions()->Clear();
DR->SetTitle("");
DR->SetMarkerStyle(25);
DR->SetMarkerSize(1.1);
DR->SetLineWidth(1);
DR->SetLineColor(kBlack);
DR->SetMarkerColor(kBlack);
DR->GetXaxis()->SetRangeUser(0.0, 1.5);
DR->GetYaxis()->SetRangeUser(0.7501, 2.501);
// ---- suppress x-axis labels on main pad ----
DR->GetXaxis()->SetTitle("");
DR->GetXaxis()->SetLabelSize(0);
DR->GetYaxis()->SetTitle("Double Ratio DR(q)");
DR->GetYaxis()->SetTitleSize(0.06);
DR->GetYaxis()->SetLabelSize(0.055);
DR->Draw("E1");

            fExp.SetLineColor(kBlue);
            fExp.SetLineWidth(1);
            fExp.Draw("SAME");

t.SetTextSize(0.04);
t.DrawLatex(0.58, 0.88, Form("%d #leq N_{trk} < %d", (int)nlow, (int)nhig));
t.DrawLatex(0.58, 0.82, kT_labels[iKt]);
t.DrawLatex(0.21, 0.88, "DR = SR/SG");
t.DrawLatex(0.21, 0.82, "Fit range SIG:");
t.DrawLatex(0.21, 0.78, Form("%.2f #leq q < %.1f GeV", qmin, qmax));

TLine line2(0.0, 1.0, 1.5, 1.0);
line2.SetLineColor(kBlack); line2.SetLineStyle(7); line2.SetLineWidth(1);
line2.Draw("SAME");

TLegend leg_DR(0.52, 0.58, 0.82, 0.80);
leg_DR.SetBorderSize(0); leg_DR.SetTextSize(0.04);
leg_DR.AddEntry(DR,      "C(q) data",                           "lep");
// Correct — shows chi2/ndf
leg_DR.AddEntry(&fExp, Form("Exp,  #chi^{2}/ndf = %.2f", chi2/ndf), "l");
//leg_DR.AddEntry(&fLevy,  Form("Levy, #chi^{2}/ndf = %.2g", chi2_levy), "l");
//leg_DR.AddEntry(&fExp,   Form("Exp,  #chi^{2}/ndf = %.2g", chi2_exp),  "l");
//leg_DR.AddEntry(&fGaus,  Form("Gauss, #chi^{2}/ndf = %.2g", chi2_gaus),"l");
leg_DR.Draw();
DrawCMSHeavyIonLabel();

// ============================================================
// PAD 2 RATIO: three pulls, one per fit model
// ============================================================
pad2_ratio->cd();


//TH1F* hPull_Levy = MakePullHist(DR, &fLevy, Form("pull_levy_%d_%d", iM, iKt));
TH1F* hPull_Exp  = MakePullHist(DR, &fExp,  Form("pull_exp_%d_%d",  iM, iKt));
//TH1F* hPull_Gaus = MakePullHist(DR, &fGaus, Form("pull_gaus_%d_%d", iM, iKt));

// Style all three
for (auto* h : {hPull_Exp}) {
    h->GetXaxis()->SetRangeUser(0.0, 1.5);
    h->GetYaxis()->SetRangeUser(-3.5, 3.5);
    h->GetXaxis()->SetTitle("q_{inv} [GeV]");
    h->GetYaxis()->SetTitle("Pull");
    h->GetXaxis()->CenterTitle(); h->GetYaxis()->CenterTitle();
    h->GetXaxis()->SetTitleSize(0.14); h->GetXaxis()->SetLabelSize(0.12);
    h->GetYaxis()->SetTitleSize(0.13); h->GetYaxis()->SetLabelSize(0.11);
    h->GetYaxis()->SetTitleOffset(0.4);
    h->GetYaxis()->SetNdivisions(504);
    h->SetMarkerSize(0.5);
}

//hPull_Levy->SetMarkerStyle(20); // filled circle like DR
hPull_Exp ->SetMarkerStyle(20);
///hPull_Gaus->SetMarkerStyle(20);

//->SetLineColor(kRed);     hPull_Levy->SetMarkerColor(kRed);
hPull_Exp ->SetLineColor(kBlue);    hPull_Exp ->SetMarkerColor(kBlue);
//Pull_Gaus->SetLineColor(kGreen+2); hPull_Gaus->SetMarkerColor(kGreen+2);

//hPull_Levy->Draw("E1");
hPull_Exp ->Draw("E1");
//hPull_Gaus->Draw("E1 SAME");

TLine zline2(0.0, 0.0, 1.5, 0.0);
zline2.SetLineColor(kBlack); zline2.SetLineStyle(7); zline2.SetLineWidth(1);
zline2.Draw("SAME");

TLine s1p2(0.0,  1.0, 1.5,  1.0); s1p2.SetLineColor(kGreen+2); s1p2.SetLineStyle(3); s1p2.Draw("SAME");
TLine s1m2(0.0, -1.0, 1.5, -1.0); s1m2.SetLineColor(kGreen+2); s1m2.SetLineStyle(3); s1m2.Draw("SAME");
TLine s2p2(0.0,  2.0, 1.5,  2.0); s2p2.SetLineColor(kOrange-3); s2p2.SetLineWidth(1); s2p2.SetLineStyle(3); s2p2.Draw("SAME");
TLine s2m2(0.0, -2.0, 1.5, -2.0); s2m2.SetLineColor(kOrange-3); s2m2.SetLineStyle(3); s2m2.Draw("SAME");


            c->Update();
            c->SaveAs(Form(
                "results_DR_doublegaus/DR_sig%.2fto%.2f_bkg%.2fto%.2f/"
                "DR_bkg%.2fto%.2f_sig%.2fto%.2f_m%d_kT%d.pdf",
                qmin, qmax, qmin_bkg, qmax_bkg,
                qmin, qmax, qmin_bkg, qmax_bkg,
                iM, iKt));

            delete c;

            // Extract fit parameters
            double N      = fExp.GetParameter(0);
            double Nerr   = fExp.GetParError(0);
            double R      = fExp.GetParameter(1);
            double Rerr   = fExp.GetParError(1);
            double d    = fExp.GetParameter(3);
            double derr = fExp.GetParError(3);
            double lambda = fExp.GetParameter(2);
            double lamerr = fExp.GetParError(2);

            // p[0]=N, p[1]=d (linear), p[2]=e (quadratic), p[3]=R, p[4]=lambda
            std::cout << Form(
                "  Exp: chi2/ndf=%.2f  N=%.3f  R=%.3f  lambda=%.3f\n",
                chi2 / ndf, N, R, lambda);

            fout << iM    << "," << iKt   << ","
                 << nlow  << "," << nhig  << ","
                 << ktlow << "," << kthig << ","
                 << N      << "," << Nerr   << ","
                 << R      << "," << Rerr   << ","
                 << lambda << "," << lamerr << ","
                 << d << "," << derr << ","
                 << chi2   << "," << ndf
                 << "\n";

        } // kT loop
    } // multiplicity loop

    fout.close();

    fhist->Close();

    std::cout << "Fit finished. CSV saved.\n";
}