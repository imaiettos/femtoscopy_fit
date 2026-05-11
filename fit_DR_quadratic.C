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

double QUAD_expSource(double* x, double* p)
{
    double q = x[0];

    // p[0] = N (normalisation)
    // p[1] = d (linear background coefficient)
    // p[2] = e (quadratic background coefficient)
    // p[3] = R (source radius, fm)
    // p[4] = lambda (correlation strength)
    double N      = p[0];
    double d      = p[1];
    double e      = p[2];
    double R      = p[3];
    double lambda = p[4];

    double BKG = 1.0 + d*q + e*q*q;
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

void fit_DR_quadratic()
{
    gStyle->SetOptStat(0);

    double qmin     = 0.02;
    double qmax     = 2.5;

    double qmin_bkg = 1.0;
    double qmax_bkg = 2.5;

    gSystem->Exec(Form(
        "mkdir -p results_DR_quadratic/DR_sig%.2fto%.2f_bkg%.2fto%.2f",
        qmin, qmax, qmin_bkg, qmax_bkg));

    // CSV output
    std::ofstream fout(Form(
        "results_DR_quadratic/DR_sig%.2fto%.2f_bkg%.2fto%.2f/fit_results_threefits.csv",
        qmin, qmax, qmin_bkg, qmax_bkg));

    fout << "iMult,iKt,Ntrk_low,Ntrk_high,kT_low,kT_high,"
         << "Exp_N,Exp_N_err,Exp_R,Exp_R_err,Exp_lambda,Exp_lambda_err,"
         << "Exp_chi2,Exp_ndf\n";

    TFile* fhist = TFile::Open("Cq_pT0p2.root", "READ");
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

            // SetRange is used to restrict it during each fit step.
            TF1 fQDDR("fQDDR", QUAD_expSource, qmin, qmax, 5);

            // Initial parameters and limits
            fQDDR.SetParameters(1.0, 0.0, 0.0, 2.0, 0.8);      
            fQDDR.SetParLimits(0,  0.5,  1.5);  // N
            fQDDR.SetParLimits(1, -0.05, 0.05);  // d
            fQDDR.SetParLimits(2, -0.01, 0.01);  // e
            fQDDR.SetParLimits(3,  0.1,  8.0);   // R
            fQDDR.SetParLimits(4,  0.0,  2.5);   // lambda

            // Step 1: Background-only pre-fit in [qmin_bkg, qmax_bkg]
            fQDDR.FixParameter(4, 1.0);
            fQDDR.SetRange(qmin_bkg, qmax_bkg);
            Cq_avg->Fit(&fQDDR, "R0Q"); // R=range, 0=don't draw, Q=quiet

            // Step 2: Full signal fit — release lambda, use full range
            fQDDR.ReleaseParameter(4);
            fQDDR.SetRange(qmin, qmax);

            // and used "R" (no "S"), which left rFQDR unpopulated.
            TFitResultPtr rQDDR = Cq_avg->Fit(&fQDDR, "RS");

            // Read chi2 and ndf from the stored fit result (consistent source)
            double chi2 = rQDDR->Chi2();
            double ndf  = rQDDR->Ndf();

            // --- Canvas & drawing ---
            TCanvas* c = new TCanvas(Form("c_%d_%d", iM, iKt), "", 700, 700);
            gPad->SetMargin(0.15, 0.05, 0.15, 0.05);

            Cq_avg->SetTitle("");
            Cq_avg->SetMarkerStyle(24);
            Cq_avg->SetMarkerSize(1.2);
            Cq_avg->SetLineColor(kBlack);
            Cq_avg->SetMarkerColor(kBlack);
            Cq_avg->GetXaxis()->SetRangeUser(0.0, 1.5);
            Cq_avg->GetYaxis()->SetRangeUser(0.701, 3.501);
            Cq_avg->GetXaxis()->SetTitle("q_{inv} [GeV]");
            Cq_avg->GetYaxis()->SetTitle("C(q)");
            Cq_avg->GetXaxis()->CenterTitle();
            Cq_avg->GetYaxis()->CenterTitle();
            Cq_avg->GetXaxis()->SetTitleSize(0.05);
            Cq_avg->GetXaxis()->SetLabelSize(0.045);
            Cq_avg->GetYaxis()->SetTitleSize(0.05);
            Cq_avg->GetYaxis()->SetLabelSize(0.045);

            Cq_avg->Draw("E1");

            fQDDR.SetLineColor(kRed);
            fQDDR.SetLineWidth(1);
            fQDDR.Draw("SAME");

            TLatex t;
            t.SetNDC();
            t.SetTextFont(42);
            t.SetTextSize(0.04);
            t.DrawLatex(0.65, 0.88, Form("%d #leq N_{trk} < %d", (int)nlow, (int)nhig));
            t.DrawLatex(0.65, 0.82, kT_labels[iKt]);
            t.SetTextSize(0.035);
            t.DrawLatex(0.565, 0.5, Form("#chi^{2}/ndf = %.2f", chi2 / ndf));
            t.DrawLatex(0.21, 0.88, "SR = N(1+#lambda e^{-(qR)})(1+dq+eq^{2})");
            t.DrawLatex(0.21, 0.82, "Fit range:");
            t.DrawLatex(0.21, 0.78, "0.02 #leq q #leq 2.5 GeV");

            TLine lRef(0.0, 1.0, 1.5, 1.0);
            lRef.SetLineColor(kGray + 2);
            lRef.SetLineStyle(7);
            lRef.Draw("SAME");

            TLegend leg(0.55, 0.55, 0.75, 0.7);
            leg.SetBorderSize(0);
            leg.SetTextSize(0.035);
            leg.AddEntry(Cq_avg, "C(q)", "lep");
            leg.AddEntry(&fQDDR, "SR with quad bkg", "l");
            leg.Draw();

            DrawCMSHeavyIonLabel(); //

            c->Update();
            c->SaveAs(Form(
                "results_DR_quadratic/DR_sig%.2fto%.2f_bkg%.2fto%.2f/"
                "DR_bkg%.2fto%.2f_sig%.2fto%.2f_m%d_kT%d.pdf",
                qmin, qmax, qmin_bkg, qmax_bkg,
                qmin, qmax, qmin_bkg, qmax_bkg,
                iM, iKt));

            delete c;

            // Extract fit parameters
            double N      = fQDDR.GetParameter(0);
            double Nerr   = fQDDR.GetParError(0);
            double d_par  = fQDDR.GetParameter(1);
            double e_par  = fQDDR.GetParameter(2);
            double R      = fQDDR.GetParameter(3);
            double Rerr   = fQDDR.GetParError(3);
            double lambda = fQDDR.GetParameter(4);
            double lamerr = fQDDR.GetParError(4);

            // p[0]=N, p[1]=d (linear), p[2]=e (quadratic), p[3]=R, p[4]=lambda
            std::cout << Form(
                "  SGDR: chi2/ndf=%.2f  N=%.3f  d=%.4f  e=%.4f  R=%.3f  lambda=%.3f\n",
                chi2 / ndf, N, d_par, e_par, R, lambda);

            fout << iM    << "," << iKt   << ","
                 << nlow  << "," << nhig  << ","
                 << ktlow << "," << kthig << ","
                 << N      << "," << Nerr   << ","
                 << R      << "," << Rerr   << ","
                 << lambda << "," << lamerr << ","
                 << chi2   << "," << ndf
                 << "\n";

        } // kT loop
    } // multiplicity loop

    fout.close();

    fhist->Close();

    std::cout << "Fit finished. CSV saved.\n";
}
