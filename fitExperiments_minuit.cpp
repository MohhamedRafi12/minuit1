#include <iostream>
#include <cmath>
#include <TROOT.h>
#include <TFile.h>
#include <TH1F.h>
#include <TF1.h>      
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TMinuit.h>
#include <TMath.h>

using namespace std;

// =======================================================
// Globals for MINUIT shared fit
TH1F *h1 = nullptr;
TH1F *h2 = nullptr;

// =======================================================
// Model components
double gauss(double xx, double A, double mu, double sigma) {
    return A * TMath::Gaus(xx, mu, sigma, true);
}
double bkg_exp(double xx, double B, double lambda) {
    return B * TMath::Exp(-xx / lambda);
}
double bkg_power(double xx, double C, double n) {
    if (xx <= 0) return 0;
    return C * TMath::Power(xx, n);
}

// =======================================================
// Composite models
double model_exp(double *x, double *p) {
    double xx = x[0];
    return gauss(xx, p[0], p[1], p[2]) + bkg_exp(xx, p[3], p[4]);
}
double model_pow(double *x, double *p) {
    double xx = x[0];
    double gaus = p[0]*TMath::Gaus(xx,p[1],p[2],true);
    double bkg  = p[5]*TMath::Power(xx,p[6]) + p[7];
    return gaus + bkg;
}

// =======================================================
double calcChi2(TH1F *h, double (*model)(double*, double*), double *p) {
    double chi2 = 0.0;
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        double obs = h->GetBinContent(i);
        if (obs <= 0) continue;
        double xx = h->GetBinCenter(i);      
        double expv = model(&xx, p);
        if (expv <= 0) expv = 1e-9;
        chi2 += pow((obs - expv), 2) / obs;
    }
    return chi2;
}

// =======================================================
void fcn(int &npar, double *grad, double &f, double par[], int flag) {
    double chi2_1 = calcChi2(h1, model_exp, par);
    double chi2_2 = calcChi2(h2, model_pow, par);
    f = chi2_1 + chi2_2;
}

// =======================================================
int main() {
    gStyle->SetOptStat(0);
    gStyle->SetTitleFontSize(0.05);

    // open file
    TFile *file = TFile::Open("experiments.root");
    if (!file || file->IsZombie()) {
        cerr << "Cannot open experiments.root" << endl;
        return 1;
    }

    h1 = (TH1F*)file->Get("hexp1");
    h2 = (TH1F*)file->Get("hexp2");
    if (!h1 || !h2) {
        cerr << "Missing histograms hexp1 or hexp2" << endl;
        return 1;
    }

    const int npar = 8;  // [A, mu, sigma, B, lambda, C, n]
    TMinuit minuit(npar);
    minuit.SetFCN(fcn);

    minuit.DefineParameter(0, "A", 400.0, 10.0, 0.0, 5000.0);
    minuit.DefineParameter(1, "mu", 78.0, 0.2, 70.0, 85.0);
    minuit.DefineParameter(2, "sigma", 4.0, 0.1, 2.0, 8.0);

    // exponential background
    minuit.DefineParameter(3, "B", 150.0, 5.0, 0.0, 5000.0);
    minuit.DefineParameter(4, "lambda", 30.0, 1.0, 5.0, 60.0);

    // power-law background
    minuit.DefineParameter(5, "C", 100.0, 5.0, 0.0, 5000.0);
    minuit.DefineParameter(6, "n", -1.5, 0.05, -3.0, -0.5);
    minuit.DefineParameter(7, "D", 0.0, 1.0, -100.0, 100.0);

    minuit.Migrad();

    double par[npar], err[npar];
    for (int i = 0; i < npar; ++i)
        minuit.GetParameter(i, par[i], err[i]);

    // Compute stats
    double chi2_1 = calcChi2(h1, model_exp, par);
    double chi2_2 = calcChi2(h2, model_pow, par);
    int ndof1 = h1->GetNbinsX() - 5;
    int ndof2 = h2->GetNbinsX() - 5;

    double prob1 = TMath::Prob(chi2_1, ndof1);
    double prob2 = TMath::Prob(chi2_2, ndof2);

    cout << "\n===== Shared Gaussian Signal =====" << endl;
    cout << "A     = " << par[0] << " ± " << err[0] << endl;
    cout << "mu    = " << par[1] << " ± " << err[1] << endl;
    cout << "sigma = " << par[2] << " ± " << err[2] << endl;

    cout << "\nExp1 Background: B=" << par[3] << " lambda=" << par[4]
         << "\nExp2 Background: C=" << par[5] << " n=" << par[6] << endl;

    cout << "\nChi2_1/NDF = " << chi2_1 << "/" << ndof1 << " = " << chi2_1 / ndof1
         << "  Prob = " << prob1 << endl;
    cout << "Chi2_2/NDF = " << chi2_2 << "/" << ndof2 << " = " << chi2_2 / ndof2
         << "  Prob = " << prob2 << endl;

    // Draw fits
    TCanvas *c = new TCanvas("c", "Simultaneous Fits", 1200, 600);
    c->Divide(2, 1);

    // exp background
    c->cd(1);
    h1->SetMarkerStyle(20);
    h1->Draw("E");
    TF1 fexp("fexp", model_exp, h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax(), npar);
    fexp.SetParameters(par);
    fexp.SetLineColor(kRed);
    fexp.SetLineWidth(3);
    fexp.Draw("same");

    TLatex t1;
    t1.SetNDC();
    t1.SetTextSize(0.035);
    t1.DrawLatex(0.15, 0.85, "Experiment 1 (Gaussian + exp(-x/lambda))");
    t1.DrawLatex(0.15, 0.78, Form("Chi2/ndf = %.2f/%d = %.2f", chi2_1, ndof1, chi2_1 / ndof1));
    t1.DrawLatex(0.15, 0.72, Form("Prob = %.3f", prob1));

    // power-law background
    c->cd(2);
    h2->SetMarkerStyle(20);
    h2->Draw("E");
    TF1 fpow("fpow", model_pow, h2->GetXaxis()->GetXmin(), h2->GetXaxis()->GetXmax(), npar);
    fpow.SetParameters(par);
    fpow.SetLineColor(kBlue);
    fpow.SetLineWidth(3);
    fpow.Draw("same");

    TLatex t2;
    t2.SetNDC();
    t2.SetTextSize(0.035);
    t2.DrawLatex(0.15, 0.85, "Experiment 2 (Gaussian + x^{n})");
    t2.DrawLatex(0.15, 0.78, Form("Chi2/ndf = %.2f/%d = %.2f", chi2_2, ndof2, chi2_2 / ndof2));
    t2.DrawLatex(0.15, 0.72, Form("Prob = %.3f", prob2));

    c->SaveAs("ex2.pdf");
    cout << "\nSaved combined fit figure as ex2.pdf\n";
    return 0;
}
