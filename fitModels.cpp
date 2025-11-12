#include <iostream>
#include <cmath>
#include <TROOT.h>
#include <TFile.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TF1.h>
#include <TVirtualPad.h>
#include <TMinuit.h>
#include <TMath.h>

using namespace std;

// Globals for MINUIT
TH1F *hdata = nullptr;
TF1  *fparam = nullptr;
bool  useNLL = false; // switch between chi2 and NLL

// =======================================================
// Model 1: Sum of two Gaussians
double twoGauss(double *x, double *par) {
    double xx = x[0];
    double A1 = par[0];
    double mu1 = par[1];
    double sigma1 = par[2];
    double A2 = par[3];
    double mu2 = par[4];
    double sigma2 = par[5];
    return A1 * TMath::Gaus(xx, mu1, sigma1, true)
         + A2 * TMath::Gaus(xx, mu2, sigma2, true);
}

// Model 2: Gumbel distribution f(x; mu,beta) = A * exp(-(x - mu)/beta - exp(-(x - mu)/beta))
double gumbel(double *x, double *par) {
    double xx = x[0];
    double A = par[0];
    double mu = par[1];
    double beta = par[2];
    double z = (xx - mu) / beta;
    return A * exp(-z - exp(-z));
}

// =======================================================
// Chi2 objective
double calcChi2(TH1F *h, TF1 *f) {
    double chi2 = 0.0;
    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        double obs = h->GetBinContent(i);
        double expv = f->Eval(h->GetBinCenter(i));
        if (obs > 0) chi2 += pow((obs - expv) / sqrt(obs), 2);
    }
    return chi2;
}

// Poisson NLL objective
double calcNLL(TH1F* h, TF1* f) {
  double nll = 0.0;                     
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    const double n  = h->GetBinContent(i);
    double       mu = f->Eval(h->GetBinCenter(i));
    if (mu < 1e-12) mu = 1e-12;          
    nll += (mu - n * TMath::Log(mu));   
  }
  return 2.0 * nll;                     
}


double calcDeviance(TH1F* h, TF1* f) {
  double D = 0.0;
  for (int i = 1; i <= h->GetNbinsX(); ++i) {
    const double n  = h->GetBinContent(i);
    double       mu = f->Eval(h->GetBinCenter(i));
    if (mu < 1e-12) mu = 1e-12;
    if (n > 0.0) {
      D += 2.0 * (mu - n + n * TMath::Log(n / mu));
    } else {
      D += 2.0 * mu; // n=0 term
    }
  }
  return D;
}

// =======================================================
void fcn(int &npar, double *grad, double &f, double par[], int flag) {
    for (int i = 0; i < npar; ++i) fparam->SetParameter(i, par[i]);
    f = useNLL ? calcNLL(hdata, fparam) : calcChi2(hdata, fparam);
}

// =======================================================
void doFit(const char *title, TF1 *model, int npar, bool nllFit, TVirtualPad *pad) {
    pad->cd();
    useNLL = nllFit;

    // --- work on a private clone so titles/states don't bleed across pads ---
    TH1F *h = (TH1F*)hdata->Clone(Form("h_%s_%s", title, nllFit ? "nll" : "chi2"));
    h->SetDirectory(0);                         // detach from file
    h->SetTitle(Form("%s (%s)", title, nllFit ? "NLL Fit" : "Chi2 Fit"));
    h->GetXaxis()->SetTitle("x");
    h->GetYaxis()->SetTitle("Counts");
    h->SetMarkerStyle(20);
    h->SetMarkerSize(0.8);
    h->Draw("E");

    // --- set up Minuit on this histogram ---
    TMinuit minuit(npar);
    minuit.SetFCN(fcn);

    std::vector<double> par(npar), step(npar);
    for (int i = 0; i < npar; ++i) {
        par[i]  = model->GetParameter(i);
        step[i] = std::fabs(par[i]) * 0.1 + 1e-3;
        minuit.DefineParameter(i, Form("p%d", i), par[i], step[i], 0, 0);
    }

    // point the global used by fcn() at our private hist and model
    fparam = model;
    hdata  = h;

    minuit.Migrad();

    std::vector<double> outpar(npar), err(npar);
    for (int i = 0; i < npar; ++i)
        minuit.GetParameter(i, outpar[i], err[i]);
    model->SetParameters(outpar.data());

    // --- draw fitted curve ---
    model->SetLineColor(kRed);
    model->SetLineWidth(3);
    model->Draw("SAME");

    // --- recompute objective & GOF from final params (more robust than fmin) ---
    const double ndof = h->GetNbinsX() - npar;
    TLatex t; t.SetNDC(); t.SetTextSize(0.038);

    if (!nllFit) {
        const double chi2   = calcChi2(h, model);
        const double redchi = chi2 / ndof;
        const double pchi   = TMath::Prob(chi2, ndof);

        t.DrawLatex(0.12, 0.86, "Chi^{2} Fit");
        t.DrawLatex(0.12, 0.80, Form("#chi^{2}/ndf = %.2f / %.0f = %.2f", chi2, ndof, redchi));
        t.DrawLatex(0.12, 0.74, Form("Prob = %.3f", pchi));
    } else {
        const double twoNLL = calcNLL(h, model);        // always >= 0
        const double red2N  = twoNLL / ndof;            // scale-only
        const double D      = calcDeviance(h, model);   // Poisson GOF
        const double pdev   = TMath::Prob(D, ndof);

        t.DrawLatex(0.12, 0.86, "Poisson NLL Fit");
        t.DrawLatex(0.12, 0.80, Form("2NLL = %.2f  (2NLL/ndf = %.2f)", twoNLL, red2N));
        t.DrawLatex(0.12, 0.74, Form("Deviance = %.2f,  Prob = %.3f", D, pdev));
    }

    pad->Modified();
    pad->Update();
}


// =======================================================
int main() {
    gStyle->SetOptStat(0);

    TFile *file = TFile::Open("distros.root");
    if (!file || file->IsZombie()) {
        cerr << "Error: cannot open distros.root\n";
        return 1;
    }

    TH1F *h1 = (TH1F *)file->Get("dist1");
    if (!h1) {
        cerr << "Error: dist1 not found\n";
        return 1;
    }
    hdata = h1;

    TCanvas *c = new TCanvas("c", "Fits", 1600, 800);
    c->Divide(2, 2);

    // ----- Two-Gaussian: Chi2 -----
    c->cd(1);
    TF1 *gauss2_chi2 = new TF1("Sum of Two Gaussians (Chi2 Fit)", twoGauss,
                            h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax(), 6);
    gauss2_chi2->SetParameters(h1->GetMaximum(), h1->GetMean() - 1, 1,
                            h1->GetMaximum() / 2, h1->GetMean() + 1, 1);
    doFit("Sum of Two Gaussians", gauss2_chi2, 6, false, gPad);

    // ----- Two-Gaussian: NLL -----
    c->cd(2);
    TF1 *gauss2_nll = new TF1("Sum of Two Gaussians (NLL Fit)", twoGauss,
                            h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax(), 6);
    gauss2_nll->SetParameters(h1->GetMaximum(), h1->GetMean() - 1, 1,
                            h1->GetMaximum() / 2, h1->GetMean() + 1, 1);
    doFit("Sum of Two Gaussians", gauss2_nll, 6, true, gPad);

    // ----- Gumbel: Chi2 -----
    c->cd(3);
    TF1 *gumb_chi2 = new TF1("Gumbel Distribution (Chi2 Fit)", gumbel,
                            h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax(), 3);
    gumb_chi2->SetParameters(h1->GetMaximum(), h1->GetMean(), 1.0);
    doFit("Gumbel Distribution", gumb_chi2, 3, false, gPad);

    // ----- Gumbel: NLL -----
    c->cd(4);
    TF1 *gumb_nll = new TF1("Gumbel Distribution (NLL Fit)", gumbel,
                            h1->GetXaxis()->GetXmin(), h1->GetXaxis()->GetXmax(), 3);
    gumb_nll->SetParameters(h1->GetMaximum(), h1->GetMean(), 1.0);
    doFit("Gumbel Distribution", gumb_nll, 3, true, gPad);

    c->SaveAs("ex1.pdf");
    cout << "\nSaved plot: ex1.pdf\n";
    return 0;
}
