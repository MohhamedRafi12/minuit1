// ex3_minuit.C
// Usage: root -l -q 'ex3_minuit.C("fitInputs.root","ex3.pdf")'

#include "TFile.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TMinuit.h"
#include "TMath.h"
#include "TString.h"
#include <iostream>

static TH2* gHData = nullptr;
static TH2* gHBkg  = nullptr;

Double_t model(Double_t x, Double_t y, const Double_t* p) {
  double A  = p[0];
  double mx = p[1];
  double sx = p[2];
  double my = p[3];
  double sy = p[4];
  double k  = p[5];
  double sig = A * TMath::Exp(-((x - mx)*(x - mx))/(sx*sx))
                 * TMath::Exp(-((y - my)*(y - my))/(sy*sy));
  double bkg = 0.0;
  if (gHBkg) bkg = k * gHBkg->Interpolate(x, y);
  return sig + bkg;
}

void fcn(Int_t& /*npar*/, Double_t* /*gin*/, Double_t& f, Double_t* par, Int_t /*iflag*/) {
  double chi2 = 0.0;
  int nx = gHData->GetNbinsX();
  int ny = gHData->GetNbinsY();
  for (int ix = 1; ix <= nx; ++ix) {
    double x = gHData->GetXaxis()->GetBinCenter(ix);
    for (int iy = 1; iy <= ny; ++iy) {
      double y = gHData->GetYaxis()->GetBinCenter(iy);
      double obs = gHData->GetBinContent(ix, iy);
      if (obs <= 0) continue;
      double expv = model(x, y, par);
      double err  = TMath::Sqrt(obs);
      chi2 += ((obs - expv)*(obs - expv)) / (err*err);
    }
  }
  f = chi2;
}

void ex3(const char* infile="fitInputs.root", const char* outfile="ex3.pdf") {
  gStyle->SetOptStat(0);

  TFile* f = TFile::Open(infile);
  if (!f || f->IsZombie()) { std::cerr << "Cannot open " << infile << "\n"; return; }
  gHData = (TH2*) f->Get("hdata");
  gHBkg  = (TH2*) f->Get("hbkg");
  if (!gHData || !gHBkg) { std::cerr << "Missing histograms!\n"; return; }

  const int npar = 6;
  TMinuit minuit(npar);
  minuit.SetFCN(fcn);

  // name, start, step, min, max
  minuit.DefineParameter(0, "A",         100.0,  1.0,   0.0, 1e5);
  minuit.DefineParameter(1, "mu_x",        2.5,  0.1,   0.0, 6.0);
  minuit.DefineParameter(2, "sigma_x",     0.8,  0.05,  0.1, 3.0);
  minuit.DefineParameter(3, "mu_y",        2.5,  0.1,   0.0, 6.0);
  minuit.DefineParameter(4, "sigma_y",     0.8,  0.05,  0.1, 3.0);
  minuit.DefineParameter(5, "bkg_scale",   1.0,  0.05,  0.0, 10.0);

  Double_t arglist[2]; Int_t ierflg = 0;
  arglist[0] = 5000; arglist[1] = 1.0;
  minuit.mnexcm("MIGRAD", arglist, 2, ierflg);

  Double_t par[npar], err[npar];
  std::cout << "\n=== Best-fit parameters ===\n";
  for (int i = 0; i < npar; ++i) {
    TString pname; Double_t v, e, b1, b2; Int_t ivar;
    minuit.mnpout(i, pname, v, e, b1, b2, ivar);
    par[i] = v; err[i] = e;
    std::cout << Form("p[%d] %s = %.6f +/- %.6f\n", i, pname.Data(), par[i], err[i]);
  }

  TH2D* hfit = (TH2D*) gHData->Clone("hfit");
  TH2D* hres = (TH2D*) gHData->Clone("hres");
  TH2D* hsub = (TH2D*) gHData->Clone("hsub");

  int nx = gHData->GetNbinsX();
  int ny = gHData->GetNbinsY();
  for (int ix = 1; ix <= nx; ++ix) {
    double x = gHData->GetXaxis()->GetBinCenter(ix);
    for (int iy = 1; iy <= ny; ++iy) {
      double y = gHData->GetYaxis()->GetBinCenter(iy);
      double obs  = gHData->GetBinContent(ix, iy);
      double fval = model(x, y, par);
      double bval = par[5] * gHBkg->Interpolate(x, y);
      hfit->SetBinContent(ix, iy, fval);
      hres->SetBinContent(ix, iy, obs - fval);
      hsub->SetBinContent(ix, iy, obs - bval);
    }
  }

  TCanvas* c = new TCanvas("c", "Templated Background Subtraction (MINUIT)", 1200, 900);
  c->Divide(2,2);
  c->cd(1); gHData->SetTitle("Data;X;Y"); gHData->Draw("lego2");
  c->cd(2); hfit ->SetTitle("Fit (Signal + Background);X;Y"); hfit->Draw("lego2");
  c->cd(3); hres ->SetTitle("Residuals (Data - Fit);X;Y");   hres->Draw("lego2");
  c->cd(4); hsub ->SetTitle("Data - Background;X;Y");        hsub->Draw("lego2");
  c->SaveAs(outfile);

  double As = par[0], sx = par[2], sy = par[4];
  double yield = As * TMath::Pi() * sx * sy;
  double relErr = TMath::Sqrt(TMath::Power(err[0]/As,2) + TMath::Power(err[2]/sx,2) + TMath::Power(err[4]/sy,2));
  double errYield = yield * relErr;

  std::cout << "\n=== Estimated Signal Yield ===\n";
  std::cout << Form("N_signal = %.3f +/- %.3f\n", yield, errYield) << std::endl;
}
