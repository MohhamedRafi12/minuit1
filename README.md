# minuit1

- expFit.cpp: example using the more general fitting interface in minuit
- expFit.ipynb: equivalent example using lmfit
- SimultaneousExps(lm).ipynb: generation of histograms with correlated signals for simultaneous fit exercise
- rootExample.cpp: just another example of using ROOT classes in a C++ program
- *.root files: various input histograms for fitting exercises

-----

Team member names and computing IDs: Mohhamed Rafi (fxc7ss)

-----

Exercise 1 comments:
--

The gambit fit for both Chi2 and NLL produces better fits compared to the Gaussian dist.

-----



-----

Exercise 2 comments:
--

A     = 262.707 ± 29.4797
mu    = 74.698 ± 0.473986
sigma = 4.39042 ± 0.430573

Exp1 Background: B=790.208 lambda=20.1009
Exp2 Background: C=5000 n=-0.920318

Chi2_1/NDF = 32.8261/45 = 0.729468  Prob = 0.911411
Chi2_2/NDF = 49.6601/45 = 1.10356  Prob = 0.292883

The the reduced chisq is approx. one so I think we got a good fit. I included a scaling offset to get a better fit.


-----

Exercise 3 comments:
--

p[0] A = 53.997760 +/- 0.521500
p[1] mu_x = 3.516810 +/- 0.005432
p[2] sigma_x = 0.989221 +/- 0.007118
p[3] mu_y = 1.903694 +/- 0.015122
p[4] sigma_y = 1.953785 +/- 0.017823
p[5] bkg_scale = 0.245090 +/- 0.002289

=== Estimated Signal Yield ===
N_signal = 327.865 +/- 4.954
Where N_signal = A pi sigma_1 sigma_2
-----
