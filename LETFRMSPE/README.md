<div style="margin: 0; padding: 0; text-align: center; border: none;">
<a href="https://quantlet.com" target="_blank" style="text-decoration: none; border: none;">
<img src="https://github.com/StefanGam/test-repo/blob/main/quantlet_design.png?raw=true" alt="Header Image" width="100%" style="margin: 0; padding: 0; display: block; border: none;" />
</a>
</div>

```
Name of Quantlet: LETFRMSPE

Published in: Leveraged ETF options implied volatility paradox

Description: Computes and plots the root mean squared prediction error of the dynamic one-step-ahead forecast of the implied volatility surface of the SPY LETF call option. The estimation of the dynamic semiparametric factor model with B-spline basis is performed for this purpose

Keywords: DSFM, dynamic, semiparametric, semiparametric model, pca, principal-component-analysis, factor, factor-model,

See also: LETFFactorFuncs, LETFStochLoads, LETFIVSurfPlot

Author: Sergey Nasekin

Submitted: 2017/01/16

Datafile: SPYFULL.mat, SSOFULL.mat, yrates1415.mat

Input: 
- Km: 'B-spline order in moneyness direction
- Kt: 'B-spline order in time-to-maturity direction
- dim_mon: number of grid points for estimation in moneyness direction
- dim_ttm: number of grid points for estimation in time-to-maturity direction
- ikmon: 'parameter for setting the number of B-spline knots in moneyness direction
- ikttm: 'parameter for setting the number of B-spline knots in time-to-maturity direction
- tol: convergence tolerance for the Newton method
- maxiter: maximal number of iterations for the Newton method
- L: number of factor functions in the model
- lower_mon: lowest moneyness value to keep in the dataset
- upper_mon: highest moneyness value to keep in the dataset

Output: plot of the root mean squared prediction error for 3 different factor dimensions: L=2,3,4

```
<div align="center">
<img src="https://raw.githubusercontent.com/QuantLet/LETF-Moneyness/master/LETFRMSPE/rmspe234.png" alt="Image" />
</div>

