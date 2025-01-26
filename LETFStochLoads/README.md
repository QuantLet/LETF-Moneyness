<div style="margin: 0; padding: 0; text-align: center; border: none;">
<a href="https://quantlet.com" target="_blank" style="text-decoration: none; border: none;">
<img src="https://github.com/StefanGam/test-repo/blob/main/quantlet_design.png?raw=true" alt="Header Image" width="100%" style="margin: 0; padding: 0; display: block; border: none;" />
</a>
</div>

```
Name of Quantlet: LETFStochLoads

Published in: Leveraged ETF options implied volatility paradox

Description: Compute and plot the dynamics of stochastic factor loadings of the SPY LETF option implied volatility surface as well as the CBOE Volatility Index (VIX) over the period Sep. 2014 - Jul. 2015. The estimation of the dynamic semiparametric factor model with B-spline basis is performed for this purpose

Keywords: DSFM, dynamic, semiparametric, semiparametric model, pca, principal-component-analysis, factor, factor-model,

See also: LETFFactorFuncs, LETFRMSPE, LETFIVSurfPlot

Author: Sergey Nasekin

Submitted: 2017/01/17

Datafile: SPYFULL.mat, yrates1415.mat, vix.mat, zetas.mat

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

Output: plot of stochastic factor loadings of the SPY LETF option implied volatility surface and the VIX index

```
<div align="center">
<img src="https://raw.githubusercontent.com/QuantLet/LETF-Moneyness/master/LETFStochLoads/LETFStochLoads.png" alt="Image" />
</div>

