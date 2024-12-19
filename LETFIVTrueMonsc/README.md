<div style="margin: 0; padding: 0; text-align: center; border: none;">
<a href="https://quantlet.com" target="_blank" style="text-decoration: none; border: none;">
<img src="https://github.com/StefanGam/test-repo/blob/main/quantlet_design.png?raw=true" alt="Header Image" width="100%" style="margin: 0; padding: 0; display: block; border: none;" />
</a>
</div>

```
Name of Quantlet: LETFIVTrueMonsc

Published in: Leveraged ETF options implied volatility paradox

Description: 'Compute and plot the true SSO LETF option and moneyness-scaling predicted IV via the dynamic

Keywords: 'DSFM, dynamic, semiparametric, semiparametric model, pca, principal-component-analysis, factor, factor-model,

See also: LETFFactorFuncs, LETFIVSurfPlot, LETFStochLoads

Author: Sergey Nasekin

Submitted: 2016/01/21

Datafile: SPYMAT.mat, SSOMAT.mat

Input: 

- Km: 'B-spline order in moneyness direction'

- Kt: 'B-spline order in time-to-maturity direction'

- dim_mon: number of grid points for estimation in moneyness direction

- dim_ttm: number of grid points for estimation in time-to-maturity direction

- ikmon: 'parameter for setting the number of B-spline knots in moneyness direction'

- ikttm: 'parameter for setting the number of B-spline knots in time-to-maturity direction'

- tol: convergence tolerance for the Newton method

- maxiter: maximal number of iterations for the Newton method

- L: number of factor functions in the model

- beta1: leverage ratio

- c: LETF expense ratio

Output: 

- Plot 1: true SSO LETF option (red dots) and moneyness-scaling predicted IV for TTM 0.8 years on 20150520

- Plot 2: true SSO LETF option (red dots) and moneyness-scaling predicted IV for TTM 0.8 years on 20150521

```
<div align="center">
<img src="https://raw.githubusercontent.com/QuantLet/LETF-Moneyness/master/LETFIVTrueMonsc/LETFIVTrueMonsc-1.png" alt="Image" />
</div>

<div align="center">
<img src="https://raw.githubusercontent.com/QuantLet/LETF-Moneyness/master/LETFIVTrueMonsc/LETFIVTrueMonsc-2.png" alt="Image" />
</div>

