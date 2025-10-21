<div style="margin: 0; padding: 0; text-align: center; border: none;">
<a href="https://quantlet.com" target="_blank" style="text-decoration: none; border: none;">
<img src="https://github.com/StefanGam/test-repo/blob/main/quantlet_design.png?raw=true" alt="Header Image" width="100%" style="margin: 0; padding: 0; display: block; border: none;" />
</a>
</div>

```
Name of Quantlet: LETFConfBands05SSO

Published in: Leveraged ETF options implied volatility paradox

Description: Calculate and plot uniform bootstrap confidence bands for the ProShares Ultra S&P500 LETF option implied volatility at the time-to-maturity 0.5 years

Keywords: confidence-bands, bandwidth, robust estimation, kernel, implied-volatility, uniform, option, leverage effekt

See also: LETFConfBands05SPY, LETFConfBands05SDS, LETFConfBands05UPRO, LETFConfBands06SPY, LETFConfBands06SSO, LETFConfBands06UPRO, LETFConfBands06SDS, LETFConfBands07UPRO, LETFConfBands07SPY, LETFConfBands07SSO, LETFConfBands07SDS

Author: Sergey Nasekin

Submitted: 2016/01/19

Datafile: mivttmdata_05_SSO.csv, mivttmdata_05_SPY.csv

Input: 
- B: number of bootstrap iterations
- alpha: '1-confidence level of the bands'
- gridn: number of grid points for estimation
- beta:  leverage ratio of the LETF

Output: plot of bootstrap uniform confidence bands around the true curve

```
<div align="center">
<img src="https://raw.githubusercontent.com/QuantLet/LETF-Moneyness/master/LETFConfBands05SSO/sso05_bands.png" alt="Image" />
</div>

