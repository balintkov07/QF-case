This folder contains a simple MATLAB demo for estimating and plotting a
GARCH model using the ES-mini futures dataset. The comments in the code
follow the notation in the case study.

Files in this folder:
- Demo.m: main script that loads the data, makes the figures, estimates the
  model, and plots the filtered volatility.
- GARCH_Filter.m: applies the GARCH recursion, returns the filtered
  variance path h_t, and computes the negative Gaussian log-likelihood.
- Logpdf_normal.m: evaluates the Gaussian log-density for one observation.
- data.mat: MATLAB data file used by Demo.m.
- data.xlsx: spreadsheet version of the distributed dataset.

Variables stored in data.mat:
- dates
- ES_open
- ES_close
- returns
- rv5_ss
- vix

Notes:
- rv5_ss stores the scaled realised variance series RV5_SS x 10^4, denoted
  rv5_ss_t in the case study.
- Demo.m estimates only the symmetric in-sample GARCH(1,1) model. It does not
  implement GJR-GARCH, real-time GARCH, real-time GJR-GARCH, HAR-RV,
  out-of-sample forecasting, or the MSE/QLIKE forecast comparison exercises.
- The optimisation uses only simple lower and upper bounds. The stationarity
  condition alpha + beta < 1 is checked inside GARCH_Filter by returning an
  infinite objective value for invalid parameter values.

To run the example, open this folder in MATLAB and run Demo.m.
