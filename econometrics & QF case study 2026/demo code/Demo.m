%% Start with a clean Matlab session
clear
close all

%% Figure style
% Force all figures to use a white background and black axes/text.
% This keeps exported figures readable in the case-study handout.
set(groot, 'defaultFigureColor', 'w');
set(groot, 'defaultAxesColor', 'w');
set(groot, 'defaultAxesXColor', 'k');
set(groot, 'defaultAxesYColor', 'k');
set(groot, 'defaultAxesZColor', 'k');
set(groot, 'defaultTextColor', 'k');

%% Notation used in this demo
% The comments below follow the notation in the case study.
% returns = r_t, the close-to-close return in percent.
% h       = h_t, the conditional variance in the GARCH model.
% rv5_ss  = RV5_SS x 10^4, the scaled realised-variance target.
% vix     = VIX_t, used to construct the rescaled VIX benchmark forecast.

%% Find the folder that contains this script
% This allows the code to work even if you move the folder to a new location.
file_name_full   = matlab.desktop.editor.getActiveFilename;
separator_id     = strfind(file_name_full,'/');
file_name_short  = file_name_full(max(separator_id)+1:end);
folder_name      = erase(file_name_full,file_name_short);

%% Make the script folder the current working folder
cd(folder_name)
disp(folder_name)

%% Load the prepared dataset
% The distributed spreadsheet is data.xlsx. The Matlab file data.mat contains
% the same variables needed in this demo and loads faster than the spreadsheet.
load('data')

%% Remove observations with missing values
% For example, the first close-to-close return is often missing.
idx_keep   = isfinite(returns) & isfinite(rv5_ss) & isfinite(vix);
dates      = dates(idx_keep,:);
returns    = returns(idx_keep,:);
rv5_ss     = rv5_ss(idx_keep,:);
vix        = vix(idx_keep,:);
if exist('ES_open','var'),  ES_open  = ES_open(idx_keep,:);  end
if exist('ES_close','var'), ES_close = ES_close(idx_keep,:); end

%% Convert the dates to Matlab datetime format
dates_str = cellstr(num2str(dates));
dates     = datetime(dates_str, 'InputFormat', 'yyyyMMdd');

%% Define x-axis ticks and limits for the available sample
tick_dates = datetime((year(dates(1)):year(dates(end))+1)',1,1);
x_axis_limits = [dates(1), dates(end)];

%% Plot the daily returns
figure('Color', 'w')
plot(dates, returns, 'k')
ylabel('$r_t$', 'Interpreter', 'latex')
hold on

% Format the x-axis
xticks(tick_dates);
xtickformat('yy');
xtickangle(0);

% Set axis limits and style
xlim(x_axis_limits);
ylim([-15, 15]);
set(gca, 'FontName', 'Times', 'FontSize', 20, 'TickDir', 'out', ...
    'Box', 'off', 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1);

%% Plot the VIX, corresponding to the left panel of Figure 1 in the case study
figure('Color', 'w')
plot(dates, vix, 'k')

% Format the x-axis
xticks(tick_dates);
xtickformat('yy');
xtickangle(0);

% Set axis limits and style
xlim(x_axis_limits);
ylim([0, 100]);
set(gca, 'FontName', 'Times', 'FontSize', 20, 'TickDir', 'out', ...
    'Box', 'off', 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1);
title('VIX (daily close)')

%% Plot the VIX-change/return relation, corresponding to the right panel of Figure 1
figure('Color', 'w')
scatter(diff(vix), returns(2:end), 'ok')

% Set axis style
set(gca, 'FontName', 'Times', 'FontSize', 20, 'TickDir', 'out', ...
    'Box', 'off', 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1);
xlabel('$\Delta \mbox{VIX}_t$', 'Interpreter', 'latex')
ylabel('$r_t$ (E-mini futures)', 'Interpreter', 'latex')

%% Estimate a simple GARCH(1,1) model
% This is the basic GARCH recursion in the case study:
% h_{t+1} = omega + alpha * (r_t - mu)^2 + beta * h_t.
% This demo estimates only the symmetric GARCH(1,1) model. It does not implement
% the GJR-GARCH, real-time GARCH, HAR-RV, out-of-sample forecast, or loss-comparison
% exercises from the full case study.
format short
clear GARCH_Filter

%% Choose starting values for the optimisation
% Order of parameters: [mu, omega, alpha, beta]
startingvalues = [mean(returns); 0.04; 0.10; 0.80];

%% Evaluate the negative log-likelihood at the starting values
% A good optimiser should find a lower value than this.
GARCH_Filter(startingvalues, returns)

%% Set optimisation options
clearvars options
options  =  optimset('fmincon');
options  =  optimset(options , 'TolFun'      , 1e-6);
options  =  optimset(options , 'TolX'        , 1e-6);
options  =  optimset(options , 'Display'     , 'on');
options  =  optimset(options , 'Diagnostics' , 'on');
options  =  optimset(options , 'LargeScale'  , 'off');
options  =  optimset(options , 'MaxFunEvals' , 10^6);
options  =  optimset(options , 'MaxIter'     , 10^6);

%% Set lower and upper bounds for the parameters
% Order of parameters: [mu, omega, alpha, beta]
% The condition alpha + beta < 1 is checked inside GARCH_Filter.
lowerbound = [-inf, 1e-8, 0, 0];
upperbound = [ inf,  inf, 1, 1];

%% Run the optimisation
% This should be very fast for the basic GARCH model.
tic
[ML_parameters, NegativeLogLikelihood1] = fmincon(@(parameter_vector) ...
    GARCH_Filter(parameter_vector, returns), ...
    startingvalues, [], [], [], [], lowerbound, upperbound, [], options);
toc
disp(ML_parameters)

%% Check the estimated GARCH persistence
alpha_ML = ML_parameters(3);
beta_ML  = ML_parameters(\nGARCH persistence alpha + beta = %.4f\n', alpha_ML + beta_ML);

if alpha_ML + beta_ML < 1
    fprintf('The estimated GARCH model satisfies alpha + beta < 1.\n');
else
    fprintf('The estimated GARCH model does not satisfy alpha + beta < 1.\n');
end
4);

fprintf('
%% Plot the filtered volatility path
% GARCH_Filter returns h_t, the filtered conditional variance.
[~, h] = GARCH_Filter(ML_parameters, returns);

figure('Color', 'w')
plot(dates, returns)
hold on
plot(dates, sqrt(h), 'r')

% Format the x-axis
xticks(tick_dates);
xtickformat('yy');
xtickangle(0);

% Set axis limits and style
xlim(x_axis_limits);
ylim([-15, 15]);
set(gca, 'FontName', 'Times', 'FontSize', 20, 'TickDir', 'out', ...
    'Box', 'off', 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1);
leg = legend('ES-mini futures returns', 'Estimated $\sqrt{h_t}$ from GARCH model', 'Interpreter', 'latex');
set(leg, 'TextColor', 'k', 'Color', 'w', 'EdgeColor', 'k');

%% Check whether the filtered shocks look approximately standard normal
% The model implies eps_t = (r_t - mu) / sqrt(h_t). Under the model, eps_t
% should have mean close to 0 and variance close to 1.
mu_ML           = ML_parameters(1);
implied_epsilon = (returns - mu_ML) ./ sqrt(h);
disp([mean(implied_epsilon); var(implied_epsilon)])

%% Compare different volatility measures
% This graph shows four series used in the case study:
% 1. squared ES-mini returns, one target variable;
% 2. rv5_ss_t = RV5_SS x 10^4, the realised-variance target;
% 3. VIX_t^2/250, the rescaled VIX benchmark forecast;
% 4. h_t, the filtered GARCH variance from this demo.
figure('Color', 'w')
plot(dates, returns.^2, 'ok')
hold on
plot(dates, rv5_ss, 'b', 'LineWidth', 3)
plot(dates, (1/250) * vix.^2, 'y', 'LineWidth', 3) % VIX_t^2/250 benchmark
plot(dates, h, 'r:', 'LineWidth', 2)

% Format the x-axis
xticks(tick_dates);
xtickformat('yy');
xtickangle(0);

% Set axis limits and style
xlim(x_axis_limits);
ylim([0, 100]);
set(gca, 'FontName', 'Times', 'FontSize', 20, 'TickDir', 'out', ...
    'Box', 'off', 'Color', 'w', 'XColor', 'k', 'YColor', 'k', 'LineWidth', 1);

% Add a legend and title
leg = legend('Squared ES-mini returns', 'RV5\_SS $\times 10^4$', 'Rescaled version of VIX$^2$', ...
    'Estimated $h_t$ from GARCH model', 'Interpreter', 'latex');
set(leg, 'TextColor', 'k', 'Color', 'w', 'EdgeColor', 'k');
title('Variance of ES-mini futures returns by different measures')

%% Compare the sample means of squared returns and realised variance
% Both are on the squared-percent scale because returns are measured in percent.
disp([mean(returns.^2), mean(rv5_ss)])


