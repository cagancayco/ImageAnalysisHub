function [model0, model1, pdfModel2, pdfModel3, outparams] = fit_logIEI (logData, identifier, varargin)
%% Fit log(IEI) logData to curves
% Usage: [model0, model1, pdfModel2, pdfModel3, outparams] = fit_logIEI (logData, identifier, varargin)
% Explanation:
%       TODO
% Side Effects:
%       Plots TODO
% Arguments:
%       logData     - a vector of log(IEI)s
%                   must be a numeric nonnegative vector
%       identifier  - an identifier for plotting logData
%                   must be a string scalar or a character vector
%       varargin    - 'XUnit': unit for IEIs
%                   must be a string scalar or a character vector
%                   default == ms
%                   - 'XLimits': x limits of histogram
%                   must be a 2-element number vector
%                   default == [min(logData) max(logData)] 
%                   - 'TruncateFlag': whether to truncate logData
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'OutFolder': directory for saving plots
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'PlotFlag': whether to plot histograms and fits
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%
% Used By:
%       /home/Matlab/Adams_Functions/ZG_fit_IEI_distributions.m
%       /media/adamX/Paula_IEIs/paula_iei4.m
%
% File History: 
% 2017-11-29 - Modified from fit_IEI.m
% 2017-11-29 - Added model0 (kernel distribution)
% 2017-11-29 - Now uses mle() instead of fistgmdist() to fit two Gaussians
% 2017-11-29 - Added lower bounds and upper bounds to mle fits
% 2017-11-30 - pdfModel1 is now a function whereas pdfPlotModel1 is a vector
% 2017-11-30 - Added spacing, threshold and void
% 2017-12-01 - threshold, void -> threshold1, void1, threshold2, void2
% 2017-12-11 - Constrain fminbnd so that mu1Model0 would not stop 
%               at a local maximum
% 2018-07-30 - Model 1 now fits truncated data based on thresholdModel0
% 2018-07-30 - Added zScoreCutoff and thresholdModel1

%% Parameters
nPoints = 500;                      % number of points for plotting pdfs
defaultLineWidth = 2;               % default line width for plots
linewidthLines = 0.5;               % line width for lines
colorHist = 'k';                    % color of histogram
lambdaInit = 0.7;                   % initial weight of Gaussian part
                                    %   in Gaussian-Exp-Exponential Fit
zScoreCutoff = 2;                   % number of standard deviations 
                                    %   away from the mean for thresholdModel1

%% Default values for optional arguments
xUnitDefault = '';                  % default unit for IEIs
xLimitsDefault = [];                % default x limits of histogram
truncateFlagDefault = false;        % whether to truncate logData by default
outFolderDefault = '';              % default directory for saving plots
plotFlagDefault = true;             % whether to plot by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = 'fit_logIEI';

% Add required inputs to the Input Parser
addRequired(iP, 'logData', ...              % a vector of IEIs
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'identifier', ...           % an identifier for plotting logData
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'XUnit', xUnitDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
addParameter(iP, 'XLimits', xLimitsDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'numel', 2}));
addParameter(iP, 'TruncateFlag', truncateFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
addParameter(iP, 'PlotFlag', plotFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, logData, identifier, varargin{:});
xUnit        = iP.Results.XUnit;
xLimits      = iP.Results.XLimits;
truncateFlag = iP.Results.TruncateFlag;
outFolder    = iP.Results.OutFolder;
plotFlag     = iP.Results.PlotFlag;

% Set dependent argument defaults
if isempty(xUnit)
    xUnit = 'log s';
    fprintf('Unit for IEIs not provided. Set to default %s!\n', xUnit);
end
if isempty(xLimits)
    xLimits = [min(logData), max(logData)];   
end
if isempty(outFolder)
    outFolder = pwd;
end

%% Fit logData
fprintf('Fitting logData for %s ... \n', identifier);
if truncateFlag
    % Eliminate logData outside of xLimits
    logData = logData(logData <= xLimits(2));
    logData = logData(logData >= xLimits(1));
end

% Total logData count and mean of logData
% Calculate statistics from logData
npts = length(logData);
logMeanIEI = log(mean(exp(logData)));
logMedianIEI = log(median(exp(logData)));
minLogIEI = min(logData);
maxLogIEI = max(logData);
meanLogIEI = mean(logData);
stdLogIEI = std(logData);

% Number of bins is square root of npts logData count
%   Cocatre-Zilgien & Delcomyn, 1992 
%   ("Identification of bursts in spike trains")
% nBins = round(sqrt(npts));
nBins = round(1.87 * (npts - 1) ^ 0.4);

% Compute histogram edges, area, and x values for fit pdf
edges = linspace(xLimits(1), xLimits(2), nBins);
binWidth = xLimits(2)/nBins;        % bin width
harea = binWidth * npts;            % total histogram area
x = linspace(xLimits(1), xLimits(2), nPoints)';

%% Fit logData to a kernel distribution using Gaussian kernels
%   This yielded many small peaks that were probably insignificant
%   model0 = fitdist(logData, 'Kernel');

%% Fit logData to a kernel distribution using Epanechnikov (parabolic) kernels
%   This was not much different from Gaussian kernels
% model0 = fitdist(logData, 'Kernel', 'Kernel', 'epanechnikov');

%% Fit logData to a kernel distribution using Gaussian kernels
%   and using a fifth of the standard deviation as the bandwidth
model0 = fitdist(logData, 'Kernel', 'Bandwidth', stdLogIEI/5);

% Get pdf and scaled pdf for the fit
pdfModel0 = @(X) pdf(model0, X);
pdfPlotModel0 = harea * pdfModel0(x);

% Find right most two peaks of the kernel distribution
%   This follows the algorithm of Selinger et al., 2007
[~, locs] = findpeaks(pdfPlotModel0);
if length(locs) > 0
    peak2Model0 = x(locs(end));
else
    peak2Model0 = maxLogIEI;
end
if length(locs) > 1
    peak1Model0 = x(locs(end-1));
else
    peak1Model0 = NaN;
end

% Find the mode of the kernel distribution
[~, idxModeApprox] = max(pdfPlotModel0);
modeModel0Approx = x(idxModeApprox);
skewness = abs(modeModel0Approx - meanLogIEI);
leftBnd = max(minLogIEI, modeModel0Approx - skewness);
rightBnd = min(maxLogIEI, modeModel0Approx + skewness);
mu1Model0 = fminbnd(@(X) -pdfModel0(X), leftBnd, rightBnd);

% Find all peaks to the right of the mode of the kernel distribution
toTheRight = x >= mu1Model0;                % whether the index of x is 
                                            % to the right of the mode
xToTheRight = x(toTheRight);                % parts of x that is 
                                            % to the right of the mode
pdfToTheRight = pdfPlotModel0(toTheRight);  % parts of the pdf that is 
                                            % to the right of the mode
[pksToTheRight, locsToTheRight] = findpeaks(pdfToTheRight);

% Select the peak that gives the largest void parameter 
%   when paired with the mode
%   This follows the algorithm of Pasquale et al., 2010
nToTheRight = length(locsToTheRight);
if nToTheRight == 0
    mu2Model0 = NaN;
    spacingModel0 = 0;
    thresholdModel0 = NaN;
    voidModel0 = 0;
else
    thresholds = zeros(nToTheRight, 1);     % stores all possible thresholds
    voids = zeros(nToTheRight, 1);          % stores all void parameters
    for iToTheRight = 1:nToTheRight
        % The threshold is the minimum between the mode and the current peak
        thresholds(iToTheRight) = ...
            fminbnd(pdfModel0, mu1Model0, ...
                    xToTheRight(locsToTheRight(iToTheRight)));

        % Compute the void parameter corresponding to this pair of peaks
        voids(iToTheRight) = ...
            1 - pdfModel0(thresholds(iToTheRight)) / ...
                sqrt(pdfModel0(mu1Model0) * pksToTheRight(iToTheRight));
    end
    % Select the peak that gives the largest void parameter
    [voidModel0, iSelected] = max(voids);   

    % Record the x value and threshold for this peak
    mu2Model0 = xToTheRight(locsToTheRight(iSelected));
    thresholdModel0 = thresholds(iSelected);    

    % Compute the spacing parameter for this peak
    spacingModel0 = mu2Model0 - mu1Model0;
end

%% Fit truncated logData to one Gaussian
% Truncate logData with thresholdModel0
if ~isnan(thresholdModel0)
    logDataTruncated = logData(logData < thresholdModel0);
else
    logDataTruncated = logData;
end

% Fit to a single Gaussian
model1 = fitgmdist(logDataTruncated, 1);

% Find mean and standard deviations for the fits
muModel1 = model1.mu;
sigmaModel1 = model1.Sigma;

% Get pdf and scaled pdf for the fit
pdfModel1 = @(X) pdf(model1, X);
pdfPlotModel1 = harea * pdfModel1(x);

% Find mean and standard deviations for the fits
thresholdModel1 = muModel1 + zScoreCutoff * sigmaModel1;

%% Fit logData to two Gaussians

%{
model2 = fitgmdist(logData, 2);

% Find mean and standard deviations for the fits
[~, origInd] = sort(model2.mu, 'ascend');
mu1Model2 = model2.mu(origInd(1));
mu2Model2 = model2.mu(origInd(2));
sigma1Model2 = model2.Sigma(origInd(1));
sigma2Model2 = model2.Sigma(origInd(2));
prop1Model2 = model2.ComponentProportion(origInd(1));
prop2Model2 = model2.ComponentProportion(origInd(2));

% Get scaled pdf for the fit
pdfModel2 = @(X) pdf(model2, X);
pdfPlotModel2 = harea * pdfModel2(x);

% Get scaled pdfs for components of the mixture fit
comp1Model2 = @(X) prop1Model2 * normpdf(X, mu1Model2, sigma1Model2);
comp1PlotModel2 = harea * comp1Model2(x);
comp2Model2 = @(X) prop2Model2 * normpdf(X, mu2Model2, sigma2Model2);
comp2PlotModel2 = harea * comp2Model2(x);
%}

mypdf2 = @(X, lambda, mu1, sigma1, mu2, sigma2) ...
                lambda * normpdf(X, mu1, sigma1) + ...
                (1 - lambda) * normpdf(X, mu2, sigma2);
lowerBound2 = [0, minLogIEI, 0.0000001, logMedianIEI, 0.0000001];
start2 = [lambdaInit, logMedianIEI, stdLogIEI, ...
            max(logMedianIEI, peak2Model0), stdLogIEI];
upperBound2 = [1, maxLogIEI, Inf, maxLogIEI, Inf];
[phat2, pci2] = mle(logData, 'pdf', mypdf2, 'start', start2, ...
                    'LowerBound', lowerBound2, 'UpperBound', upperBound2);

% Find parameters for the fits; make sure mu1 <= mu2
if phat2(2) <= phat2(4)
    lambdaModel2 = phat2(1);
    lambdaCIModel2 = pci2(:, 1);
    mu1Model2 = phat2(2);
    mu1CIModel2 = pci2(:, 2);
    sigma1Model2 = phat2(3);
    sigma1CIModel2 = pci2(:, 3);
    mu2Model2 = phat2(4);
    mu2CIModel2 = pci2(:, 4);
    sigma2Model2 = phat2(5);
    sigma2CIModel2 = pci2(:, 5);
else
    lambdaModel2 = 1 - phat2(1);
    lambdaCIModel2 = flipud(1 - pci2(:, 1));
    mu1Model2 = phat2(4);
    mu1CIModel2 = pci2(:, 4);
    sigma1Model2 = phat2(5);
    sigma1CIModel2 = pci2(:, 5);
    mu2Model2 = phat2(2);
    mu2CIModel2 = pci2(:, 2);
    sigma2Model2 = phat2(3);
    sigma2CIModel2 = pci2(:, 3);
end

% Get model pdf
pdfModel2 = @(X) mypdf2(X, lambdaModel2, mu1Model2, sigma1Model2, ...
                        mu2Model2, sigma2Model2);

% Get scaled pdf for the fit
pdfPlotModel2 = harea * pdfModel2(x);

% Get scaled pdfs for components of the mixture fit
comp1Model2 = @(X) lambdaModel2 * normpdf(X, mu1Model2, sigma1Model2);
comp1PlotModel2 = harea * comp1Model2(x);
comp2Model2 = @(X) (1 - lambdaModel2) * normpdf(X, mu2Model2, sigma2Model2);
comp2PlotModel2 = harea * comp2Model2(x);

% Compute the "time scale spacing parameter",
%   see Selinger et al., 2007
spacingModel2 = mu2Model2 - mu1Model2;
spacingErrModel2 = (mu2CIModel2(2) - mu2Model2) + ...
                    (mu1CIModel2(2) - mu1Model2);
spacingCIModel2 = [spacingModel2 - spacingErrModel2; ...
                    spacingModel2 + spacingErrModel2];

% Compute threshold #1 that separates the two components, 
%   defined as the intersection between the pdfs of each component
diffCompModel2 = @(X) comp1Model2(X) - comp2Model2(X);
if diffCompModel2(mu1Model2) * diffCompModel2(mu2Model2) < 0
    threshold1Model2 = fzero(diffCompModel2, [mu1Model2, mu2Model2]);
elseif diffCompModel2(mu1Model2) * ...
            diffCompModel2(mu2Model2 + 2 * sigma2Model2) < 0
    threshold1Model2 = fzero(diffCompModel2, ...
                        [mu1Model2, mu2Model2 + 2 * sigma2Model2]);
elseif diffCompModel2(mu1Model2) * diffCompModel2(maxLogIEI) < 0
    threshold1Model2 = fzero(diffCompModel2, [mu1Model2, maxLogIEI]);
else
    threshold1Model2 = NaN;
end
% TODO: threshold1CIModel2

% Compute threshold #2 that separates the two components, 
%   defined as the minimum between the peaks in the overall pdf
threshold2Model2 = fminbnd(pdfModel2, mu1Model2, mu2Model2);
% TODO: threshold2CIModel2

% Compute the "void parameter" #1 or #2 using threshold #1 or #2 
%   as the minimum in the equation in Selinger et al., 2007
void1Model2 = 1 - pdfModel2(threshold1Model2)/ ...
                    sqrt(pdfModel2(mu1Model2) * pdfModel2(mu2Model2));
% TODO: void1CIModel2
void2Model2 = 1 - pdfModel2(threshold2Model2)/ ...
                    sqrt(pdfModel2(mu1Model2) * pdfModel2(mu2Model2));
% TODO: void2CIModel2

%% Fit logData to a mixture of Gaussian and Exp-Exponential
mypdf3 = @(X, lambda, mu1, sigma, mu2) ...
                lambda * normpdf(X, mu1, sigma) + ...
                (1 - lambda) * (1/mu2) * exp(X - exp(X)/mu2);
lowerBound3 = [0, minLogIEI, 0.0000001, exp(logMedianIEI)];
upperBound3 = [1, maxLogIEI, Inf, exp(maxLogIEI)];
try
    start3 = [lambdaInit, mu1Model2, sigma1Model2, exp(mu2Model2)];
    [phat3, pci3] = mle(logData, 'pdf', mypdf3, 'start', start3, ...
                    'LowerBound', lowerBound3, 'UpperBound', upperBound3);
catch
    start3 = [lambdaInit, logMedianIEI, stdLogIEI, exp(peak2Model0)];
    [phat3, pci3] = mle(logData, 'pdf', mypdf3, 'start', start3, ...
                    'LowerBound', lowerBound3, 'UpperBound', upperBound3);
end

% Find parameters for the fits
lambdaModel3 = phat3(1);
lambdaCIModel3 = pci3(:, 1);
mu1Model3 = phat3(2);
mu1CIModel3 = pci3(:, 2);
sigmaModel3 = phat3(3);
sigmaCIModel3 = pci3(:, 3);
mu2Model3 = phat3(4);
mu2CIModel3 = pci3(:, 4);
pdfModel3 = @(X) mypdf3(X, lambdaModel3, mu1Model3, sigmaModel3, mu2Model3);

% Get scaled pdf for the fit
pdfPlotModel3 = harea * pdfModel3(x);

% Get scaled pdfs for components of the mixture fit
comp1Model3 = @(X) lambdaModel3 * normpdf(X, mu1Model3, sigmaModel3);
comp1PlotModel3 = harea * comp1Model3(x);
comp2Model3 = @(X) (1 - lambdaModel3) * ...
                    (1/mu2Model3) * exp(X - exp(X)/mu2Model3);
comp2PlotModel3 = harea * comp2Model3(x);

% Compute the "time scale spacing parameter",
%   see Selinger et al., 2007
spacingModel3 = log(mu2Model3) - mu1Model3;
spacingErrModel3 = (mu2CIModel3(2) - mu2Model3)/mu2Model3 + ...
                    (mu1CIModel3(2) - mu1Model3);
spacingCIModel3 = [spacingModel3 - spacingErrModel3; ...
                    spacingModel3 + spacingErrModel3];

% Compute threshold #1 that separates the two components, 
%   defined as the intersection between the pdfs of each component
diffCompModel3 = @(X) comp1Model3(X) - comp2Model3(X);
if diffCompModel3(mu1Model3) * diffCompModel3(log(mu2Model3)) < 0
    threshold1Model3 = fzero(diffCompModel3, [mu1Model3, log(mu2Model3)]);
elseif diffCompModel3(mu1Model3) * ...
        diffCompModel3(log(mu2Model3) + 2 * sigma2Model2) < 0
    threshold1Model3 = fzero(diffCompModel3, ...
                [mu1Model3, log(mu2Model3) + 2 * sigma2Model2]);
elseif diffCompModel3(mu1Model3) * diffCompModel3(maxLogIEI) < 0
    threshold1Model3 = fzero(diffCompModel3, [mu1Model3, maxLogIEI]);
else
    threshold1Model3 = NaN;
end
% TODO: threshold1CIModel3

% Compute threshold #2 that separates the two components, 
%   defined as the minimum between the peaks in the overall pdf
%   or the minimum after mu1
if log(mu2Model3) >= mu1Model3
    threshold2Model3 = fminbnd(pdfModel3, mu1Model3, log(mu2Model3));
else
    threshold2Model3 = fminbnd(pdfModel3, mu1Model3, maxLogIEI);
end

% TODO: threshold2CIModel3

% Compute the "void parameter" #1 or #2 using threshold #1 or #2 
%   as the minimum in the equation in Selinger et al., 2007
void1Model3 = 1 - pdfModel3(threshold1Model3)/ ...
                    sqrt(pdfModel3(mu1Model3) * pdfModel3(log(mu2Model3)));
% TODO: void1CIModel3
void2Model3 = 1 - pdfModel3(threshold2Model3)/ ...
                    sqrt(pdfModel3(mu1Model3) * pdfModel3(log(mu2Model3)));
% TODO: void2CIModel3

%% Plot histograms and fits
if plotFlag
    fprintf('Plotting logData for %s ... \n', identifier);

    % Change default values
    set(groot, 'defaultLineLineWidth', defaultLineWidth);

    % Plot histograms with fits
    h = figure(100);
    clf(h);
    hold on;
    histogram(logData, edges, 'FaceColor', colorHist);
    xlim(xLimits);
    yLimits = get(gca, 'YLim');
    line(logMeanIEI * ones(1, 2), yLimits, 'Color', 'g', 'LineStyle', '--', ...
            'LineWidth', linewidthLines, ...
            'DisplayName', ['logMean = ', num2str(logMeanIEI, 3)]);
    line(logMedianIEI * ones(1, 2), yLimits, 'Color', 'y', 'LineStyle', '--', ...
            'LineWidth', linewidthLines, ...
            'DisplayName', ['logMedian = ', num2str(logMedianIEI, 3)]);
    legend('location', 'northeast');
    xlabel(['Log of Inter-event intervals (', xUnit, ')']);
    ylabel('Count');
    title(['Raw logData for ', strrep(identifier, '_', '\_')]);
    saveas(h, fullfile(outFolder, [identifier, '_HistOnly']), 'png');

    if ~isempty(model0)       
        % Plot histograms with Kernel distribution
        h = figure(101);
        clf(h);
        hold on;
        histogram(logData, edges, 'FaceColor', colorHist);
        xlim(xLimits);
        plot(x, pdfPlotModel0, 'c', 'Displayname', 'kernel fit');
        yLimits = get(gca, 'YLim');
        line(logMeanIEI * ones(1, 2), yLimits, 'Color', 'g', 'LineStyle', '--', ...
                'LineWidth', linewidthLines, ...
                'DisplayName', ['logMean = ', num2str(logMeanIEI, 3)]);
        line(logMedianIEI * ones(1, 2), yLimits, 'Color', 'y', 'LineStyle', '--', ...
                'LineWidth', linewidthLines, ...
                'DisplayName', ['logMedian = ', num2str(logMedianIEI, 3)]);
        line(peak1Model0 * ones(1, 2), yLimits, 'Color', 'm', 'LineStyle', '--', ...
                'LineWidth', linewidthLines, ...
                'DisplayName', ['2nd right most peak = ', num2str(peak1Model0, 3)]);
        line(peak2Model0 * ones(1, 2), yLimits, 'Color', 'r', 'LineStyle', '--', ...
                'LineWidth', linewidthLines, ...
                'DisplayName', ['Right most peak = ', num2str(peak2Model0, 3)]);
        line(mu1Model0 * ones(1, 2), yLimits, 'Color', 'm', 'LineStyle', '-', ...
                'LineWidth', linewidthLines, ...
                'DisplayName', ['Mode = ', num2str(mu1Model0, 3)]);
        line(mu2Model0 * ones(1, 2), yLimits, 'Color', 'r', 'LineStyle', '-', ...
                'LineWidth', linewidthLines, ...
                'DisplayName', ['Best peak = ', num2str(mu2Model0, 3)]);
        line(thresholdModel0 * ones(1, 2), yLimits, 'Color', 'b', 'LineStyle', '--', ...
                'LineWidth', linewidthLines, ...
                'DisplayName', ['Threshold = ', num2str(thresholdModel0, 3)]);
        text(0.05, 0.975, ...
                ['Spacing Parameter = ', num2str(spacingModel0, 3)], ...
                'Units', 'normalized')
        text(0.05, 0.95, ...
                ['Void Parameter = ', num2str(voidModel0, 3)], ...
                'Units', 'normalized')
        legend('location', 'northeast');
        xlabel(['Log of Inter-event intervals (', xUnit, ')']);
        ylabel('Count');
        title(['Kernel Fit for ', strrep(identifier, '_', '\_')]);
        saveas(h, fullfile(outFolder, [identifier, '_Kernel']), 'png');
    end

    if ~isempty(model1) || ~isempty(pdfModel2)
        % Plot histograms with Gaussian fits
        h = figure(102);
        clf(h);
        hold on;
        histogram(logData, edges, 'FaceColor', colorHist);
        xlim(xLimits);
        if ~isempty(model1)
            plot(x, pdfPlotModel1, 'g', 'Displayname', 'single Gaussian');
        end
        if ~isempty(pdfModel2)
            plot(x, pdfPlotModel2, 'c', 'Displayname', 'double Gaussian');
            plot(x, comp1PlotModel2, 'm', 'Displayname', 'Component #1 of double Gaussian');
            plot(x, comp2PlotModel2, 'r', 'Displayname', 'Component #2 of double Gaussian');
        end
        yLimits = get(gca, 'YLim');
        if ~isempty(model1)
            line(muModel1 * ones(1, 2), yLimits, 'Color', 'g', 'LineStyle', '--', ...
                    'LineWidth', linewidthLines, ...
                    'DisplayName', ['Mean of Single Gaussian = ', num2str(muModel1, 3)]);
        end
        if ~isempty(pdfModel2)
            line(mu1Model2 * ones(1, 2), yLimits, 'Color', 'm', 'LineStyle', '-', ...
                    'LineWidth', linewidthLines, ...
                    'DisplayName', ['Mean of Component #1 = ', num2str(mu1Model2, 3)]);
            line(mu2Model2 * ones(1, 2), yLimits, 'Color', 'r', 'LineStyle', '-', ...
                    'LineWidth', linewidthLines, ...
                    'DisplayName', ['Mean of Component #2 = ', num2str(mu2Model2, 3)]);
            line(threshold1Model2 * ones(1, 2), yLimits, 'Color', 'c', 'LineStyle', '--', ...
                    'LineWidth', linewidthLines, ...
                    'DisplayName', ['Threshold #1 = ', num2str(threshold1Model2, 3)]);
            line(threshold2Model2 * ones(1, 2), yLimits, 'Color', 'b', 'LineStyle', '--', ...
                    'LineWidth', linewidthLines, ...
                    'DisplayName', ['Threshold #2 = ', num2str(threshold2Model2, 3)]);
            text(0.05, 0.975, ...
                    ['Spacing Parameter = ', num2str(spacingModel2, 3)], ...
                    'Units', 'normalized')
            text(0.05, 0.95, ...
                    ['Void Parameter #1 = ', num2str(void1Model2, 3)], ...
                    'Units', 'normalized')
            text(0.05, 0.925, ...
                    ['Void Parameter #2 = ', num2str(void2Model2, 3)], ...
                    'Units', 'normalized')
        end
        legend('location', 'northeast');
        xlabel(['Log of Inter-event intervals (', xUnit, ')']);
        ylabel('Count');
        title(['Gaussian Fits for ', strrep(identifier, '_', '\_')]);
        saveas(h, fullfile(outFolder, [identifier, '_GaussOnly']), 'png');
    end
    
    if ~isempty(pdfModel3)
        % Plot histograms with Gaussian-Exp-Exponential fits
        h = figure(103);
        clf(h);
        hold on;
        histogram(logData, edges, 'FaceColor', colorHist);
        xlim(xLimits);
        plot(x, pdfPlotModel3, 'c', 'Displayname', 'Gaussian + Exp-Exponential');
        plot(x, comp1PlotModel3, 'm', 'Displayname', 'Gaussian part');
        plot(x, comp2PlotModel3, 'r', 'Displayname', 'Exp-Exponential part');
        yLimits = get(gca, 'YLim');
        line(mu1Model3 * ones(1, 2), yLimits, 'Color', 'm', 'LineStyle', '-', ...
                'LineWidth', linewidthLines, ...
                'DisplayName', ['Mean of Gaussian part = ', num2str(mu1Model3, 3)]);
        line(log(mu2Model3) * ones(1, 2), yLimits, 'Color', 'r', 'LineStyle', '-', ...
                'LineWidth', linewidthLines, ...
                'DisplayName', ['log(mu) of Exp-Exponential part = ', num2str(log(mu2Model3), 3)]);
        line(threshold1Model3 * ones(1, 2), yLimits, 'Color', 'c', 'LineStyle', '--', ...
                'LineWidth', linewidthLines, ...
                'DisplayName', ['Threshold #1 = ', num2str(threshold1Model3, 3)]);
        line(threshold2Model3 * ones(1, 2), yLimits, 'Color', 'b', 'LineStyle', '--', ...
                'LineWidth', linewidthLines, ...
                'DisplayName', ['Threshold #2 = ', num2str(threshold2Model3, 3)]);
        text(0.05, 0.975, ...
                ['Spacing Parameter = ', num2str(spacingModel3, 3)], ...
                'Units', 'normalized')
        text(0.05, 0.95, ...
                ['Void Parameter #1 = ', num2str(void1Model3, 3)], ...
                'Units', 'normalized')
        text(0.05, 0.925, ...
                ['Void Parameter #2 = ', num2str(void2Model3, 3)], ...
                'Units', 'normalized')
        legend('location', 'northeast');
        xlabel(['Log of Inter-event intervals (', xUnit, ')']);
        ylabel('Count');
        title(['Gaussian-Exp-Exponential Fit for ', strrep(identifier, '_', '\_')]);
        saveas(h, fullfile(outFolder, [identifier, '_GaussExp']), 'png');
    end

    % Remove default values
    set(groot, 'defaultLineLineWidth', 'remove');
end

%% Package parameters into outparams
outparams.logMeanIEI = logMeanIEI;
outparams.logMedianIEI = logMedianIEI;

if ~isempty(model0)
    outparams.model0 = model0;
    outparams.pdfModel0 = pdfModel0;
    outparams.peak1Model0 = peak1Model0;
    outparams.peak2Model0 = peak2Model0;
    outparams.mu1Model0 = mu1Model0;
    outparams.mu2Model0 = mu2Model0;
    outparams.spacingModel0 = spacingModel0;
    outparams.thresholdModel0 = thresholdModel0;
    outparams.voidModel0 = voidModel0;
end

if ~isempty(model1)
    outparams.model1 = model1;
    outparams.pdfModel1 = pdfModel1;
    outparams.muModel1 = muModel1;
    outparams.sigmaModel1 = sigmaModel1;
    outparams.thresholdModel1 = thresholdModel1;
end

if ~isempty(pdfModel2)
    outparams.mypdf2 = mypdf2;
    outparams.pdfModel2 = pdfModel2;
    outparams.comp1Model2 = comp1Model2;
    outparams.comp2Model2 = comp2Model2;
    outparams.lambdaModel2 = lambdaModel2;
    outparams.lambdaCIModel2 = lambdaCIModel2;
    outparams.mu1Model2 = mu1Model2;
    outparams.mu1CIModel2 = mu1CIModel2;
    outparams.sigma1Model2 = sigma1Model2;
    outparams.sigma1CIModel2 = sigma1CIModel2;
    outparams.mu2Model2 = mu2Model2;
    outparams.mu2CIModel2 = mu2CIModel2;
    outparams.sigma2Model2 = sigma2Model2;
    outparams.sigma2CIModel2 = sigma2CIModel2;
    outparams.spacingModel2 = spacingModel2;
    outparams.spacingCIModel2 = spacingCIModel2;
    outparams.threshold1Model2 = threshold1Model2;
%    outparams.threshold1CIModel2 = threshold1CIModel2;
    outparams.threshold2Model2 = threshold2Model2;
%    outparams.threshold2CIModel2 = threshold2CIModel2;
    outparams.void1Model2 = void1Model2;
%    outparams.void1CIModel2 = void1CIModel2;
    outparams.void2Model2 = void2Model2;
%    outparams.void2CIModel2 = void2CIModel2;

%    outparams.mu2Model2 = mu2Model2;
%    outparams.sigma1Model2 = sigma1Model2;
%    outparams.sigma2Model2 = sigma2Model2;
%    outparams.prop1Model2 = prop1Model2;
%    outparams.prop2Model2 = prop2Model2;
end

if ~isempty(pdfModel3)
    outparams.mypdf3 = mypdf3;
    outparams.pdfModel3 = pdfModel3;
    outparams.comp1Model3 = comp1Model3;
    outparams.comp2Model3 = comp2Model3;
    outparams.lambdaModel3 = lambdaModel3;
    outparams.lambdaCIModel3 = lambdaCIModel3;
    outparams.mu1Model3 = mu1Model3;
    outparams.mu1CIModel3 = mu1CIModel3;
    outparams.sigmaModel3 = sigmaModel3;
    outparams.sigmaCIModel3 = sigmaCIModel3;
    outparams.mu2Model3 = mu2Model3;
    outparams.mu2CIModel3 = mu2CIModel3;
    outparams.spacingModel3 = spacingModel3;
    outparams.spacingCIModel3 = spacingCIModel3;
    outparams.threshold1Model3 = threshold1Model3;
%    outparams.threshold1CIModel3 = threshold1CIModel3;
    outparams.threshold2Model3 = threshold2Model3;
%    outparams.threshold2CIModel3 = threshold2CIModel3;
    outparams.void1Model3 = void1Model3;
%    outparams.void1CIModel3 = void1CIModel3;
    outparams.void2Model3 = void2Model3;
%    outparams.void2CIModel3 = void2CIModel3;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:
    if exist('mu1Model2', 'var') == 1

skipModel3Default = false;          % whether to skip model 3 by default
addParameter(iP, 'SkipModel3', skipModel3Default, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
skipModel3   = iP.Results.SkipModel3;

try
catch
    model0 = [];
end

try
catch
    model1 = [];
end

try
catch
    pdfModel3 = [];
end

% Use the mean plus 2 standard deviations as the threshold
thresholdModel1 = muModel1 + 2 * sigmaModel1

%% Fit logData to one Gaussian
model1 = fitgmdist(logData, 1);

% Find mean and standard deviations for the fits
muModel1 = model1.mu;
sigmaModel1 = model1.Sigma;

% Get pdf and scaled pdf for the fit
pdfModel1 = @(X) pdf(model1, X);
pdfPlotModel1 = harea * pdfModel1(x);

%}
