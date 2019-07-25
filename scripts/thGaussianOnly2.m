function [thMinimum, thIntersection, mu1ModelGaussian, sigma1ModelGaussian, mu2ModelGaussian, plotGaussian] = thGaussianOnly2(periods, logTransformed, plotTh)

if ~logTransformed
    medianPeriods = median(periods);
else
    medianPeriods = log(median(exp(periods)));
end

stdPeriods = std(periods);
minPeriods = min(periods);
maxPeriods = max(periods);
npts = length(periods);
nBins = round(1.87 * (npts - 1) ^ 0.4);
edges = linspace(minPeriods, maxPeriods, nBins);
binWidth = max(periods)/nBins;
harea = binWidth * npts;
x = linspace(minPeriods, maxPeriods, 500);

modelKernel   = fitdist(periods, 'Kernel', 'Bandwidth', stdPeriods/5);
pdfModelKernel = @(X) pdf(modelKernel, X);
pdfPlotModelKernel = harea * pdfModelKernel(x);

[~, locs] = findpeaks(pdfPlotModelKernel);
if ~isempty(locs)
    peakEndModelKernel = x(locs(end));
else
    peakEndModelKernel = maxPeriods;
end


pdfGaussian = @(X, lambda, mu1, sigma1, mu2, sigma2) lambda * normpdf(X, mu1, sigma1) + (1 - lambda) * normpdf(X, mu2, sigma2);
lowerBnd = [0, minPeriods, 0.0000001, medianPeriods, 0.0000001];
if ~logTransformed
    start = [0.7, medianPeriods, stdPeriods, max([medianPeriods, mean([peakEndModelKernel, medianPeriods])]), stdPeriods];
else
    start = [0.7, medianPeriods, stdPeriods, max(medianPeriods, peakEndModelKernel), stdPeriods];
end
upperBnd = [1, maxPeriods, Inf, maxPeriods, Inf];
[phat, ~] = mle(periods, 'pdf', pdfGaussian, 'start', start, 'LowerBound', lowerBnd, 'UpperBound', upperBnd);


if phat(2) <= phat(4)
    lambdaModelGaussian = phat(1);
    mu1ModelGaussian = phat(2);
    sigma1ModelGaussian = phat(3);
    mu2ModelGaussian = phat(4);
    sigma2ModelGaussian = phat(5); 
else
    lambdaModelGaussian = 1 - phat(1);
    mu1ModelGaussian = phat(4);
    sigma1ModelGaussian = phat(5);
    mu2ModelGaussian = phat(2);
    sigma2ModelGaussian = phat(3);
end

pdfModelGaussian = @(X) pdfGaussian(X, lambdaModelGaussian, mu1ModelGaussian, sigma1ModelGaussian, mu2ModelGaussian, sigma2ModelGaussian);
pdfPlotModelGaussian = harea * pdfModelGaussian(x);

comp1ModelGaussian = @(X) lambdaModelGaussian * normpdf(X, mu1ModelGaussian, sigma1ModelGaussian);
comp1PlotModelGaussian = harea * comp1ModelGaussian(x);

comp2ModelGaussian = @(X) (1 - lambdaModelGaussian) * normpdf(X, mu2ModelGaussian, sigma2ModelGaussian);
comp2PlotModelGaussian = harea * comp2ModelGaussian(x);

diffCompModelGaussian = @(X) comp1ModelGaussian(X) -  comp2ModelGaussian(X);

if diffCompModelGaussian(mu1ModelGaussian) * diffCompModelGaussian(mu2ModelGaussian) < 0
    if ~logTransformed
        thIntersection = fzero(diffCompModelGaussian, [mu1ModelGaussian, mu2ModelGaussian]);
    else
        thIntersection = exp(fzero(diffCompModelGaussian, [mu1ModelGaussian, mu2ModelGaussian]));
    end
elseif diffCompModelGaussian(mu1ModelGaussian) * diffCompModelGaussian(mu2ModelGaussian + 2*sigma2ModelGaussian) < 0
    if ~logTransformed
        thIntersection = fzero(diffCompModelGaussian, [mu1ModelGaussian, mu2ModelGaussian + 2*sigma2ModelGaussian]);
    else
        thIntersection = exp(fzero(diffCompModelGaussian, [mu1ModelGaussian, mu2ModelGaussian + 2*sigma2ModelGaussian]));
    end
elseif diffCompModelGaussian(mu1ModelGaussian) * diffCompModelGaussian(maxPeriods) < 0
    if ~logTransformed
        thIntersection = fzero(diffCompModelGaussian, [mu1ModelGaussian, maxPeriods]);
    else
        thIntersection = exp(fzero(diffCompModelGaussian, [mu1ModelGaussian, maxPeriods]));
    end
else
    thIntersection = 0;
end

if ~logTransformed
    thMinimum = fminbnd(pdfModelGaussian, mu1ModelGaussian, mu2ModelGaussian);
else
    thMinimum = exp(fminbnd(pdfModelGaussian, mu1ModelGaussian, mu2ModelGaussian));
end



if plotTh
    plotGaussian = figure('Visible', 'off');
    hold on;
    histogram(periods, edges, 'FaceColor', 'k')
    xlim([minPeriods, maxPeriods])
    
    plot(x, pdfPlotModelGaussian, 'c', 'DisplayName', 'double Gaussian');
    plot(x, comp1PlotModelGaussian, 'm', 'DisplayName', 'Component #1 of double Gaussian');
    plot(x, comp2PlotModelGaussian, 'r', 'DisplayName', 'Component #2 of double Gaussian');
    
    yLimits = get(gca, 'YLim');
    
    line(mu1ModelGaussian * ones(1,2), yLimits, 'Color', 'm', 'LineStyle', '-', ...
        'LineWidth', 0.5, ...
        'DisplayName', ['Mean of Component #1 = ', num2str(mu1ModelGaussian, 3)]);
    line(mu2ModelGaussian * ones(1,2), yLimits, 'Color', 'r', 'LineStyle', '-', ...
        'LineWidth', 0.5, ...
        'DisplayName', ['Mean of Component #2 = ', num2str(mu2ModelGaussian, 3)]);
    line(thIntersection * ones(1,2), yLimits, 'Color', 'c', 'LineStyle', '--', ...
        'LineWidth', 0.5, ...
        'DisplayName', ['Threshold Intersection = ', num2str(thIntersection, 3)]);
    line(thMinimum * ones(1,2), yLimits, 'Color', 'b', 'LineStyle', '--', ...
        'LineWidth', 0.5, ...
        'DisplayName', ['Threshold Minimum = ', num2str(thMinimum, 3)]);
    legend('location', 'northeast')
    xlabel('Periods (s)')
    ylabel('Count')
else
    plotGaussian = NaN;
end

end