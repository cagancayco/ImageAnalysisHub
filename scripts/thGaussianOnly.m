function [thMinimum, thIntersection] = thGaussianOnly(periods, logTransformed)

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

comp1ModelGaussian = @(X) lambdaModelGaussian * normpdf(X, mu1ModelGaussian, sigma1ModelGaussian);
comp2ModelGaussian = @(X) (1 - lambdaModelGaussian) * normpdf(X, mu2ModelGaussian, sigma2ModelGaussian);

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

end