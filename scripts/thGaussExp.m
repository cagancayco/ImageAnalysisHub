function [thMinimum, thIntersection] = thGaussExp(periods)

medianPeriods = median(periods);
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
start = [0.7, medianPeriods, stdPeriods, max([medianPeriods, mean([peakEndModelKernel, medianPeriods])]), stdPeriods];
upperBnd = [1, maxPeriods, Inf, maxPeriods, Inf];
[phat, ~] = mle(periods, 'pdf', pdfGaussian, 'start', start, 'LowerBound', lowerBnd, 'UpperBound', upperBnd);


if phat(2) <= phat(4)
    sigma2ModelGaussian = phat(5); 
else
    sigma2ModelGaussian = phat(3);
end

mypdf = @(X, lambda, mu1, sigma, mu2) lambda * normpdf(X, mu1, sigma) + (1 - lambda) * exppdf(X, mu2);
lowerBnd = [0, minPeriods, 0.0000001, medianPeriods];
upperBnd = [1, maxPeriods, Inf, maxPeriods];

start = [0.7, medianPeriods, stdPeriods, max([medianPeriods, mean([peakEndModelKernel, medianPeriods])])];
[phat, ~] = mle(periods, 'pdf', mypdf, 'start', start, 'LowerBound', lowerBnd, 'UpperBound', upperBnd);

lambdaModel = phat(1);
mu1Model = phat(2);
sigmaModel = phat(3);
mu2Model = phat(4);
pdfModel = @(X) mypdf(X, lambdaModel, mu1Model, sigmaModel, mu2Model);

comp1Model = @(X) lambdaModel * normpdf(X, mu1Model, sigmaModel);
comp2Model = @(X) (1 - lambdaModel) * exppdf(X, mu2Model);
diffCompModel = @(X) comp1Model(X) - comp2Model(X);

if diffCompModel(mu1Model) * diffCompModel(mu2Model) < 0
    thIntersection = fzero(diffCompModel, [mu1Model, mu2Model]);
elseif diffCompModel(mu1Model) * diffCompModel(mu2Model + 2*sigma2ModelGaussian) < 0
    thIntersection = fzero(diffCompModel, [mu1Model, mu2Model + 2*sigma2ModelGaussian]);
elseif diffCompModel(mu1Model) * diffCompModel(maxPeriods) < 0
    thIntersection = fzero(diffCompModel, [mu1Model, maxPeriods]);
else
    thIntersection = 0;
end


if mu2Model >= mu1Model
    thMinimum = fminbnd(pdfModel, mu1Model, mu2Model);
else
    thMinimum = fminbnd(pdfModel, mu1Model, maxPeriods);
end



end