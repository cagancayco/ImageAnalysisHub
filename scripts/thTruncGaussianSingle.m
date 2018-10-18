function threshold = thTruncGaussianSingle(periods, logTransformed)

% Determining threshold for truncating Gaussian

if ~logTransformed
    meanPeriods = mean(periods);
else
    meanPeriods = log(mean(exp(periods)));
end
stdPeriods = std(periods);
minPeriods = min(periods);
maxPeriods = max(periods);
npts = length(periods);
nBins = round(1.87 * (npts - 1) ^ 0.4);
binWidth = max(periods)/nBins;
harea = binWidth * npts;
x = linspace(minPeriods, maxPeriods, 500)';



modelKernel   = fitdist(periods, 'Kernel', 'Bandwidth', stdPeriods/5);
pdfModelKernel = @(X) pdf(modelKernel, X);
pdfPlotModelKernel = harea * pdfModelKernel(x);



[~, idxModeApprox] = max(pdfPlotModelKernel);
modeModelKernelApprox = x(idxModeApprox);
skewness = abs(modeModelKernelApprox - meanPeriods);
leftBnd = max(minPeriods, modeModelKernelApprox - skewness);
rightBnd = min(maxPeriods, modeModelKernelApprox + skewness);
mu1ModelKernel = fminbnd(@(X) -pdfModelKernel(X), leftBnd, rightBnd);



toTheRight = x >= mu1ModelKernel;
xRight = x(toTheRight);
pdfRight = pdfPlotModelKernel(toTheRight);
[pksRight, locsRight] = findpeaks(pdfRight);

nToTheRight = length(locsRight);

if nToTheRight == 0
    model = fitgmdist(periods, 1);
    mu = model.mu;
    sigma = model.Sigma;

    if ~logTransformed
        threshold = mu + 2*sigma;  
    else
        threshold = exp(mu + 2*sigma);
    end
else
    thresholds = zeros(nToTheRight, 1);
    voids = zeros(nToTheRight, 1);
    
    for i = 1:nToTheRight
        thresholds(i) = fminbnd(pdfModelKernel, mu1ModelKernel, xRight(locsRight(i)));
        voids(i) = 1 - pdfModelKernel(thresholds(i))/sqrt(pdfModelKernel(mu1ModelKernel)*pksRight(i));  
    end
    
    [~, iSelected] = max(voids);
    thresholdModelKernel = thresholds(iSelected);
    
    periodsTrunc = periods(periods < thresholdModelKernel);
    
    model = fitgmdist(periodsTrunc, 1);
    mu = model.mu;
    sigma = model.Sigma;

    if ~logTransformed
        threshold = mu + 2*sigma;  
    else
        threshold = exp(mu + 2*sigma);
    end 
    
end

end