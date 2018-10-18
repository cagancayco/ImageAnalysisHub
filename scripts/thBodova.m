function threshold = thBodova(periods)

medianPeriods = median(periods);
stdPeriods = std(periods);
minPeriods = min(periods);
maxPeriods = max(periods);

mypdf = @(X, pSQ, pQS, mu1, sigma1, mu2, sigma2) ...
        (1 - pSQ) * normpdf(X, mu1, sigma1) +    ...
        pSQ * pQS * normpdf(X, mu1 + mu2, sqrt(sigma1^2 + sigma2^2)) + ...
        pSQ * (1 - pQS) * pQS * normpdf(X, mu1 + 2*mu2, sqrt(sigma1^2 + 2*sigma2^2)) + ...
        pSQ * (1 - pQS)^2 * pQS * normpdf(X, mu1 + 3*mu2, sqrt(sigma1^2 + 3*sigma2^2)) + ...
        pSQ * (1 - pQS)^3 * pQS * normpdf(X, mu1 + 4*mu2, sqrt(sigma1^2 + 4*sigma2^2));
    
lwrBnd = [0, 0, minPeriods, 0.0000001, 0, 0.0000001];
uprBnd = [1, 1, maxPeriods, Inf, maxPeriods, Inf];

start  = [0.5, 0.5, medianPeriods, stdPeriods, medianPeriods, stdPeriods];
[phat, ~] = mle(periods, 'pdf', mypdf, 'start', start, 'LowerBound', lwrBnd, 'UpperBound', uprBnd);

mu1Model = phat(3);
mu2Model = phat(5);


threshold = mu1Model + mu2Model;

end