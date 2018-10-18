function threshold = thresholdGaussianSingle(AllPeriods)

model = fitgmdist(AllPeriods, 1);
mu = model.mu;
sigma = model.Sigma;

threshold = mu + 2*sigma;
end