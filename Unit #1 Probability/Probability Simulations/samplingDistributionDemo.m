% samplingDistributionDemo.m

n = [5,10,50,100,500,1000];

stdOfSD = zeros(size(n));
meanOfSEM = zeros(size(n));

nSims = 1000;
for k = 1:length(n)
    % draw nSims random samples of size n(k)
    allSamp = randn(n(k),nSims);
    
    % calculate the SEM of each sample using the formula:
    allSEM = std(allSamp) ./ sqrt(n(k));
    % store the mean of all of the SEMs:
    meanOfSEM(k) = mean(allSEM);
    
    % calculate the sampling distribution for the mean
    sdMean = mean(allSamp);
    % calculate the standard deviation of the sampling distribution
    stdOfSD(k) = std(sdMean);
end

figure, plot(stdOfSD,meanOfSEM,'k*');
hold on
ax = axis;
line([ax(1),ax(2)],[ax(1),ax(2)]);