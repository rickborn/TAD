% GibbsSamplingDemo.m
%
% Simple demonstration of Gibbs sampling from a bivariate normal
% distribution (i.e. a distribution where we can sample directly from the
% join distribution, so we know the answer and can compare).
%
% Based on a video coding demo in Python by 'ritvikmath' on YouTube:
% https://www.youtube.com/watch?v=MNHIbOqH3sk
%
% RTB wrote it, early morning darkness, Oct. 29, 2025

%% First sample directly from the bivariate normal distribution

% # of samples (eventually will be passed to function)
nSamp = 10000;

% Paramters:
mu = [0,0];             % means
sigma = [1,0.5;0.5,1];  % covariance matrix

% Draw random samples from the distribution
R=mvnrnd(mu,sigma,nSamp);

% plot it
figure, subplot(2,1,1);
scatter(R(:,1),R(:,2));
axis('equal')
title('Sampling from the Joint Distribution')

%% Now we do it with Gibbs Sampling

% create a variable to hold the samples
samples = zeros(nSamp,2);
% initialize the first entry:
samples(1,1) = 1;
samples(1,2) = -1;

for k = 2:nSamp
    currY = samples(k-1,2);

    % Sample from the conditional probability distribution p(X|Y)
    % Derived here: https://www.youtube.com/watch?v=7LB1VHp4tLE
    % formula for the conditional distribution:
    % p(x|y) = Norm(y/2, 1 - rho^2)
    newX = normrnd(currY/2,sqrt(3/4));

    % Now sample from p(Y|X)
    newY = normrnd(newX/2,sqrt(3/4));

    % add these to the samples array:
    samples(k,:) = [newX,newY];
end

% Plot it:
subplot(2,1,2);
scatter(samples(:,1),samples(:,2));
axis('equal')
title('Gibbs Sampling')
% Looks similar!

%% Verify by looking at marginal distributions

% histograms of the X distributions:
figure, subplot(2,1,1)
histogram(R(:,1),25);
hold on
histogram(samples(:,1),25);
xlabel('X')
ylabel('#')
legend({'Joint','Gibbs'})

% histograms of the Y distributions:
subplot(2,1,2)
histogram(R(:,2),25);
hold on
histogram(samples(:,2),25);
xlabel('Y')
ylabel('#')
legend({'Joint','Gibbs'})

%% Show the 2 methods are different by looking at the correlations between consecutive samples

% With the random sampling from the true joint distribution, consecutive
% random samples are uncorrelated. The Gibbs sampling, by virtue of the
% fact that we are building successive samples off of the previous ones
% introduces a weak correlation.

figure, subplot(2,1,1)
scatter(R(1:end-1,1),R(2:end,1));
rJoint = corrcoef(R(1:end-1,1),R(2:end,1));
rJoint = rJoint(2,1);
axis('equal');
tstr = ['Joint Samples, \rho = ', num2str(rJoint,2)];
title(tstr);

subplot(2,1,2)
scatter(samples(1:end-1,1),samples(2:end,1));
rGibbs = corrcoef(samples(1:end-1,1),samples(2:end,1));
rGibbs = rGibbs(2,1);
axis('equal');
tstr = ['Gibbs Samples, \rho = ', num2str(rGibbs,2)];
title(tstr);