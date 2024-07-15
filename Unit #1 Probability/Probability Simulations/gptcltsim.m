% gptcltsim.m
%
% On 14 June 2023, I asked Chat-GPT4 the following:
% Write a MATLAB script that uses a simulation to illustrate the core idea
% of the Central Limit Theorem.
%
% It replied:
% The Central Limit Theorem (CLT) is a fundamental concept in statistics.
% It states that if you have a population with a mean μ and standard
% deviation σ and take sufficiently large random samples from the
% population with replacement, then the distribution of the sample means
% will be approximately normally distributed. This will hold true
% regardless of the shape of the population distribution.
%
% Here is a simple MATLAB script that illustrates the central limit
% theorem. It takes random samples of increasing size from a uniform
% distribution and calculates their means. The histograms of these means
% should become increasingly close to a normal distribution.

% Define parameters
num_samples = 1000; % number of samples to draw for each sample size
sample_sizes = [1 5 10 30 50]; % different sample sizes to use

% Create figure
figure

% Loop over different sample sizes
for i = 1:length(sample_sizes)
    
    sample_size = sample_sizes(i);
    sample_means = zeros(1, num_samples);
    
    % Draw samples and calculate their means
    for j = 1:num_samples
        sample = rand(sample_size, 1); % draw sample from uniform distribution
        sample_means(j) = mean(sample); % calculate mean of sample
    end
    
    % Create subplot for this sample size
    subplot(length(sample_sizes), 1, i)
    
    % Plot histogram of sample means
    hist(sample_means, 30) % you can adjust the number of bins here
    xlim([0 1]) % limit x-axis to range of possible means
    
    % Label subplot
    title(['Sample size = ', num2str(sample_size)])
    xlabel('Sample mean')
    ylabel('Frequency')
    
end

% Chat-GPT4 explanation:
% This script will generate a separate histogram for each sample size. As
% the sample size increases, the distribution of the sample means becomes
% increasingly bell-shaped, illustrating the central limit theorem. The
% script uses the uniform distribution for simplicity, but you could
% replace this with any other distribution to show that the CLT applies
% regardless of the shape of the original distribution.