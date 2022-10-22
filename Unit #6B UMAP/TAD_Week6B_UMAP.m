% TAD_Week6B_UMAP.m

% TAD Week 6: Dimensionality Reduction To include in the mini-lecture
% 
% geometric description of multi-d data intuitive explanation of PCA
% mathematical equation for PCA reconstruction Today is a bit of a change
% from the previous weeks. Instead of focusing on statistical tests or
% models, we will focus instead on another important technique in data
% analysis: dimensionality reduction.
% 
% As you saw in today's mini-lecture, principal component analysis (PCA) is
% a powerful tool for visualizing your data in a low-dimensional space for
% two reasons: it is guaranteed to capture the axes of maximum variance in
% your data (i.e., it is optimal and objective), and it is a linear
% transform (i.e. it is easy to compute and understand). UMAP is also a
% powerful tool, but for very different reasons: it attempts to find some
% low-dimensional sub-space along which your data are scattered, and to
% reduce that sub-space to a few dimensions while preserving the
% relationship of the original data. UMAP is neither optimal nor objective,
% opting instead for flexibility -- as we will see, this has its pros and
% cons.
% 
% To illustrate these basic properties of PCA and UMAP, we start with
% MNIST, a classic machine learning dataset. Then we do a deep dive into
% UMAP with some scRNA-seq data, looking at its dependence on key
% parameters, some of its failure modes, and whether or not it really does
% what it claims to do!

% Jonah Pearl wrote it for Python

%% Loading the data 

% We will start with the MNIST dataset for the first few
% exercises. This is a set of 70,000 low-res images of handwritten digits.
% Each image is 28x28, and each image is unwrapped in the data, so we'll be
% working with a 70,000 x 28^2 (784) matrix. MNIST is a classic machine
% learning dataset -- many ML researchers spend years trying to shave off
% tenths of percentage points of accuracy on it! See eg:
% https://benchmarks.ai/mnist. It is generally considered to be an easy
% problem in machine learning; very simple classifiers perform >90%
% accurate.
% 
% Here, we won't bother with classification, but the semantic meaning of
% the data points as numbers will help us to understand what PCA is doing
% with the data.
load MNIST_Data

%% Make some sample plots

% NOTE: We do the 'k+1' thing to plot the same digits as JP
% They are in different order because of the stupid subplot thing.
plot_MNIST_sample(X,y);

%% PCA Visualization 

% This PCA tutorial is partially borrowed / entirely inspired by NeuroMatch
% Academy's Dimensionality Reduction tutorial:
% https://compneuro.neuromatch.io/tutorials/W1D4_DimensionalityReduction/...
%       student/W1D4_Tutorial3.html
% 
% Let's run PCA on the MNIST dataset and see what we get. Run the following
% cells and then answer questions 1-3.

% Normalize each variable (ie columnn) (ie pixel position) by subtracting
% its mean across the dataset
X = double(X);
X_norm = X - mean(X);

% Perform PCA on the normalized data
% In MATLAB-speak: [coeff,score,latent] = pca(X_norm);
% coef: the coefficients for the principal components; 1st PC is column 1, etc.
% score: the representations of X in the principal component space
% latent: the eigenvalues of the covariance matrix: variance explained by each PC
%
% We will give these same variables more ML-friendly names:
[components,X_transform,explained_variance] = pca(X_norm);

%% Plot the data in the new coordinates
x_pc = 1;
y_pc = 2;

figure('Name','PCA Magic');
scatter(X_transform(:,x_pc),X_transform(:,y_pc),2,uint8(y)-1);
xlabel('PC1')
ylabel('PC2')
title('MNIST in the reduced PCA space')
hc = colorbar('Ticks',[0:9]);
colormap tab10(10)
caxis([-0.5,9.5]);

%% Plot the first few principal components

figure('Name','MNIST Principal Components');
for k = 1:4
    subplot(2,2,k);
    imshow(reshape(components(:,k),28,28));
    title(['MNIST PC' num2str(k)]);
    hc = colorbar('Ticks',[-.15, -.1, -.05, 0, .05, .1, .15]);
    colormap msl
    caxis([-0.15,0.15]);
end

%% Plot the explained variance vs the i-th PC

figure('Name','Explained variance');

% all PCs
subplot(2,1,1);
hp = plot(explained_variance);
set(hp, 'LineWidth', 2);
xlabel('PC number')
ylabel('Variance explained')
title('Scree plot for MNIST PCA (variance explained by i-th component)')
ax = axis;
axis([-10, ax(2), -10000, ax(4)]);

% zoom in to 1st 20 PCs
subplot(2,1,2);
hp = plot(explained_variance(1:20,1));
set(hp, 'LineWidth', 2);
xlabel('PC number')
ylabel('Variance explained')
title('Scree plot for MNIST PCA (variance explained by i-th component)')
ax = axis;
axis([0, ax(2), -10000, ax(4)]);

%% TODO (Q2): What percentage of the total variance does the first PC explain

percent_explained_PC1 = (explained_variance(1) / sum(explained_variance)) * 100;
disp(percent_explained_PC1);

%% TODO: Reconstruct using PCs

% Using the number of PCs from question 5, reconstruct the orginal
% MNIST dataset from the transformed data and the principal components.
% Then determine the mean squared error for your reconstruction.
% 
% Hint: you'll need to use matrix multiplication, and don't forget to add
% back the column means!

n_pcs_to_use = 20;
t_components = components';
X_reconstructed = (X_transform(:,1:n_pcs_to_use) * t_components(1:n_pcs_to_use,:)) + mean(X);
reconstruction_pca_mse = mean(mean((X - X_reconstructed).^2));
disp(reconstruction_pca_mse)

plot_MNIST_sample(X_reconstructed,y)

%% Comparing PCA to random projections 

% Another common technique for performing dimensionality reduction
% (although it is becoming less common as computers get faster) is using
% random projections. Instead of calculating the optimal set of projection
% axes (i.e. the principal components), we just use random ones that happen
% to have decent statistical properties, on average.
% 
% Here, we compare one type of random projections to PCA, based on variance
% explained by the axes of the projections.



%% Helper functions

function [] = plot_MNIST_sample(Z,q)
figure('Name','Sample digits from the MNIST data set');
for k = 1:9
    subplot(3,3,k);
    imshow(reshape(Z(k+1, :),28,28),[]);
    % Need to subtract '1' because the category for '0' is '1', etc.
    title(num2str(uint8(q(k+1)) - 1));
end
end