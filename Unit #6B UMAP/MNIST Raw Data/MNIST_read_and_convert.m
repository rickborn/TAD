% MNIST_read_and_convert.m
%
% RTB wrote it, 22 September 2022, for use with Jonah Pearl's UMAP exercise

% Download the MNIST files from http://yann.lecun.com/exdb/mnist/ and load
% the data set into the workspace. To load the data from the files as
% MATLAB arrays, place the files in the working directory, then use the
% helper functions processImagesMNIST and processLabelsMNIST

% NOTE: I had to first use 7zip to uncompress the gz files, then use the
% names of the uncompressed files. I also had to re-write the file that
% reads in and converts the data files.

% Names of data files:
filenameImagesTrain = 'train-images.idx3-ubyte';
filenameLabelsTrain = 'train-labels.idx1-ubyte';
filenameImagesTest = 't10k-images.idx3-ubyte';
filenameLabelsTest = 't10k-labels.idx1-ubyte';

% Open, read and convert to MATLAB arrays:
XTrain = myProcessImagesMNIST(filenameImagesTrain);
YTrain = myProcessLabelsMNIST(filenameLabelsTrain);
XTest = myProcessImagesMNIST(filenameImagesTest);
YTest = myProcessLabelsMNIST(filenameLabelsTest);

% Reality check: Plot sample of digits
for k = 1:9
    subplot(3,3,k);
    imshow(reshape(XTrain(k+1, :),28,28),[]);
    title(num2str(double(YTrain(k+1)) - 1));
end

% Concatenate train & test data
X = [XTrain; XTest];
y = [YTrain; YTest];

% Save it as a single MATLAB data file:
save MNIST_Data X y
        

