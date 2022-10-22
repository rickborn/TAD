function Y = processLabelsMNIST(filename)
% The processLabelsMNIST function operates similarly to the
% processImagesMNIST function. After opening the file and reading the magic
% number, it reads the labels and returns a categorical array containing
% their values.

[fileID,errmsg] = fopen(filename,'r','b');

if fileID < 0
    error(errmsg);
end

magicNum = fread(fileID,1,'int32',0,'b');
if magicNum == 2049
    fprintf('\nRead MNIST label data...\n')
end

numItems = fread(fileID,1,'int32',0,'b');
fprintf('Number of labels in the dataset: %6d ...\n',numItems);

Y = fread(fileID,inf,'unsigned char');

% I may want to comment this next line out. Not sure yet.
Y = categorical(Y);

fclose(fileID);
end