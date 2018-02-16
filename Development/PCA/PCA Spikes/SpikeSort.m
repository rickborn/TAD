% SpikeSort.m: Sorting spikes from noisy multi-unit extracellular recordings.
%
% This example was modified from: Gabbiani & Cox, Mathematics for Neuroscientists
%
% MATLAB concepts covered:
% 1. thresholding to find peaks in data
% 2. using 'diff' to find events in continuous data
% 3. PCA/SVD to reduce dimensionality of data
% 4. 3D plotting
% 5. ISI histograms

%
%  Align the active events into a data matrix, compute the associated 
%  singular value decomposition and finally plot the triplet of scores 
%  for each event. You should be able to dynamically rotate the axes of 
%  figure 3 and arrive at 3 relevant clusters (as in Fig 14.2B).
%

% Load in some data recording using an extracellular electrode.
load SpikeSortData

% This file contains 3 variables:
% 'Time' is the time-base, from 0 to 3000 ms in 0.02 ms steps.
% 'Vtotal' is the digitized voltage trace from the amplifier.
% 'Vspikes' is the digitized voltage trace for each of the 3 neurons
plot(Time,Vtotal,'k');
xlabel('time (ms)'); ylabel('mV');

% The data are very condensed, so use the zoom tool to expand parts of the
% trace. You should be able to see action potentials ('spikes') of
% three different sizes. In this case we know that there are three, because
% we generated them ourselves. We will use this knowledge below to see how
% well our spike-sorting algorithm has performed. The data from the
% individual neurons is in the variable 'Vspikes', with each column
% representing the voltage trace of the spikes generated by one neuron. You
% can prove this by plotting their sum and comparing it to Vtotal:
figure;
for k = 1:3
    subplot(3,1,k); plot(Time,Vspikes(:,k),'k');
end
figure; plot(Time,sum(Vspikes'),'k');

% What do you notice about the variability in spike heights comparing
% figures 1 and 2? Why might this be?

%% The first thing we need to do is break up the "continuous" voltage trace
% into discrete events containing single spikes. Looking at the raw voltage
% trace, we see that spike events are generally much bigger than the noise.
% We also notice that there are some very large events that are clearly
% bigger than any of the spikes we generated. What is going on here? These
% are superpositions of spikes, which, are thankfully rare in real
% recordings. How might we select the "true" spikes? As a first step, we
% might choose data points that exceed the noise:
Vlo = 1; Vhi = 4;
Vwin = find(Vtotal > Vlo)';
% Remember that Vwin contains indices to an array, not data values!

% What does this look like superimposed on our original voltage trace?
figure(1); hold on;
plot(Time(Vwin),Vtotal(Vwin),'g.');

%% So now we've selected data points that mostly belong to spikes (We'll
% still need to do some cleaning up later, removing those that exceed Vhi).
% But how can we break these up into discrete events? If you zoom in on
% this new trace, you'll see that all of the green dots within a single
% spike are neighbors, and that there are gaps between the spikes. Another
% way to see what's going on is to actually plot our indices:
figure;
plot(Vwin);

% From a distance, this looks like a smooth line, but if we zoom in, we'll
% see that it is jagged, with smooth, nearly horizontal lines broken by
% sudden vertical jumps. We can use the MATLAB 'diff' function to make
% these jumps explicit. 'Diff' simply returns the difference between
% successive members of a vector:
dVwin = diff(Vwin);
plot(dVwin);

% This is kind of handy: The values "within" a spike will now be 1s, and
% those between spikes will be very much larger. We can now use a threshold
% on this vector to find the gaps between spikes:
GapIndex = find(dVwin > 50);
Ev = Vwin(GapIndex+1);    % These are the indices to events likely to be spikes
EvNum = length(Ev);      % number of events
% Let's see if we've got the right thing:
figure(1); plot(Time(Ev),Vtotal(Ev),'r*');

%% Now we need to plow through this list of events and grab some actual data
% on both sides of the index. From looking at the raw traces, we see that
% our spikes are generally less than 10 ms in width, so let's grab 3 ms
% before and 9 ms after each event marker:
ts = linspace(0,12,601)';   % time vector for plotting single spike
% How else might we generate this vector?
% ts = 0:0.02:12;

% We're going to create a big array of spikes for the SVD analysis. Each
% row will be one spike's worth of data. We'll need a counter for the
% spikes:
SpikeCtr = 0;

% During this same 'for' loop, we'll plot each spike so we can see if the
% data are good. Since we know where the data came from, we can use our
% Vspike array to identify which neuron produced which spike and color code
% them. Neuron 1 will be red, 2 blue and 3 green:
cstrs = ['r' 'b' 'g'];      % different color for each spike

figure;
for j=1:EvNum
    % Grap a block of data
    Vrec = Vtotal(Ev(j)-150:Ev(j)+450);     % back 3 ms, ahead 9 ms
    [Vpeak,Ipeak] = max(Vrec);              % top of spike, for time
    SpIndex = Ipeak + Ev(j) - 150 + 1;      % time of occurrence of peak
    
    % This is where we'll exclude events that are too large.
    if Vpeak < Vhi
       SpikeCtr = SpikeCtr + 1;
       SpikeData(SpikeCtr,:) = Vrec;        % This holds all of the spikes
       SpikeIndices(SpikeCtr) = SpIndex;
       
       % Find out which spike we've got here.
       v = [Vspikes(Ev(j),1) Vspikes(Ev(j),2) Vspikes(Ev(j),3)];
       % The largest value identifies which spike it is. We don't care
       % about the actual value; we want the corresponding index.
       [foo,clab(SpikeCtr)] = max(v);
       
       plot(ts(1:4:end),Vrec(1:4:end),cstrs(clab(SpikeCtr))); hold on;
    end
     
end
xlabel('t  (ms)','fontsize',14); ylabel('mV','fontsize',14);

figure(1)
plot(Time(SpikeIndices),Vtotal(SpikeIndices),'m*');

% This looks pretty good. We've not got our spikes aligned in one big
% array. Now we can use singular value decomposition to help with the
% sorting. Basically, you can think of each each spike as a vector of 601
% dimensions. That is, each spike is uniquely defined by 601 numbers
% consisting of the voltage value at each of 601 time points. But perhaps
% we can effectively describe our spikes with fewer dimensions. We'll use
% SVD to find new coordinates that maximize the variability accounted for
% in our spikes.

%% Perform SVD
% The details of how SVD works are well beyond the scope of this course. 
% For now, just use it and see what it does! The basic "formula" for PCA is:
% 1. Organize a data set as an m � n matrix, where m is the number of 
% measurement types and n is the number of trials. In our case, we want to
% treat each time-point as a separate measurement type and each spike as
% a different "trial". Hence, a 601 x 98.
% 2. Subtract off the mean for each measurement type or row xi.
% 3. Calculate the SVD or the eigenvectors of the covariance.
SpikeData = SpikeData';     % now 601 x 98
SpikeNum = size(SpikeData,2);
SpikeData = SpikeData - repmat(mean(SpikeData,2),1,SpikeNum);    % subtract mean
[Y,Sig,X] = svd(SpikeData);     % All the heavy lifting done here.
sig = diag(Sig);
figure; semilogy(sig(sig>1),'kx-')      % plot the significant singular values
xlabel('index','fontsize',14); ylabel('singular value','fontsize',14)

% This is interesting. The singular value, in this case, gives us a measure
% of how much of the variance each new vector is accounting for. We can see
% that the first 3 or 4 account for a lot, and subsequent ones not so much.

%% In order to see where each spike falls in our new space, we plot each
% spike according to its projection on the first three principle components.
% To do this, we just multiply our SpikeData matrix by the vector defining
% each of the first three principle components:
% SpikeData' is 98 x 601 (98 spikes by 601 time points)
% Y(:,1) is 601 x 1
% so the matrix product is 98 x 1
% Essentially, the 601-vector defining the 1st principle component is
% multiplied by the time series for each trace to get one value (the
% projection of that trace onto the 1st principle component) for each
% spike. Spikes that have a similar projection onto the 1st component
% will be plotted near each other.
sc1 = SpikeData'*Y(:,1);
sc2 = SpikeData'*Y(:,2);
sc3 = SpikeData'*Y(:,3);
figure; plot3(sc1,sc2,sc3,'k+','markersize',10);
hold on; grid
xlabel('s_1','fontsize',14); ylabel('s_2','fontsize',14); zlabel('s_3','fontsize',14)

% Use the 'Rotate 3D' tool on the figure to see where each cell falls in
% this new space. You can see pretty clearly that the data fall into
% clusters (with a few stragglers). In most modern "Spike-sorting GUIs",
% one can use the mouse to draw circles around clusters of points in the
% PCA-space and thus define the range of values that will be accepted as a
% spike belonging to one neuron.

%% In this case, we know which spikes belong to which neurons, because we
% generated them ourselves. We can thus go back in and calculate the
% projections for each spike from each known neuron and then color code
% them.

% These are the scores for each spike generated by neuron #1
sc11 = SpikeData(:,find(clab==1))'*Y(:,1);
sc21 = SpikeData(:,find(clab==1))'*Y(:,2);
sc31 = SpikeData(:,find(clab==1))'*Y(:,3);
plot3(sc11,sc21,sc31,'ro','markersize',14);

% These are the scores for each spike generated by neuron #2
sc12 = SpikeData(:,find(clab==2))'*Y(:,1);
sc22 = SpikeData(:,find(clab==2))'*Y(:,2);
sc32 = SpikeData(:,find(clab==2))'*Y(:,3);
plot3(sc12,sc22,sc32,'bs','markersize',14)

% These are the scores for each spike generated by neuron #3
sc13 = SpikeData(:,find(clab==3))'*Y(:,1);
sc23 = SpikeData(:,find(clab==3))'*Y(:,2);
sc33 = SpikeData(:,find(clab==3))'*Y(:,3);
plot3(sc13,sc23,sc33,'gd','markersize',14)
legend('All Spikes','Spike #1','Spike #2','Spike #3','Location','NorthEast');
hold off

%% Exercise:
% One reality check for spike sorting is to look at the interspike
% intervals. How might this help?
SpikeTimes = Time(SpikeIndices);

% First, plot ISI histogram of all spikes:
figure; hist(diff(SpikeTimes));

% Next, plot ISI histogram for each spike separately
ISI = ones(length(SpikeTimes),3) .* NaN;
for k = 1:3
    isi = diff(SpikeTimes(find(clab==k)));
    ISI(1:length(isi),k) = isi;
end
figure; hist(ISI);

% Exercise 2: What would the ISI of a truly random process look like?

% Exercise 3: What do the principle components look like? Plot them:
figure;
for k = 1:3
    plot(ts(1:4:end),Y(1:4:end,k),cstrs(k));
    hold on;
end

% What does this tell you about the way the spikes were generated?

% NOTE: The original file "spikepca.m" both generated the fake spike data
% and then analyzed it with SVD. For teaching purposes, I am going to
% separate the two parts:
%   spikegen: generate then save the fake spike data as a mat-file 
%   (SpikeData.mat, with time and voltage)
%
%   spikesort: do the svd analysis stuff
%
% RTB 3/6/11