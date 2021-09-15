% genSeminarQuestionData.m
%
% reverse engineer seminar question data to get simulated raw data
%
% RTB wrote it, Friday, Sept. 3, 2021

% NOTE: I pulled the values below from the graphic in the piece by The
% Economist, so they might not be exactly correct
%
% Each datum represents a value derived from one academic seminar. Values
% are percentage of questions from women minus the percentage of seminar
% attendees who were women. Positive values indicate that women asked more
% questions; negative values indicate that men asked more questions. The
% two variables are the values for each seminar when a woman asked the 1st
% question ('womanFirst') vs. when a man asked the 1st question
% ('manFirst').
%
% womanFirst = [60,44,36,36,28,28,28,28,28,24,20,16,16,12,8,8,8,8,8,8,8,8,...
%     4,4,4,4,4,4,0,0,0,0,0,0,0,-4,-4,-4,-4,-4,-4,-8,-8,-8,-8,-8,...
%     -12,-12,-12,-12,-12,-12,-12,-12,-12,-16,-16,-16,-16,-16,-16,-16,-16,-16,...
%     -20,-20,-20,-20,-20,-20,-20,-20,-24,-24,-28,-32,-36,-40];
% 
% manFirst = [24,24,20,20,16,16,16,16,16,16,16,12,8,4,4,4,0,0,0,0,0,...
%     -4,-4,-4,-4,-4,-4,-4,-4,-4,-8,-8,-8,-8,-8,-8,-12,-12,-12,-12,-12,...
%     -12,-12,-12,-12,-12,-12,-12,-12,-12,-12,-16,-16,-16,-16,-20,-20,...
%     -20,-20,-20,-20,-20,-24,-24,-24,-24,-24,-24,-24,-24,-24,-24,-24,-24,...
%     -24,-24,-24,-24,-24,-24,-24,-24,-24,-24,-28,-28,-28,-28,-28,-28,-28,-28,-28,...
%     -32,-32,-32,-32,-32,-32,-32,-32,-32,-32,-36,-36,-36,-36,-36,-36,-36,...
%     -36,-36,-36,-40,-40,-40,-40,-40,-40,-40,-40,-40,-40,-40,-44,-44,-44,-44,...
%     -44,-44,-44,-44,-44,-44,-44,-44,-48,-48,-48,-48,-48,-48,-48,-48,-52,-52,...
%     -52,-52,-52,-52,-52,-56,-56,-56,-56,-56,-56,-56,-56,-60,-60,-60,-60,...
%     -64,-64,-64,-64,-68,-72,-72,-76];

%% Load the processed data:
load seminarQuestionData

% a few useful numbers:
nWF = length(womanFirst);
nMF = length(manFirst);
nTotal = nWF + nMF;

%% Generate the 'raw' data: each row is a question and each column is a seminar

% NOTE: Under the assumption that all the seminars were attended by 50%
% women, the largest effect sizes we can get are +50% (women asked all
% questions) and -50% (men asked all questions). So, to account for the
% values outside of this range, we will need to vary the attendance
% percentages.

nQperSeminar = 100;
allWFQuestions = zeros(nQperSeminar,nWF);
allWFQuestions(1,:) = ones(1,nWF);

% Variable percentage of women attending
percentAttendeesWomen = randi([30,50],nWF,1);

for k = 1:nWF
    thisSeminar = zeros(nQperSeminar,1);
    percentQfromWomen = womanFirst(k) + percentAttendeesWomen(k);
    
    % These correct for the values the are >50 or <-50
    while percentQfromWomen > 100
        percentAttendeesWomen(k) = percentAttendeesWomen(k) - 10;
        percentQfromWomen = womanFirst(k) + percentAttendeesWomen(k);
    end
    
    while percentQfromWomen < 0
        percentAttendeesWomen(k) = percentAttendeesWomen(k) + 10;
        percentQfromWomen = womanFirst(k) + percentAttendeesWomen(k);
    end
    
    nQfromWomen = round((percentQfromWomen/100) * nQperSeminar);
    thisSeminar(2:nQfromWomen,1) = ones(nQfromWomen-1,1);
    % shuffle every row but the first one:
    myShuffle = randperm(nQperSeminar-1) + 1;
    allWFQuestions(myShuffle,k) = thisSeminar(2:end); 
end

%% Reality check: Do we get the right answer back out?

simWF = ((sum(allWFQuestions) ./ nQperSeminar) * 100) - percentAttendeesWomen';
