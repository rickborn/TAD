% effectSizeDemo.m 

% How likely are we to get a given effect size under H0?
% 
% Inspired by a failure-to-listen of Slater Sharp
% 
% RTB wrote it, 19 Oct. 2017, gorgeous autumn day
% RTB added Type S & Type M error calculations, 19 May 2023, Swan's Island
% at the Hermit Hut, courtesy of East Point Mike
% 
% For explanation of Type S and Type M errors, see:
% Gelman, A. and Carlin, J., 2014. Beyond power calculations: Assessing
% type S (sign) and type M (magnitude) errors. Perspectives on
% Psychological Science, 9(6), pp.641-651.
% 
% see also: retrodesign.m

%% Distribution of d' under H0

% Assume H0 is true. Draw 'n' samples from same distribution and measure
% dPrime. How is it distributed? What is the probability of getting a value
% of dPrime > 0.57 by chance?

dPrimeCrit = 0.57;
nSims = 10000;
allN = [5,10,20,50];
pGreaterThanCrit = zeros(size(allN));
powerSim = zeros(size(allN));
yCrit = -log10(0.05);

rng shuffle
figure('Name','Distribution of d-prime under H0');
for k = 1:length(allN)
    subplot(2,2,k);
    allSamples = randn(allN(k),2,nSims);
    allDprime = squeeze(diff(mean(allSamples)));
    
    % histogram of d-prime values
    yyaxis left
    histogram(allDprime);
    hold on
    xlabel('d-Prime'); ylabel('#');
    title(['N=' num2str(allN(k))]);
    %set(gca,'FontSize',12);
    
    % plot of p-values for a 2-sample t-test
    yyaxis right
    [h,p] = ttest2(squeeze(allSamples(:,1,:)),squeeze(allSamples(:,2,:)));
    h = logical(h);
    plot(allDprime(~h),-log10(p(~h)),'k.');
    plot(allDprime(h),-log10(p(h)),'r.');
    ylabel('-log10(p-Value)');
        
    if k == 1
        %maxAx = axis;
        maxAx = [-2.5,2.5,0,750];
    end
    ax = axis;
    axis([maxAx(1),maxAx(2),ax(3),ax(4)]);
    line([maxAx(1),maxAx(2)],[yCrit,yCrit],'Color','r');
        
    % What is the probability of getting a d-prime greater than a given value?
    pGreaterThanCrit(k) = sum(abs(allDprime) >= dPrimeCrit) / nSims;
    %text(maxAx(1)+1,(0.9*(ax(4)-ax(3))+ax(3)),...
        %['p(d''>' num2str(dPrimeCrit) ') = ' num2str(pGreaterThanCrit(k),2)]);
        
    % What is the probability of getting a false positive?
    powerSim(k) = sum(h) / nSims;
%     text(maxAx(1)+1,(0.9*(ax(4)-ax(3))+ax(3)),...
%         ['p(FP) = ' num2str(pFP(k),2)]);
    
    % What will be the median effect size published?
    pubES = median(allDprime(h==1 & allDprime' > 0));
    pubESleft = median(allDprime(h==1 & allDprime' < 0));
%     line([pubESright,pubESright],[ax(3),ax(4)],'Color','b','LineStyle','--');
%     line([pubESleft,pubESleft],[ax(3),ax(4)],'Color','b','LineStyle','--');
    plot(pubES,yCrit,'rs','MarkerFaceColor','r');
    plot(pubESleft,yCrit,'rs','MarkerFaceColor','r');
    line([pubES,pubES],[0,yCrit],'Color','r');
    line([pubESleft,pubESleft],[0,yCrit],'Color','r');
end

%% What if there is a real effect?

% Ans. It will still be inflated in the published literature if power is
% low.

% define better colors:
red = [0.6350, 0.0780, 0.1840];
gold = [0.9290, 0.6940, 0.1250];

% Cohen suggested that d=0.2 be considered a 'small' effect size, 0.5
% represents a 'medium' effect size and 0.8 a 'large' effect size.
alpha = 0.05;   % criterion for significance
trueEffectSize = 0.5;

% Make plots?
pFlag = 1;

if pFlag
    allN = [5,10,25,50];
    s = sprintf('Distribution of d'' with effect size of %0.2f',trueEffectSize);
    figure('Name',s);
else
    allN = [5,10,20,50,100,200,500];
end

nSims = 10000;
powerSim = zeros(size(allN));
pGreaterThanCrit = powerSim;
pubES = powerSim;
TypeS = powerSim;
TypeM = powerSim;
yCrit = -log10(alpha);

rng default
for k = 1:length(allN)
    
    % Idealized simulation using standard normal distribution
    % For a different approach, using raw effect size and standard error
    % (and the corresponding t-distribution), see retrodesign.m
    allSamples = randn(allN(k),2,nSims);
    
    % To simulate a real d-prime, we just add our mu-offset to one of the
    % samples (already normalized to units of s.d.).
    allSamples(:,2,:) = allSamples(:,2,:) + trueEffectSize;
    allDprime = squeeze(diff(mean(allSamples)));    % SD = 1
    
    % Do a one-sided t-test:
    %[h,p] = ttest2(squeeze(allSamples(:,2,:)),squeeze(allSamples(:,1,:)),'Tail','Right');
    % Do a 2-sided t-test:
    [h,p] = ttest2(squeeze(allSamples(:,2,:)),squeeze(allSamples(:,1,:)),'Alpha',alpha);
    h = logical(h);
    
    % What is the power?
    powerSim(k) = sum(h) / nSims;
    
    % What will be the median effect size published?
    % Following Gelmin & Carlin (2014), we'll now add calculation of their
    % 'Type S' error, which is the probability that the published effect
    % will be of the wrong sign. To do this, we need to know the
    % probability density in the right (= correct) and left (= wrong) tail
    pubES(k) = mean(abs(allDprime(h))); % Gelmin & Carlin's 'Type M' or 'exaggeration'
    %pubES(k) = median(allDprime(h==1 & allDprime' > 0));   % my original way
    % powRight = sum(h == 1 & allDprime' > 0) / nSims;
    powLeft = sum(h == 1 & allDprime' < 0) / nSims;
    TypeS(k) = powLeft / powerSim(k);
    TypeM(k) = pubES(k) / trueEffectSize;
    
    if pFlag
        subplot(2,2,k);
        % histogram of d-prime values
        yyaxis left
        histogram(allDprime);
        hold on
        xlabel('Cohen''s d'''); ylabel('#');
%         s = sprintf('N = %d, true Effect Size (ES) = %0.2f',...
%             allN(k),trueEffectSize);
%         title(s)
        title(['N = ' num2str(allN(k),2), ', true Effect Size = ', num2str(trueEffectSize,2),...
            ', \alpha = ', num2str(alpha,3)]);
        set(gca,'FontSize',12);
                
        % plot of p-values for a 2-sample t-test
        yyaxis right
        plot(allDprime(~h),-log10(p(~h)),'k.');
        %plot(allDprime(h),-log10(p(h)),'r.');
        plot(allDprime(h),-log10(p(h)),...
            'Color',red,...
            'LineStyle','none',...
            'Marker','.');
        ylabel('-log10(p-Value)');
        
        if k == 1
            %maxAx = axis;
            maxAx = [-2.5,2.5,0,750];
        end
        ax = axis;
        axis([maxAx(1),maxAx(2),ax(3),ax(4)]);
        line([maxAx(1),maxAx(2)],[yCrit,yCrit],'Color',gold);
                
        text(maxAx(1)+0.2,(0.9*(ax(4)-ax(3))+ax(3)),...
            ['Power = ' num2str(powerSim(k),2)],'FontSize',12);
        
        text(maxAx(1)+0.2,(0.8*(ax(4)-ax(3))+ax(3)),...
            ['Mean Sig. ES = ' num2str(pubES(k),2)],'FontSize',12);
        
        text(maxAx(1)+0.2,(0.7*(ax(4)-ax(3))+ax(3)),...
            ['Type S error = ' num2str(TypeS(k),2)],'FontSize',12);
        
        text(maxAx(1)+0.2,(0.6*(ax(4)-ax(3))+ax(3)),...
            ['Type M error = ' num2str(TypeM(k),2)],'FontSize',12);
        
        %pubESleft = median(allDprime(h==1 & allDprime' < 0));
        %     line([pubESright,pubESright],[ax(3),ax(4)],'Color','b','LineStyle','--');
        %     line([pubESleft,pubESleft],[ax(3),ax(4)],'Color','b','LineStyle','--');
        plot(pubES(k),yCrit,'Color',gold,...
            'LineStyle','none',...
            'Marker','square',...
            'MarkerFaceColor',gold);
        %plot(pubESleft,yCrit,'rs','MarkerFaceColor','r');
        line([pubES(k),pubES(k)],[0,yCrit],'Color',gold);
        %line([pubESleft,pubESleft],[0,yCrit],'Color','r');
        
        % vertical line for true effect size:
        ax = axis;
        hl = line([trueEffectSize,trueEffectSize],[ax(3),ax(4)]);
        set(hl,'Color',red,'LineStyle','--','LineWidth',1.5);
    end
end

%% Make a plot of "effect size inflation" as a function of power

relBiasFlag = 1;
blue = [0, 0.4470, 0.7410];
figure('Name','Effect Size Inflation');
if relBiasFlag
    hp=plot(powerSim .* 100,((pubES - trueEffectSize) ./ trueEffectSize),...
        'Color',blue,...
        'LineStyle','-',...
        'LineWidth',2,...
        'Marker','o',...
        'MarkerFaceColor',blue);
    ylabel('Relative bias');
    ax = axis;
    line([ax(1),ax(2)],[0,0],'Color',red,'LineStyle','--');
else
    hp=semilogy(powerSim,pubES ./ trueEffectSize,'bo-');
    set(hp,'MarkerFaceColor','b','LineWidth',2);
    hold on
    ylabel('Effect size inflation');
    axis([0.1, 0.85, 0.9, 2.8]);
    ax = axis;
    line([ax(1),ax(2)],[1,1],'Color','r','LineStyle','--');
end
xlabel('Statistical power of study (%)');

% Hah! I independently discovered a cool thing. See:
% https://www.nature.com/articles/nrn3475/figures/5
% They express their y-axis as "relative bias", which, in my example would
% be: (pubESright - dPrimeSim) ./ dPrimeSim
%
% Reference:
% Button, K., Ioannidis, J., Mokrysz, C. et al. Power failure: why small
% sample size undermines the reliability of neuroscience. Nat Rev Neurosci
% 14, 365–376 (2013).
