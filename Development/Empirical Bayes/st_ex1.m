% st_ex1.m
%
% Example #1 from "small telescopes" paper

% To find d_33% we identify the effect size that gives a two-sample difference of
% means test, with n=30 and 33% power
n = 30; 
power = 1/3; 

% Compute effect size using an anonymous function
% 'sampsizepwer' calculates the power for a given effect size, 'x'
% 'fzero' finds the zero value of our function--note that we subtract the
% power, so it will return the value of 'x' the returns a value equal to
% 'power' (= 0.333)
% [0 2] is the range over which fzero will vary x
effect_size = fzero(@(x) sampsizepwr('t2',[0 1],x,[],n) - power, [0 2]);
disp(effect_size);

% A more intuitive way to see this is to just run through a range of
% possible effect sizes and plot the resulting power:
effectSizeRange = 0:0.01:1;
n = 30;
pwr = sampsizepwr('t2',[0 1],effectSizeRange,[],n);
hp = plot(effectSizeRange,pwr,'k.');
xlabel('Effect Size');
ylabel('Power')

hold on;

% now plot the line for 33% power
ax = axis;
h1 = line([ax(1),ax(2)],[0.33,0.33]);
set(h1,'Color','b','LineStyle','--');

% We can estimate the answer as follows:
t1 = find(pwr > 1/3);
t2 = find(pwr < 1/3);
estEffectSize = mean([effectSizeRange(min(t1)),effectSizeRange(max(t2))]);
disp(estEffectSize);