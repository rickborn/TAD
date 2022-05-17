% fun_with_norm.m
%
% fun with the normal distribution

% Plot the standard normal distribution:
zVals = [-4:0.1:4];
pVals = normpdf(zVals,0,1);
plot(zVals,pVals);
xlabel('# of standard deviations');
ylabel('normal probability density');
hold on

% draw some relevant lines
x = [-2:2];
y = normpdf(x,0,1);
for k = 1: length(x)
    line([x(k),x(k)],[0,y(k)],'Color','b','LineStyle','--');
end
   
% fill in the region from 0 to 1 s.d.
x = [0:0.1:1];
y = normpdf(x,0,1);
area(x,y);
% How large is this area?
ar = normcdf(1,0,1) - normcdf(0,0,1);
text(0.1, max(pVals)/2, num2str(ar,3));

% fill in the region from 1 to 2 s.d.
x = [1:0.1:2];
y = normpdf(x,0,1);
area(x,y);
% How large is this area?
ar = normcdf(2,0,1) - normcdf(1,0,1);
text(1.1, max(pVals)/6, num2str(ar,3));

% fill in the region from 2 to 3 s.d.
x = [2:0.1:3];
y = normpdf(x,0,1);
area(x,y);
% How large is this area?
ar = normcdf(3,0,1) - normcdf(2,0,1);
text(2.25, max(pVals)/10, num2str(ar,3));