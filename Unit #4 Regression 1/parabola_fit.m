% parabola_fit.m
%
% Fit a line to a parabola: even though there is a perfect relationship
% between x & y, the linear fit gives slope & Rsquared of 0

% Define a parabola over 0:1
x = 0:0.01:1;
y = (3.*x).*(1 - x);

% plot it
figure, plot(x,y,'b-');
title('f(x) = 3x(1-x)');
xlabel('x'); ylabel('f(x)');
hold on

% perform regression
mdl1 = fitglm(x,y,'linear','Distribution','normal','Link','identity');
b0 = mdl1.Coefficients{1,'Estimate'};
b1 = mdl1.Coefficients{2,'Estimate'};

% plot the regression line:
y_fit_mdl1 = b0 + b1.*x;
plot(x, y_fit_mdl1, 'k--');

% HOWEVER, we can add an x^2 term to our model:
X = [x', x'.^2];
mdl2 = fitglm(X,y','linear','Distribution','normal','Link','identity');
b0 = mdl2.Coefficients{1,'Estimate'};
b1 = mdl2.Coefficients{2,'Estimate'};
b2 = mdl2.Coefficients{3,'Estimate'};

y_fit_mdl2 = b0 + b1.*x + b2.*(x.^2);
plot(x, y_fit_mdl2, 'r--');

% Now we get a perfect fit!