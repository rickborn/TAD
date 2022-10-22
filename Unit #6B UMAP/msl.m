function c = msl

% MSL.M: luminance-corrected color map, invented by Marge Livingstone
%
%   MSL returns a 256-by-3 matrix containing Marge's custom colormap.
%
%   For example, to reset the colormap of the current figure:
%
%             colormap(msl)
%
%   See also HSV, GRAY, HOT, BONE, COPPER, PINK, FLAG, 
%   COLORMAP, RGBPLOT.

it=jet(256); c=jet(256);

for t=1:170
    c(t,3)=it(t+32,2);
end

for l=129:202
    c(l,3)=(202-l)/73;
end

for t=44:160
    c(t,1)=it(t-26,2);
end

for l=225:256
    c(l,1)=(256-l)/96+65/96;
end

for l=1:128
    c(l,1)=l/128;
end

for l=25:100
    c(l,2)=(l-25)/75;
end

for l=144:249
    c(l,2)=(249-l)/105;
end

% RTB addition for saccade maps
c([1:3],:) = 0;