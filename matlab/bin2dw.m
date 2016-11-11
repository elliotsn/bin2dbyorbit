%
% Function to bin 2D weighted data onto a 2D grid. Operates by looping over
% every grid cell. Most efficient when the number of bins is greater than 
% the number of points needed to be binned. Incompatible with updating 
% formulae, i.e. for datasets that cannot be entirely stored in volatile 
% memory.
%
% Inputs:
%  xedges, yedges - Vectors of bin edges. Must be monotonic ascending. 
%                   There are numel(xedges)-1 bins in the x-direction.
%  x,y            - Coordinates of each data point.
%  t              - Data point.
%  w              - Weight of each data point.
%
% Outputs:
%   avg_w  - weight mean of points in each bin
%   cnt    - number of points in each bin, irrespective of weight
%   err_w  - standard deviation on the weighted mean
%   sum_w  - total weight for all points in each bin
%   sum_wt - total weight * value for all points in each bin, this is
%            returned because it may be needed in the calculation of a 
%            rolling weighted average.
%   tmin   - minimum value of t in each bin
%   tmax   - maximum value of t in each bin
%
function [avg_w, cnt, err_w, sum_w, sum_wt, tmin, tmax] = ...
    bin2dw(x, y, t, w, xedges, yedges)

n=numel(t);

if (n ~= numel(y)) || (n ~= numel(t)) || (numel(y) ~= numel(t))
    error('x, y and t vectors must be the same length');
end

nx = numel(xedges)-1;
ny = numel(yedges)-1;

avg_w  = nan(ny,nx);
cnt    = nan(ny,nx);
err_w  = nan(ny,nx);
sum_w  = nan(ny,nx);
sum_wt = nan(ny,nx);
tmin   = nan(ny,nx);
tmax   = nan(ny,nx);

% Work in 1D coordinates for bins using ind2sub. Easiest to create point 
% mask, also more future-proof for ND binning.

% Calculate xindex and yindex for every point.
xind = getIndLinSpace(xedges, x);
yind = getIndLinSpace(yedges, y);

% Save some memory
clear x y

% Determine linear index of every point on grid.
i = sub2ind([ny nx], yind, xind);

% Discard those points that are off the grid to speed things up.
onGridMask = ~isnan(i);
i = i(onGridMask);
t = t(onGridMask);
w = w(onGridMask);

% Loop over each grid cell.
for ib = 1:(nx*ny)
    
    % Mask of points in this grid cell.
    mask = i == ib;
    
    if any(mask)
    
        % Total points
        cnt(ib) = numel(w(mask));
        
        % Sum of weight
        sum_w(ib) = sum(w(mask));

        % Sum of weight * t
        sum_wt(ib) = sum(w(mask).*t(mask));
        
        % Weighted mean
        avg_w(ib) = sum_wt(ib)./sum_w(ib);

        % Weighted standard deviation
        err_w(ib) = sqrt( sum( w(mask).*(t(mask) - avg_w(ib) ).^2) / sum_w(ib));
        
        % Extremes for this bin
        tmin(ib) = min(t(mask));
        tmax(ib) = max(t(mask));
    end 
end