%
% Function to bin a list of cartesian data points and return a 2D array.
%
% This function was written to process point-based (rather than 
% image-based) spacecraft data obtained by an instrument operating in
% pushbroom mode. The purpose is to fill in data gaps on gridded data 
% products where point density is low, but to first create contiguous orbit
% tracks on the target grid, before stacking those.
%
% Binning is performed intially within orbit tracks. For output grids with 
% the suffix 'Interp', Delauney triangulation is used to fill bins between 
% points in the interior of each orbit track. This minimizes the number of 
% empty bins the output grids. Grids entitles 'cnt', and 'err' represent
% the raw, non-interpolated data.
%
% A point density map on the target grid is calculated for points in each 
% orbit using the parameter, 'neighbourhood', a positive odd integer
% describing the length of the square side that defines a neighbourhood on
% the grid. Typically, 3 <= neighbourhood <= 7 for best results.
%
% Where there are empty points on the target grid and the density exceeds
% the density threshold (passed as densityThresh) grid cells are
% interpolated. 'neighbourhood' and 'densityThresh' should therefore be 
% tuned according to the expected point density and the resolution of the 
% target grid.
% 
% Grids produced:
%
%     Interpolated average        (avgInterp)
%     number of points per bin    (cnt)
%     standard deviation per bin  (err)
%     minimum value               (zminInterp)
%     maximum value               (zmaxInterp)
%
% Dependancies:
%       bin2d.m
%       nhoodDensity.m
%       interpWithDensity.m
%   
% Elliot Sefton-Nash 2013/12/19
%
function [avgInterp,cnt,err,zminInterp,zmaxInterp] = ...
    bin2dbyorbit(x,y,z,xedges,yedges,orbit, neighbourhood, densityThresh)

% Sanity check
n=numel(x);
if (n ~= numel(y)) || (n ~= numel(z)) || (numel(y) ~= numel(z))
    error('x, y and z vectors must be the same length');
end

% Setup arrays. Large grids could consume all the memory here, be careful
% on the grid size. Future update could use sparse arrays.
avgInterp=zeros(numel(yedges)-1,numel(xedges)-1);
totInterp=zeros(numel(yedges)-1,numel(xedges)-1);
nInterp=zeros(numel(yedges)-1,numel(xedges)-1);

avg=zeros(numel(yedges)-1,numel(xedges)-1);
tot=zeros(numel(yedges)-1,numel(xedges)-1);
cnt=zeros(numel(yedges)-1,numel(xedges)-1);
err=zeros(numel(yedges)-1,numel(xedges)-1);
sum_ximxbar=zeros(numel(yedges)-1,numel(xedges)-1);

orbavgInterp=zeros(numel(yedges)-1,numel(xedges)-1);
orbcnt=zeros(numel(yedges)-1,numel(xedges)-1);

% Arrays of largest and smallest storable numbers to compare to.
zminInterp=repmat(realmax,[numel(yedges)-1,numel(xedges)-1]);
zmaxInterp=repmat(realmin,[numel(yedges)-1,numel(xedges)-1]);

% Find the unique orbit numbers
orblist = unique(orbit);

for iorb = 1:numel(orblist)
    
    % Find all the grid cells involved this orbit.
    orbinds = find(orbit == orblist(iorb));
    
    % In a temporary grid, bin one orbit track at a time
    [orbavg,orbcnt,~] = bin2d(x(orbinds),y(orbinds),z(orbinds),xedges,yedges);
    
    oldinds = find(~isnan(orbcnt));
    
    % Fill in this orbit based on a density map
    density = nhood_density(orbavg, neighbourhood);
    orbavgInterp = interpWithDensity(orbavg, density, densityThresh);
        
     % For the non-interpolated grids
    cnt(oldinds) = cnt(oldinds) + orbcnt(oldinds);
    tot(oldinds) = tot(oldinds) + orbavg(oldinds);
    avg(oldinds) = tot(oldinds)./cnt(oldinds);
    sum_ximxbar(oldinds) = sum_ximxbar(oldinds) + (orbavg(oldinds) - avg(oldinds)).^2;
    
    % Find the places in this grid that contain new data.
    newinds = find(~isnan(orbavgInterp));
    
    % For the interpolated grids
    nInterp(newinds) = nInterp(newinds) + 1;
    totInterp(newinds) = totInterp(newinds) + orbavgInterp(newinds);
    avgInterp(newinds) = totInterp(newinds)./nInterp(newinds);
    
    zmininds = find(orbavgInterp(newinds) < zminInterp(newinds));
    zminInterp(newinds(zmininds)) = orbavgInterp(newinds(zmininds));
    
    zmaxinds = find(orbavgInterp(newinds) > zmaxInterp(newinds));
    zmaxInterp(newinds(zmaxinds)) = orbavgInterp(newinds(zmaxinds));

end
err=sqrt(sum_ximxbar./cnt);
zminInterp(zminInterp == realmax) = NaN;
zmaxInterp(zmaxInterp == realmin) = NaN;
