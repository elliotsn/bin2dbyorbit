%
% Function to bin a list of values of a variable t at cartesian data points
% (x, y) and return 2D arrays.
%
% This function was written to process point-based (rather than 
% image-based) spacecraft data obtained by an instrument operating in
% pushbroom mode. The purpose is to calculate various quantities for a vector of
% data points, t, whose  location is defined by 2D coordinates, x,y, and weight,
% w. 
% 
% The weighted mean and error on the weighted mean is first calculated for data 
% within each orbit track. At this point the routine offers the possibility to
% 1) Remove disconnected outlying data points that have a low weight from each 
%    orbit track.
% 2) Interpolate only for areas where data density is sufficient, within each 
%    orbit track.
% After these steps, orbit tracks are combined to produce a mapped data product.

% For output grids with the suffix 'Interp', Delauney triangulation is used to 
% fill bins between points in the interior of each orbit track. This minimizes 
% the number of empty bins the output grids. Grids without 'Interp' in the name
% represent non-interpolated data.
%
% To perform this quality control a data density map is calculated. The 
% parameter, 'neighbourhood', is a positive odd integer describing the length of
% the square side that defines a neighbourhood on the grid.
% Typically, 3 <= neighbourhood <= 7 for best results.
%
% Where there are empty points on the target grid and the observation 
% density exceeds the density threshold (passed as densityThresh) bins are
% interpolated. 'neighbourhood' and 'densityThresh' should therefore be 
% tuned according to the expected point density and the resolution of the 
% target grid.
% 
% Grids produced:
%
%     Interpolated average        (tmeanInterp)
%     number of points per bin    (tcnt)
%     sum of weight per bin       (tsumw)
%     standard deviation per bin  (terr)
%     minimum value               (tminInterp)
%     maximum value               (tmaxInterp)
%
% Dependancies:
%       bin2dw.m
%       nhoodDensity.m
%       interpWithDensity.m
%   
% Elliot Sefton-Nash 2016/09/02
%
function [tmeanInterp, tcnt, sumw, terr, tmin, tmax] = ...
    bin2dwbyorbit(x,y,t,w,xedges,yedges,orbit,neighbourhood,densityThresh)

% Sanity check
n=numel(x);
if (n ~= numel(y)) || (n ~= numel(t)) || (numel(y) ~= numel(t))
    error('x, y and t vectors must be the same length');
end

% Setup arrays. Large grids could consume all the memory here, be careful
% on the grid size. Future update could use sparse arrays.
nx = numel(xedges)-1;
ny = numel(yedges)-1;

tmeanInterp=zeros(ny,nx);
tcnt=zeros(ny,nx);
sumw=zeros(ny,nx);
sumwt=zeros(ny,nx);
terr=zeros(ny,nx);
sumtInterp=zeros(ny,nx);

% Updating statistics for calculation of weighted error

% Find the unique orbit numbers
orblist = unique(orbit);

%% First pass:
%  Data from each orbit is binned individually.
%  Rolling totals for calculation of statistics in the 2nd pass.
%  Calculate weighted mean.

for iorb = 1:numel(orblist)
    
    % Find all the grid cells involved this orbit.
    orbinds = find(orbit == orblist(iorb));
    
    % In a temporary grid, bin one orbit track at a time
    [orbtmeanw, orbtcnt, ~, orbsumw, orbsumwt, ~, ~] = ...
        bin2dw(x(orbinds),y(orbinds),t(orbinds),w(orbinds),xedges,yedges);
    
    % Indices of bins in the grid affected by this orbit.
    %inds = find(~isnan(orbtcnt));
    
    % Fill in this orbit based on a density map of the observations per bin
    density = nhooddensity(orbsumw, neighbourhood);
    
    % Density is the mean number of observations per bin in the
    % neighbourhood of each pixel.
    [orbtmeanInterp, maskOrbInterp] = interpWithDensity(orbtmeanw, density, densityThresh);

    % Here we must combine the interpolated orbit grid, counts and weight
    % onto the master grid. Actions depends on whether or not each pixel 
    % in the orbit grid contains interpolated or real data.
    
    maskOrbReal = ~isnan(orbtmeanw); % Where there is real data in this orbit
    
    % CONDENSED CONDITIONAL STRUCTURE
    %
    % If the orbit grid contains real data
    %   Reset counters to remove interpolated data - faster than checking if interpolated data exists. Real data always takes precedent.
    %   Add to the counts
    %   Add to the total weight
    %   Update the statistics for the weighted mean
    %
    % If the orbit grid contains interpolated data
    %   If the existing grid DOES NOT contain real data
    %       Add to interpolated data sum
    %       Add to interpolated data count
    
    % The way to tell if there is interpolated data or not in each bin is 
    % if the weight is NaN.
    maskReal = tsumw > 0;
    maskInterp = isnan(tsumw);
    
    % Places where real data from this orbit now trumps interpolated data.
    % Reset sumw from NaN to 0 ready to be added to
    maskNewReal = maskInterp & maskOrbReal;
    tsumw(maskNewReal) = 0;
    
    % Add to sum number of data points per bin
    tcnt(maskOrbReal) = tcnt(maskOrbReal) + orbtcnt(maskOrbReal);
    
    % Add to sum of weight per bin
    sumw(maskOrbReal) = sumw(maskOrbReal) + orbsumw(maskOrbReal);
    
    % Add to sum of weight * observation per bin
    sumwt(maskOrbReal) = sumwt(maskOrbReal) + orbsumwt(maskOrbReal);
        
    % RULE
    %  A bin with real data can never be converted to a bin with
    %  interpolated data.
    
    % Bins with interpolated data where the master grid does not contain
    % original data.
    maskNewInterp = maskOrbInterp & ~maskReal;
    % Set sumw to NaN to indicate it now contains interpolated data. It may
    % already.
    sumw(maskNewInterp) = NaN;
    % For interpolated bins, a regular mean is calculated, so we need to 
    % update the total (sumtInterp) and N (tcnt).
    sumtInterp(maskNewInterp) = sumtInterp(maskNewInterp) + orbtmeanInterp(maskNewInterp);
    % Note that the average is the average of the interpolated bins, and 
    % not the average of the observations that contribute to the bin 
    % - because there is one interpolated bin per orbit track.
    tcnt(maskNewInterp)=tcnt(maskNewInterp) + 1;
   
    % FULL CONDITIONAL STRUCTURE
    %
    % If the existing bin contains real data
    %
    %   If the orbit grid contains real data
    %       Add to the counts
    %       Add to the total weight
    %       Update the statistics for the weighted mean
    
    %   If the orbit grid contains interpolated data
    %       Discard the interpolated data
    % 
    % If the existing grid contains interpolated data
    %   
    %   If the orbit grid contains real data
    %       Reset counters to remove interpolated data
    %       Add to the counts
    %       Add to the total weight
    %       Update the statistics for the weighted mean
    %
    %   If the orbit grid contains interpolated data
    %       Add to interpolated data sum
    %       Add to interpolated data count
    %
    % If the existing grid contains no data
    %   If the orbit grid contains real data
    %       Add to the counts
    %       Add to the total weight
    %       Update the statistics for the weighted mean
    %
    %   If the orbit grid contains interpolated data
    %       Add to interpolated data sum
    %       Add to interpolated data count
    
    % Where there is a new bin with an interpolated point the
    % interpolated point is added to the total.
    
    % Where there is a bin with a real data point already in
    % place, the interpolated bin is disregarded. Real data always have
    % priority.

end

% At this point, we can calculate the weighted mean in the places where
% there is real data.
maskRealData = sumw > 0;
tmeanInterp(maskRealData) = sumwt(maskRealData)./sumw(maskRealData);

% And the regular mean, and error in places where is interpolated data.
maskInterpData = isnan(sumw);
tmeanInterp(maskInterpData) = sumtInterp(maskInterpData)./tcnt(maskInterpData);

%% Second pass: Calculate error on weighted mean. This is performed by a 
% straightforward binning of the entire point cloud on the grid, not by 
% orbit - because interpolation is already done.
% For each bin in the whole grid where there is data the weighted 
% mean and its error is calculated.
% Also here we retrieve tmin and tmax from the whole point cloud - these
% are not calculated for empty or interpolated bins.
[~, ~, err_w, ~, ~, tmin, tmax] = bin2dw(x, y, t, w, xedges, yedges);
terr(maskRealData) = err_w(maskRealData);