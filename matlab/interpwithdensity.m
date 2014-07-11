%
% Function to perform bilinear interpolation to fill in grid cells with a
% value of zero, only where the corresponding cells in density have a value
% greater than density_thresh.
%
% density and in must be 2D arrays of the same dimensions.
%
% Elliot Sefton-Nash 20131215
%
function out = interpwithdensity(in,density,density_thresh)

% NaN any values that are below density_thresh
in( density < density_thresh ) = NaN;

% Make a tri-scattered interpolant using the non-NaN elements in in.
indsref = find(~isnan(in));
[x,y] = ind2sub(size(in),indsref);
F = TriScatteredInterp(x,y,in(indsref));

% Evaulate the interpolant at the empty locations.

% Find where there are cells == NaN in in, where density >= density_thresh
indsfill = find( (isnan(in)) & (density >= density_thresh));

[mx,my] = ind2sub(size(in), indsfill);

out = in;
out(indsfill) = F(mx,my);
