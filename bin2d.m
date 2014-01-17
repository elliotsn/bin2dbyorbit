%
% Function to bin xy data onto a 2D grid.
%
function [avg,cnt,err] = bin2d(x,y,z,xedges,yedges)

n=numel(x);

if (n ~= numel(y)) || (n ~= numel(z)) || (numel(y) ~= numel(z))
    error('x, y and z vectors must be the same length');
end

avg=nan(numel(yedges)-1,numel(xedges)-1);
total=nan(numel(yedges)-1,numel(xedges)-1);
cnt=nan(numel(yedges)-1,numel(xedges)-1);
err=nan(numel(yedges)-1,numel(xedges)-1);
sum_ximxbar=nan(numel(yedges)-1,numel(xedges)-1);

for i = 1:n
    
    % mask of one element describing where this point should go
    xmask=find(diff(x(i) < xedges));
    ymask=find(diff(y(i) < yedges));
    
    if ~isempty(xmask) && ~isempty(ymask)
        % To separate NaNs from zeros that represent actual data.
        if isnan(total(ymask,xmask)) % First time this cell has been addressed.
            total(ymask,xmask) = z(i);
            cnt(ymask,xmask) = 1;
            avg(ymask,xmask)=total(ymask,xmask)/cnt(ymask,xmask);
            sum_ximxbar(ymask,xmask) = sum_ximxbar(ymask,xmask) + (z(i) - avg(ymask,xmask)).^2; 
        else
            total(ymask,xmask) = total(ymask,xmask) + z(i);
            cnt(ymask,xmask) = cnt(ymask,xmask) + 1;
            avg(ymask,xmask)=total(ymask,xmask)/cnt(ymask,xmask);
            sum_ximxbar(ymask,xmask) = sum_ximxbar(ymask,xmask) + (z(i) - avg(ymask,xmask)).^2; 
        end
    end
end
err=sqrt(sum_ximxbar./cnt);
