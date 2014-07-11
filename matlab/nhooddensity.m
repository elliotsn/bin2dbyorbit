%
% Function to return a neighbourhood fraction of non-zero grid cells in a
% 2D array. Neighbourhood must be odd.
%
% Elliot Sefton-Nash 24/10/2013
%
function density = nhooddensity(in, neighbourhood)

    [sx, sy]=size(in);
    
    density=zeros(sx, sy);
    
    % Only loop between the ULLR bounds of non-NaN cells.
    % p and q are vectors containing x and y indices of all non-NaN elements.
    [p,q] = ind2sub(size(in),find(~isnan(in)));
    
    for ix=min(p):max(p)
        for iy=min(q):max(q)
            
            % Nearest neighbour density for each pixel
            buffer = floor(neighbourhood/2);
            
            % Boundary conditions
            tmp=ix-buffer;
            if (tmp < 1)
                xmin=1;
            else
                xmin=tmp;
            end
            tmp=ix+buffer;
            if (tmp > sx)
                xmax=sx;
            else
                xmax=tmp;
            end
            
            tmp=iy-buffer;
            if (tmp < 1)
                ymin=1;
            else
                ymin=tmp;
            end
            tmp=iy+buffer;
            if (tmp > sy)
                ymax=sy;
            else
                ymax=tmp;
            end
                       
            density(ix,iy)=numel(find( ~isnan( in(xmin:xmax,ymin:ymax)) ));
            
        end
    end
    
    density = density./(neighbourhood^2);
    
end