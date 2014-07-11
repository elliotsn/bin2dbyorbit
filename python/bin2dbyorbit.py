# -*- coding: utf-8 -*-
"""
Created on Thu Jun  5 14:55:23 2014

@author: elliotsefton-nash
"""

def usage():
    print """
    
    ------------
    bin2dbyorbit
    ------------
    
     Program to bin a list of cartesian data points and return 2D arrays. 
     Intended for mapping point-based data from instruments aboard spacecraft.
    
     These functions was written to process point-based (rather than 
     image-based) spacecraft data obtained by an instrument operating in
     pushbroom mode. The purpose is to fill in data gaps on gridded data 
     products where point density is low, but to first create contiguous orbit
     tracks on the target grid, before stacking those.
    
     Binning is performed intially within orbit tracks. For output grids with 
     the suffix 'Interp', Delauney triangulation is used to fill bins between 
     points in the interior of each orbit track. This minimizes the number of 
     empty bins the output grids. Grids entitles 'cnt', and 'err' represent
     the raw, non-interpolated data.
    
     A point density map on the target grid is calculated for points in each 
     orbit using the parameter, 'neighbourhood', a positive odd integer
     describing the length of the square side that defines a neighbourhood on
     the grid. Typically, 3 <= neighbourhood <= 7 for best results.
    
     Where there are empty points on the target grid and the density exceeds
     the density threshold (passed as densityThresh) grid cells are
     interpolated. 'neighbourhood' and 'densityThresh' should therefore be 
     tuned according to the expected point density and the resolution of the 
     target grid.
     
     Usage:
    
        bin2dbyorbit(x,y,z,xedges,yedges,orbit, neighbourhood, densityThresh)
    
     Grids produced:
    
         Interpolated average        (avgInterp)
         number of points per bin    (cnt)
         standard deviation per bin  (err)
         minimum value               (zminInterp)
         maximum value               (zmaxInterp)
    
    """

def bin2d(x,y,z,xedges,yedges,res):

    import numpy as np    
    
    n=x.size

    if (n != y.size) | (n != z.size) | (y.size != z.size):
        error('x, y and z vectors must be the same length')
    
    nx=xedges.size-1
    ny=yedges.size-1

    # Setup arrays. Large grids could consume all the memory here.
    # Assumes a responsible human-being is the operator.
    # Initialize as linear arrays first, then reshape at the end.
    avg=np.zeros((nx,ny)).ravel()
    tot=np.zeros((nx,ny)).ravel()
    cnt=np.zeros((nx,ny)).ravel()
    err=np.zeros((nx,ny)).ravel()
    sum_ximxbar=np.zeros((nx,ny)).ravel()

    # For every point, an x,y pair.
    xinds = np.floor((x-xmin)*1/res)
    yinds = np.floor((y-ymin)*1/res)
    
    # Convert sub to ind to find only grid cells that are occupied.
    inds = np.ravel_multi_index((xinds.astype('int'),yinds.astype('int')),(nx,ny))
    
    # For every unique linear index find all points in the grid cell.
    indsuniq = np.unique(inds)
    for ib in indsuniq:

        # All elements in the grid corresponding to this bin
        theseinds = np.where(inds == ib)[0]
        
        # Gather them into the bin and do stats
        tot[ib] += np.sum(z[theseinds])
        cnt[ib] += len(theseinds)
        avg[ib] =  tot[ib]/cnt[ib]
        sum_ximxbar[ib] += np.sum((z[theseinds] - avg[ib])**2)
        err[ib] = np.sqrt(sum_ximxbar[ib]/cnt[ib])

    avg=np.reshape(avg, (nx,ny))
    #tot=np.reshape(tot, (nx,ny))
    cnt=np.reshape(cnt, (nx,ny))
    err=np.reshape(err, (nx,ny))

    return avg,cnt,err


# Function to return a neighbourhood fraction of non-zero grid cells in a
# 2D array. Neighbourhood must be odd.
def nhood_density(grid, nhood):
    
    import numpy as np    
    
    (sx,sy)=grid.shape
    density=np.zeros((sx, sy))
    
    # Only loop between the ULLR bounds of non-zero cells.
    # p and q are vectors containing x and y indices of all non-NaN elements.
    p,q =  np.unravel_index(np.where(grid.ravel()!=0)[0], grid.shape)
    
    for ix in range(min(p),max(p)+1):
        for iy in range(min(q),max(q)+1):
            # Nearest neighbour density for each pixel
            buff = np.floor(nhood/2)

            # Boundary conditions
            tmp=ix-buff
            if (tmp < 1):
                xmin=1
            else:
                xmin=tmp

            tmp=ix+buff
            if (tmp > sx):
                xmax=sx
            else:
                xmax=tmp
                
            tmp=iy-buff
            if (tmp < 1):
                ymin=1
            else:
                ymin=tmp
                
            tmp=iy+buff
            if (tmp > sy):
                ymax=sy
            else:
                ymax=tmp

            # Number of cells in this moving window that are non-zero.
            density[ix,iy] = len(np.where(grid[xmin:xmax+1,ymin:ymax+1].ravel()!=0)[0])

    return density/(nhood**2)

# Function to perform bilinear interpolation to fill in grid cells with a
# value of zero, only where the corresponding cells in density have a value
# greater than density_thresh.
#
# density and grid must be 2D arrays of the same dimensions.
#
def interpWithDensity(grid,density,density_thresh):
    
    # NaN any values that are below density_thresh
    in( density < density_thresh ) = NaN;
    
    # Make a tri-scattered interpolant using the non-NaN elements in in.
    indsref = find(~isnan(in));
    [x,y] = ind2sub(size(in),indsref);
    F = TriScatteredInterp(x,y,in(indsref));
    
    # Evaulate the interpolant at the empty locations.
    
    # Find where there are cells == NaN in in, where density >= density_thresh
    indsfill = find( (isnan(in)) & (density >= density_thresh));
    
    [mx,my] = ind2sub(size(in), indsfill);
    
    out = in;
    out(indsfill) = F(mx,my);

    
def bin2dbyorbit(x,y,z,orbit,xedges,yedges,res,nhood,densitythresh): 
    
    import numpy as np
    
    nx=xedges.size-1
    ny=yedges.size-1

    # Setup arrays. Large grids could consume all the memory here.
    avgInterp=np.zeros((nx,ny))
    totInterp=np.zeros((nx,ny))
    nInterp=np.zeros((nx,ny))
    
    avg=np.zeros((nx,ny))
    tot=np.zeros((nx,ny))
    cnt=np.zeros((nx,ny))
    err=np.zeros((nx,ny))
    sum_ximxbar=np.zeros((nx,ny))

    orbavgInterp=np.zeros((nx,ny))
    orbcnt=np.zeros((nx,ny))

    # Arrays of largest and smallest storable numbers to compare to.
    realmax=np.finfo('d').max
    realmin=np.finfo('d').min
    zminInterp=np.full((nx,ny),realmin)
    zmaxInterp=np.full((nx,ny),realmax)

    orblist = np.unique(orbit)

    for thisorb in orblist:
        
        # Return a bool that's true for points in this orbit.
        orbinds = np.where(orbit == thisorb)[0]
        
        # In a temporary grid, bin one orbit track at a time
        orbavg,orbcnt,_ = bin2d(x[orbinds],y[orbinds],z[orbinds],xedges,yedges,res)
        
        # 1D indices of the grid cells that actually have data in them
        oldinds = np.where(orbcnt.ravel() > 0)[0]
        
        # The density of points only considers that cells are non-zero.
        density = nhood_density(orbcnt, nhood)

        # Fill in this orbit based on a density map       
        orbavgInterp = interpWithDensity(orbavg, density, densityThresh)

        # For the non-interpolated grids
#        cnt[oldinds] = cnt[oldinds] + orbcnt[oldinds]
#        tot[oldinds] = tot[oldinds] + orbavg[oldinds]
#        avg[oldinds] = tot[oldinds]/cnt[oldinds]
#        sum_ximxbar[oldinds] = sum_ximxbar[oldinds] + (orbavg[oldinds] - avg[oldinds])^2
#        
#        # Find the places in this grid that contain new data.
#        newinds = np.where(orbavgInterp == np.NaN)
#        
#        # For the interpolated grids
#        nInterp[newinds] = nInterp[newinds] + 1
#        totInterp[newinds] = totInterp[newinds] + orbavgInterp[newinds]
#        avgInterp[newinds] = totInterp[newinds]/nInterp[newinds]
#        
#        zmininds = find(orbavgInterp[newinds] < zminInterp[newinds])
#        zminInterp(newinds[zmininds]) = orbavgInterp(newinds[zmininds])
#        
#        zmaxinds = find(orbavgInterp[newinds] > zmaxInterp[newinds])
#        zmaxInterp(newinds[zmaxinds]) = orbavgInterp(newinds[zmaxinds])
    
    #err=sqrt(sum_ximxbar/cnt)
    #zminInterp(zminInterp == realmax) = np.NaN
    #zmaxInterp(zmaxInterp == realmin) = np.NaN

    return avgInterp,cnt,err,zminInterp,zmaxInterp


# Function to write out a binary column grid to a path.       
#def writegrid(fpath,ingrid):

        
#    return success

    
# Main function - handle cmdline args and call appropriate functions.
if __name__ == '__main__':

    import sys
    import os
    from dlreutil import readpipe,error

    # the program name is always the first argument at sys.argv[0].
    #argv = sys.argv[1:]
    
    #argv=sys.argv
    
    #argv[1] = 
    
    # Parse cmdargs
    
    
    fpath='/Users/elliotsefton-nash/workspace/dlre/test_geom_basic_fds_10_20_20_30.out'
    
    # Read pipe file    
    dat=readpipe(fpath)

    x=dat.clon.values
    y=dat.clat.values
    z=dat.tb.values
    orbit=dat.orbit.values

    xmin=10
    xmax=20
    ymin=20
    ymax=30
    res=1./128
    xstep=res
    ystep=res

    nhood=5
    densitythresh=0.6

    # Exclude points off grid
    mask = (x >= xmin) & (x < xmax) & (y >= ymin) & (y < ymax)
    x = x[mask]
    y = y[mask]
    z = z[mask]
    
    # Build bin edge vectors.
    xedges=np.arange(xmin,xmax+xstep,xstep)
    yedges=np.arange(ymin,ymax+ystep,ystep)
    
    avgInterp,cnt,err,zminInterp,zmaxInterp = bin2dbyorbit(x,y,z,orbit,xedges,yedges,res,nhood,densitythresh)    
    