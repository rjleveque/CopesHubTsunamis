from numpy import ma
import os, glob, re
import matplotlib.pyplot as plt
import numpy as np

def interp(xi,yi,x,y,z):
    """
    Given 2d arrays x,y,z of data (on rectangular grid with uniform spacing),
    Interpolate to points specified by xi,yi and return zi.
    Uses bilinear interpolation, extrapolating if xi,yi outside limits of data

    If xi and yi are both floats, returns a float.
    If xi is a float and yi a 1-dimensional array (or vice versa), 
        returns a 1-d array
    If xi and yi both have the same shape (1-d or 2-d arrays), returns
        an array of the same shape.
    """

    X = x; Y = y; Z = z;
    if y[1,0] == y[0,0]:
        X = X.T
        Y = Y.T
        Z = Z.T
    if Y[1,0] < Y[0,0]:
        Y = np.flipud(Y)
        Z = np.flipud(Z)
    my,mx = X.shape
    xlow,ylow = X[0,0], Y[0,0]
    dx = X[0,1]-X[0,0]
    dy = Y[1,0]-Y[0,0]
    
    xoff = (xi-xlow)/dx
    yoff = (yi-ylow)/dy

    ii = xoff.astype(int)
    ii = np.where(ii > 0, ii, 0)
    ii = np.where(ii < mx-1, ii, mx-1)

    ji = yoff.astype(int)
    ji = np.where(ji > 0, ji, 0)
    ji = np.where(ji < my-1, ji, my-1)
    
    alphax = (xi - X[ji,ii]) / dx
    alphay = (yi - Y[ji,ii]) / dy
    w00 = (1-alphax)*(1-alphay)
    w01 = alphax*(1-alphay)
    w10 = (1-alphax)*alphay
    w11 = alphax*alphay
    zi = w00*Z[ji,ii] + w01*Z[ji,ii+1] + w10*Z[ji+1,ii] + w11*Z[ji+1,ii+1]
        
    return zi


def transect(xy1,xy2,x,y,z,npts=1001):
    xi = np.linspace(xy1[0], xy2[0], npts)
    yi = np.linspace(xy1[1], xy2[1], npts)
    zi = np.array([interp(xii,yii,x,y,z) for xii,yii in zip(xi,yi)])
    zi = interp(xi,yi,x,y,z)
    return xi,yi,zi

