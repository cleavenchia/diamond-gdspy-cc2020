# -*- coding: utf-8 -*-
"""
Created on Sun Feb 17 23:33:58 2019

@author: CCHIA

Contains functions related to photonic crystal holes creation
"""

import numpy as np
import gdspy as gp

def taper_function(x, func_name='cubic'):
    """returns array for taper function, given index and function name
    
    x: index between 0 and 1
    func_name: 'cubic','quadratic','linear',None
    """
    try:
        if func_name == 'quadratic':
            f = lambda x: 1 - x**2
        elif func_name == 'cubic':
            f = lambda x: 1 - 3*x**2 + 2*x**3
        elif func_name == 'linear':
            f = lambda x: 1 - x
        elif func_name is None:
            f = lambda x: np.ones(np.size(x))
        return f(x)
    except:
        print('invalid taper function specified')

def calc_phc_holes(a,hx,hy,nholes,maxdef,oblong,taper_func,
                     x0 = 0,y0 = 0,**kwargs):
    """returns array of PhC hole sizes and positions for half of cavity
    
    arguments:
    a: nominal lattice constant; hx, hy: nominal hole x-/y-widths
    nholes: no. of holes; maxdef: max lattice constant defect %
    oblong: hole ellipticity factor; 
    taper_func: taper function ('cubic','quadratic','linear',None)
    x0, y0: x-/y-shifts for holes
    
    kwargs: 
    hxT, hyT, aT: predefined hole sizes, lattice constant to taper to
    
    outputs: 
    hole_params: dictionary with following keys: 
        hx_hole, hy_hole, x_hole, y_hole, a_hole
        - arrays of hole x-/y-widths, x-/y-coordinates, lattice constants
        x_exts - x-min and max extents for hole array
    
    for constant hole size, set maxdef=0 and taper_func=None
    """
    
    # hxT, hyT, aT - predefined hole sizes and lattice constant to taper to,
    # otherwise taper using maxdef and oblong by default
    if 'hxT' in kwargs and 'hyT' in kwargs and 'aT' in kwargs:
        hxT = kwargs.get('hxT')
        hyT = kwargs.get('hyT')
        aT = kwargs.get('aT')
    else:
        hxT = 0
        hyT = 0
        aT = 0
    
    # create arrays of hole sizes and lattice constants
    # for symmetric cavity with hole in center
    # 1st hole is central hole
    # lattice constant defined using hole separation, i.e. a_i = x_i - x_(i-1)
    # --> a_hole will have 1 less element compared to hx_hole, hy_hole
    x = np.arange(nholes)/nholes
    xa = np.arange(nholes-1)/(nholes-1)
    hx_hole = hxT + (hx-hxT)*(1-maxdef*taper_function(x,taper_func))**(1-oblong)
    hy_hole = hyT + (hy-hyT)*(1-maxdef*taper_function(x,taper_func))**(1+oblong)
    a_hole  = aT  + (a-aT)  *(1-maxdef*taper_function(xa,taper_func))
    
    # set first lattice constant to 0 for definition of x-positions
    # --> now a_hole will have same no. of elements as hx_hole, hy_hole
    a_hole = np.concatenate((np.array([0]),a_hole))

#    # old implementation:    
#    # lattice constant undefined for 1st hole - replaced with 0 for taper, 
#    #    or with fixed lattice constant if taper_func = None
#    # - only for generating x-coordinates
#    if taper_func is None:
#        a_1 = np.array([a])
#    else:
#        a_1 = np.array([0])
    
    # generate (x,y)-coordinates of holes
    x_hole = x0 + np.cumsum(a_hole)
    y_hole = y0 + np.zeros(np.size(hy_hole))
    
    # reset first lattice constant for subsequent geometry construction
    a_hole[0] = aT  + (a-aT)  *(1-maxdef*taper_function(0,taper_func))
    
    # get extents
    x_exts = np.array([x_hole[0]  - 0.5*a_hole[0],
                       x_hole[-1] + 0.5*a_hole[-1]])
    
    # create dictionary for hole parameters
    hole_params = {'hx_hole': hx_hole,
                   'hy_hole': hy_hole,
                   'x_hole' : x_hole,
                   'y_hole' : y_hole,
                   'a_hole' : a_hole,
                   'x_exts' : x_exts}
    
    return hole_params

def calc_wgm_holes(a,hx,hy,nholes,taper_func,
                     x0 = 0,y0 = 0,**kwargs):
    """returns array of waveguide mirror taper hole sizes and positions
    
    arguments:
    a: nominal lattice constant; hx, hy: nominal hole x-/y-widths
    nholes: no. of holes; 
    taper_func: taper function ('cubic','quadratic','linear',None)
    x0, y0: x-/y-shifts for holes
    
    kwargs: 
    hxT, hyT, aT: predefined hole sizes, lattice constant to taper to
    
    outputs: 
    hole_params: dictionary with following keys: 
        hx_hole, hy_hole, x_hole, y_hole, a_hole
        - arrays of hole x-/y-widths, x-/y-coordinates, lattice constants
        x_exts - x-min and max extents for hole array
    
    * output hole array does not include nominal mirror hole!
    """
    
    # hxT, hyT, aT - predefined hole sizes and lattice constant to taper to,
    # otherwise taper to zero by default
    if 'hxT' in kwargs and 'hyT' in kwargs and 'aT' in kwargs:
        xT = 1;
        hxT = kwargs.get('hxT')
        hyT = kwargs.get('hyT')
        aT = kwargs.get('aT')
    else:
        hxT = 0
        hyT = 0
        aT = 0
        xT = 0;
    
    # create arrays of hole sizes and lattice constants
    # for symmetric cavity with hole in center
    # 1st hole is central hole
    # lattice constant defined using hole separation, i.e. a_i = x_i - x_(i-1)
    x = np.arange(1,nholes+1)/(nholes-xT+1)
    hx_hole = hxT + (hx-hxT)*taper_function(x,taper_func)
    hy_hole = hyT + (hy-hyT)*taper_function(x,taper_func)
    a_hole  = aT  + (a-aT)  *taper_function(x,taper_func)
    
    # generate (x,y)-coordinates of holes
    x_hole = x0 + np.cumsum(a_hole)
    y_hole = y0 + np.zeros(np.size(hy_hole))
    
    # get extents
    x_exts = np.array([x_hole[0]  - 0.5*a_hole[0],
                       x_hole[-1] + 0.5*a_hole[-1]])
    
    # create dictionary for hole parameters
    hole_params = {'hx_hole': hx_hole,
                   'hy_hole': hy_hole,
                   'x_hole' : x_hole,
                   'y_hole' : y_hole,
                   'a_hole' : a_hole,
                   'x_exts' : x_exts}
    
    return hole_params

def calc_phc_holes_lr(left_holes, right_holes, holeinctr=1, s=None):
    """assemble full list of holes given segments of phc
    
    arguments:
    left_holes: tuple of dictionaries with hole parameters 
        for left half of photonic crystal 
        - will invert sign of x-coordinates
        - left_holes = None will create symmetric crystal using parameters 
          from right_holes
    right_holes: tuple of dictionaries with hole parameters 
        for right half of photonic crystal 
    holeinctr: 1 for hole in center of photonic crystal 
        (will remove center defect in left_holes)
        0 adds dielectric in center of photonic crystal
    s: s-parameter =) (see e.g. McCutcheon OE 2007)
        
    outputs: 
    hole_params: dictionary with following keys: 
        hx_hole, hy_hole, x_hole, y_hole
        - arrays of hole x-/y-widths, x-/y-coordinates
        x_exts - x-min and max extents for hole array
    
    create left_holes, right_holes with calc_phc_holes or calc_wgm_holes
    order of holes in left_holes will be reversed (due to sign flip in x)
    e.g. specifying left_holes = (h1,h2,h3) and right_holes = (h4,h5,h6),
        where hi = segments of holes, will give final hole array as
        (h3, h2, h1, h4, h5, h6) going from left to right of phc
    x = 0 will be position of center defect hole in right_holes
    """
    
    # assemble full array of right_holes
    hxR = np.array([])
    hyR = np.array([])
    aR = np.array([])
    yR = np.array([])
    for right_hole_params in right_holes:
        hxR = np.concatenate((hxR,right_hole_params['hx_hole']))
        hyR = np.concatenate((hyR,right_hole_params['hy_hole']))
        aR  = np.concatenate((aR, right_hole_params['a_hole']))
        yR  = np.concatenate((yR, right_hole_params['y_hole']))
    
    # calculating x-coordinates: set center hole to x = 0
    aR_init = aR[0]
    aR[0] = 0
    xR = np.cumsum(aR)
    
    # assemble full array of left_holes
    if left_holes is None:
        hxL = hxR
        hyL = hyR
        aL = aR
        yL = yR
    else:
        hxL = np.array([])
        hyL = np.array([])
        aL = np.array([])
        yL = np.array([])
        for left_hole_params in left_holes:
            hxL = np.concatenate((hxL,left_hole_params['hx_hole']))
            hyL = np.concatenate((hyL,left_hole_params['hy_hole']))
            aL  = np.concatenate((aL, left_hole_params['a_hole']))
            yL  = np.concatenate((yL, left_hole_params['y_hole']))
    
    # calculating x-coordinates: set center hole to x = 0
    # if hole not in center, shift center hole to left by a lattice constant
    aL_init = aL[0]
    aL[0] = 0
    if s is not None and holeinctr == 0:
        xL = -1*np.cumsum(aL) - s
    else:
        xL = -1*np.cumsum(aL) + (holeinctr-1)*0.5*(aL_init+aR_init)
    
    # assemble full array of holes
    hx_hole = np.concatenate((np.flipud(hxL[holeinctr:]),hxR))
    hy_hole = np.concatenate((np.flipud(hyL[holeinctr:]),hyR))
    x_hole  = np.concatenate((np.flipud(xL[holeinctr:]), xR))
    y_hole  = np.concatenate((np.flipud(yL[holeinctr:]), yR))
    
    # get extents
    x_exts = np.array([xL[-1]-0.5*aL[-1],xR[-1]+0.5*aR[-1]])
    
    # create dictionary for hole parameters
    hole_params = {'hx_hole': hx_hole,
                   'hy_hole': hy_hole,
                   'x_hole' : x_hole,
                   'y_hole' : y_hole,
                   'x_exts' : x_exts}
    
    return hole_params
    
def add_holes_to_cell(cell,hole_params,layer=2,npts=49):
    """add holes to cell given generated hole parameters
    
    inputs:
    cell: gdspy.Cell object to add holes to
    hole_params: hole parameters generated by calc_..._holes_... functions
    layer: GDSII layer to add holes to (default: 2)
    npts: no. of vertices to use for holes (should be 4n+1 for some integer n)
    outputs:
    """
    
    for hx_h, hy_h, x_h, y_h in zip(hole_params['hx_hole'],
                                    hole_params['hy_hole'],
                                    hole_params['x_hole'],
                                    hole_params['y_hole']):
        cell.add(gp.Round((x_h,  y_h),
                          (hx_h*0.5, hy_h*0.5), 
                          number_of_points = npts,
                          layer = layer))
    
