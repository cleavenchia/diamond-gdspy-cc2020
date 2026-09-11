# -*- coding: utf-8 -*-
"""
Created on Sun Feb 17 23:41:19 2019

@author: CCHIA

Contains functions related to nanobeam geometry construction
"""
import numpy as np
import gdspy as gp

def add_tapsup_to_beam(device_wg,w_out,w_sup,w_in,l_out,l_tap,l_ctr,l_in,
                       layer=1,datatype=0):
    """add taper support (tapsup) structure to beam defined by path object
    
    segments of tapsup: 
    -> straight segment over l_out 
    -> quadratic tapered segment from w_out to w_sup over l_tap
    -> straight segment over l_ctr
    -> quadratic tapered segment from w_sup to w_in over l_tap
    -> straight segment over l_in
    
    arguments:
    device_wg: device waveguide (gdspy.Path object) to add tapsup to
    w_out: width of tapsup from output end (e.g. PhC to waveguide taper)
    w_sup: widest width of tapsup
    w_in: width of tapsup from input end (e.g. where fiber couples to device)
    l_out: length of straight segment before output taper
    l_tap: length over which width of tapsup tapers
    l_ctr: length of straight segment between tapers
    l_in: length of straight segment after input taper
    
    """
    def tapsup_out(t):
        if t < 0.5:
            y = w_out +     t**2 * 2 * (w_sup-w_out) 
        else:
            y = w_sup + (t-1)**2 * 2 * (w_out-w_sup) 
        return y
    
    def tapsup_in(t):
        if t < 0.5:
            y = w_sup +    t**2 * 2 * ( w_in-w_sup) 
        else:
            y = w_in + (t-1)**2 * 2 * (w_sup-w_in ) 
        return y
    
    spec = {'layer': layer, 'datatype': datatype}
    
    # 1. straight segment over l_out 
    device_wg.segment(l_out, '+x', final_width=w_out, **spec) 
    
    # 2. quadratic tapered segment from w_out to w_sup over l_tap
    device_wg.parametric(curve_function=lambda t: (l_tap*t,0),
                         number_of_evaluations=201,max_points=102,
                         final_width=tapsup_out,**spec)
    
    # get center x-coordinate of tapered support
#    bbox = device_wg.get_bounding_box()
#    tapsup_xctr = bbox[1,1] + 0.5*l_ctr
        
    # 3. straight segment over l_ctr
    # use segment of length 0 to match between segment and parametric curve of 
    # different widths, since width of segment tags to that of previous segment
    # then define actual segment
    device_wg.segment(0, '+x', final_width=w_sup, **spec)
    device_wg.segment(l_ctr, '+x', final_width=w_sup, **spec)
    
    # define small segment with final taper input width
    # for some weird reason, width of final straight segment over l_in tags 
    # to width of straight segment over l_ctr
    device_wg.segment(0, '+x', final_width=w_in, **spec)
    
    # 4. quadratic tapered segment from w_sup to w_in over l_tap
    device_wg.parametric(curve_function=lambda t: (l_tap*t,0),
                         number_of_evaluations=201,max_points=102,
                         final_width=tapsup_in,**spec)
    
    # 5. straight segment over l_in
    # code defining segment of zero length removed, as it will form undesired
    # small segment tapering quickly from w_sup to w_in
    device_wg.segment(l_in, '+x', final_width=w_in, **spec)
    
#    if return_xctr is True:
#        return tapsup_xctr
#    else:
#        pass


def add_tethers_to_tapsup(cell,x_ctr,n_tets,w_tets,l_tets,g_tets,tets_side=0,
                          layer=1,datatype=0):
    """ add tethers to taper support (tapsup) structure for rectangular 
        cross-section devices
    
    arguments:
    cell: gdspy cell to add tethers to
    x_ctr: x-coordinate of center of tethers 
    n_tets, w_tets, l_tets: no., width, full length of tethers
    g_tets: gap between tethers
    tets_side: 1/-1 to add tethers on top/bottom (+/-y) half of device;
               0 to add to both sides
    """
    spec = {'layer': layer, 'datatype': datatype}
    
    if n_tets > 0:
        # calc top and bottom y-coordinates of end of tethers
        # for tets_side = -1/0/+1, (y_bot,y_top) = (-l,0)/(-l/2,l/2)/(0,l)
        y_top = (tets_side+1)*0.5*l_tets
        y_bot = (tets_side-1)*0.5*l_tets
        
        # add tethers
        # x-coordinate calulation: shift to x_ctr, move by n tethers, then +/- 
        # half of tether width
        [cell.add(gp.Rectangle(
                (x_ctr + (n-0.5*(n_tets-1))*g_tets - 0.5*w_tets, y_bot),
                (x_ctr + (n-0.5*(n_tets-1))*g_tets + 0.5*w_tets, y_top),
                **spec)) for n in np.arange(n_tets)]
    else:
        pass
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    