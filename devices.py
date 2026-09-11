# -*- coding: utf-8 -*-
"""
Created on Sun Dec 22 23:55:52 2019

@author: CCHIA
"""


import numpy as np
import gdspy as gp

# import self-defined functions
import holes as h
from geom import add_tapsup_to_beam
from layerdefs import wgspec, clmpspec

# units are in um

def dev_clp_phc_tsp_twg(cellname,gdslib,a,hx,hy,maxdef,oblong,taper_func,
                        ndef,nmirL,nmirR,nwgm,aW,hxW,hyW,wg_wid,
                        clp_wid,clp_strlen,clp_taplen,clp_outlen,
                        out_len,tapsup_len,tapsup_wid,ctr_len,in_wid,in_len,
                        end_len,end_wid,
                        phc_len=None,wid_corr=0,phc_scl=1,holeinctr=1):
    
    """ function to create device with following components (L-to-R):
        square + trapezoidal clamp; photonic crystal; tapered support; 
        tapered waveguide
        
        arguments:
        cellname: (string) name of cell
        gdslib: gdspy library class object to add cell to
        a: nominal lattice constant; hx, hy: nominal hole x-/y-widths
        maxdef: max lattice constant defect %; oblong: hole ellipticity factor
        taper_func: taper function ('cubic','quadratic','linear',None)
        ndef,nmirL,nmirR,a,hx,hy,maxdef,oblong,taper_func,nwgm,aW,hxW,hyW,
            wg_wid: define photonic crystal width and holes
        clp_wid,clp_strlen,clp_taplen,clp_outlen: define clamp
        out_len,tapsup_len,tapsup_wid,ctr_len,in_wid,in_len: define tapered support
        end_len,end_wid: define end waveguide taper
        wid_corr: width correction increment for angled etch lateral erosion
        phc_scl: local photonic crystal scaling
        holeinctr: 1/0 for hole/dielectric at center of photonic crystal
    """
    
    # create cell for device
    devcell = gdslib.new_cell(cellname, overwrite_duplicate=True, update_references=True)
#    devcell = gp.Cell(cellname)
    
    # layer specifications
#    wgspec = {'layer': 1, 'datatype': 0}    # waveguides
#    clmpspec = {'layer': 3, 'datatype': 0}  # clamps
    
    # apply local scaling to photonic crystal region
    a,hx,hy,wg_wid = np.array([a,hx,hy,wg_wid])*phc_scl
    
    # apply width correction to waveguide
    wg_wid,tapsup_wid,in_wid,end_wid = np.array([
            wg_wid,tapsup_wid,in_wid,end_wid])+wid_corr
    
    # generate hole geometry
    defH = h.calc_phc_holes(a,hx,hy,ndef,maxdef,oblong,taper_func)
    mirL = h.calc_phc_holes(a,hx,hy,nmirL,0,oblong,None)
    mirR = h.calc_phc_holes(a,hx,hy,nmirR,0,oblong,None)
    wgmR = h.calc_wgm_holes(a,hx,hy,nwgm,taper_func,
                          aT = aW, hxT = hxW, hyT = hyW)
    allholes = h.calc_phc_holes_lr((defH,mirL), (defH,mirR,wgmR), holeinctr)
    h.add_holes_to_cell(devcell,allholes)
    
    # get x-extents of holes
    x_exts = allholes['x_exts']
    
    # find x-start of clamp
    x0 = x_exts[0] - clp_strlen - clp_taplen - clp_outlen
    
    # length of photonic crystal section
    if phc_len is None:
        phc_len = x_exts[1] - x_exts[0]
    
    # start building device from clamp
    dev_wg = gp.Path(clp_wid,(x0,0))
    dev_wg.segment(clp_strlen, '+x', **clmpspec)
    dev_wg.segment(clp_taplen, '+x', final_width=wg_wid, **clmpspec)
    dev_wg.segment(clp_outlen, '+x', **wgspec)
    
    # waveguide
    dev_wg.segment(phc_len, '+x', **wgspec)
    
    # tapered support
    add_tapsup_to_beam(dev_wg,wg_wid,tapsup_wid,in_wid,out_len,tapsup_len,
                       ctr_len,in_len,**wgspec)
    
    # waveguide taper
    dev_wg.segment(0, '+x', final_width=2*dev_wg.w, **wgspec) 
    dev_wg.segment(end_len, '+x', final_width=end_wid, **wgspec)
    
    # get bounding box
    dev_box = dev_wg.get_bounding_box()  #[[x_min, y_min], [x_max, y_max]]
#    dev_wid = dev_bndbox[1,1] - dev_bndbox[0,1]
#    dev_len = dev_bndbox[1,0] - dev_bndbox[0,0]
    
    # add device to cell
    devcell.add(dev_wg)
    
    # add identifier text
#    label = gp.Label('dev_L'+str(nmirL)+'R'+str(nmirR),(0,2*wg_width),magnification=30)
#    devcell.add(label)
    
    outdev = {'cell': devcell,
              'box': dev_box,
              'phc_len': phc_len}
    
    return outdev

def cantilever(cellname,gdslib,wg_wid,beam_len,clp_wid,clp_strlen,clp_taplen,
               wid_corr=0):
    
    """ function to create device with following components (L-to-R):
        square + trapezoidal clamp; photonic crystal; tapered support; 
        tapered waveguide
        
        arguments:
        cellname: (string) name of cell
        gdslib: gdspy library class object to add cell to
        wg_wid: define width of waveguide
        clp_wid,clp_strlen,clp_taplen: define clamp
        wid_corr: width correction increment for angled etch lateral erosion
    """
    
    # create cell for device
    devcell = gdslib.new_cell(cellname, overwrite_duplicate=True, update_references=True)
    
    # apply width correction to waveguide
    wg_wid = wg_wid+wid_corr
    
    # find x-start of clamp
    x0 = -0.5*beam_len - clp_strlen - clp_taplen
    
    # start building device from clamp
    dev_wg = gp.Path(clp_wid,(x0,0))
    dev_wg.segment(clp_strlen, '+x', **clmpspec)
    dev_wg.segment(clp_taplen, '+x', final_width=wg_wid, **clmpspec)
    
    # waveguide
    dev_wg.segment(beam_len, '+x', **wgspec)
    
    # get bounding box
    dev_box = dev_wg.get_bounding_box()  #[[x_min, y_min], [x_max, y_max]]
    
    # add device to cell
    devcell.add(dev_wg)
    
    outdev = {'cell': devcell,
              'box': dev_box}
    
    return outdev