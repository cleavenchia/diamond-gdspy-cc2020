# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 13:00:22 2020

@author: CCHIA
"""
import gdspy as gp
from layerdefs import markspec

def sq4_markers(cellname,gdslib,l_trbl=2,l_tlbr=0.5):
    
    """ function to create 4-square alignment markers on GDS layout
        with small gap of 50x50nm at center
        
        arguments:
        cellname: (string) name of cell
        gdslib: gdspy library class object to add cell to
        l_trbl: side length of top-right and bottom-left squares in um
        l_tlbr: side length of top-left and bottom-right squares in um
    """
    
    # create cell to place markers in
    mkcell = gdslib.new_cell(cellname, overwrite_duplicate=True, update_references=True)
    
    # create square markers at: top-right (large), bottom-right (small),
    #                           bottom-left (large), top-left (small)
    # with small gap of 50x50nm at center
    rectTR = gp.Rectangle(( 0.025, 0.025),(   l_trbl,   l_trbl),**markspec)
    rectBR = gp.Rectangle(( 0.025,-0.025),(   l_tlbr,-1*l_tlbr),**markspec)
    rectBL = gp.Rectangle((-0.025,-0.025),(-1*l_trbl,-1*l_trbl),**markspec)
    rectTL = gp.Rectangle((-0.025, 0.025),(-1*l_tlbr,   l_tlbr),**markspec)
    
    mkcell.add(rectTR).add(rectBR).add(rectBL).add(rectTL)
    mkbox = mkcell.get_bounding_box()
    
    outmk = {'cell': mkcell,
             'box': mkbox}
    
    return outmk