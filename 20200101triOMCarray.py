# -*- coding: utf-8 -*-
"""
Created on Tue Feb 19 02:59:04 2019

@author: CCHIA
"""

import numpy as np
import gdspy as gp
import os
import copy


# import self-defined functions
from devices import dev_clp_phc_tsp_twg as omc
from markers import sq4_markers
from layerdefs import lablspec

# units are in um

# initialize library
gp.current_library = gp.GdsLibrary()    # clear current library
                                        # (method to be deprecated in future ver of gdspy)
gdslib = gp.GdsLibrary(name='gdslib')
gdsname = '20200101triOMCarraytest.gds'

""" geometry definitions
    ====================
    wg_wid: PhC waveguide width, clp_wid: clamp width,
    clp_strlen: clamp straight length, clp_taplen: clamp taper length
    clp_outlen: straight length from clamp to PhC
    out_len: straight length from PhC to tapered support
    tapsup_len: length of taper from widest to narrowest section
    tapsup_wid: max width of tapered support
    ctr_len: tapered support center straight length
    in_wid: width of tapered support at input waveguide
    in_len: straight length from tapered support to input waveguide
    end_len: length of end input waveguide, end_wid: width of tip 
    wid_corr: width correction for angled etch lateral erosion
    phc_scl: local photonic crystal scaling
"""

# nominal geometry definitions in dict geom0
wg_wid = 0.873
geom0 = {'gdslib':gdslib,
         'wg_wid':wg_wid,      'clp_wid':5,                 'clp_strlen':5,
         'clp_taplen':3,       'clp_outlen':2.5,            'out_len':2.5,
         'tapsup_len':12.5,    'tapsup_wid': wg_wid*1.3,    'ctr_len':0,
         'in_wid':wg_wid,      'in_len':0,                  'end_len':40,
         'end_wid':0.1,        'wid_corr':0.2,              'phc_scl':1.02,
         'a':0.637,         'hx':0.318,     'hy':0.613,     'maxdef':0.16269,
         'oblong':1.4802,   'holeinctr':1,  'taper_func':'cubic',
         'ndef':9,          'nmirL':10,     'nmirR':4,      'nwgm':5,
         'aW':0.486,        'hxW':0.248,    'hyW':0.171}

# device variations
nmirR_all = np.array([2,3,4])  #from min, min+1, ..., max-1
phc_scl_all = np.array([0.98,1.00,1.02,1.04])

""" array layout rationale in this script:
    ======================================
    columns L to R: increasing global device scaling
    rows B to T: replicates of device with fixed input mirror, 
                 then groups of replicates with increasing input mirrors
    devices in adjacent columns offset by 7.5um ( = row spacing/no. of cols)
    to allow for fiber coupling to tapered waveguide
"""

# find length of longest device, then use this to fix length of 
# all photonic crystals andset array spacing
geomtest = copy.deepcopy(geom0)
geomtest['nmirR'] = np.max(nmirR_all)
geomtest['phc_scl'] = np.max(phc_scl_all)
devtest = omc('omctest',**geomtest)
devtestbox = devtest['box']
dev_xmin = devtestbox[0,0]
maxlen = devtestbox[1,0] - devtestbox[0,0]
geom0['phc_len'] = devtest['phc_len']

# array parameters
col_sp = np.ceil(maxlen/5.0)*5.0    # spacing between columns, ceil to nearest 5um
row_sp = 30                         # spacing between rows
row_offset = 7.5                    # vertical offset between adjacent columns
nreps = 4                           # no. of replicates for each device
ncols = phc_scl_all.shape[0]
nrows = nmirR_all.shape[0]
yctr = 0.5*(nrows*nreps-1)*row_sp

# generate device array
devarr = gdslib.new_cell('devarr', overwrite_duplicate=True, update_references=True)
geom1 = copy.deepcopy(geom0)
for psi,ps in enumerate(phc_scl_all):
    pscellname = 'omc_ps'+str("{:03d}".format(int(ps*100))) #cell name
    pscell = gdslib.new_cell(pscellname, overwrite_duplicate=True, update_references=True)
    for nri,nr in enumerate(nmirR_all):
        # generate devices for each input mirror no. and PhC scaling
        geom1['nmirR'] = nr
        geom1['phc_scl'] = ps
        devname = 'omc_nr'+str(int(nr))+'_ps'+str("{:03d}".format(int(ps*100))) #cell name
        outdev = omc(devname,**geom1)
        
        # replicate same device going up within column
        [pscell.add(gp.CellReference(outdev['cell'], 
                                     (0, (nri*nreps + rep)*row_sp - yctr)))
            for rep in np.arange(nreps)]
        
    devarr.add(gp.CellReference(pscell,((psi-1.5)*col_sp,(psi-2)*row_offset)))

# add text labels
for nri,nr in enumerate(nmirR_all):
    # add mirror text labels to left of devices, with 20um spacing from clamp
    mirlbl1 = gp.Text('M'+str(int(nr)),15,
                      (dev_xmin-1.5*col_sp,-2*row_offset+(nri*nreps)*row_sp-yctr),**lablspec) 
    mirlblbox = mirlbl1.get_bounding_box()
    mirlbldx = mirlblbox[1,0] - mirlblbox[0,0]
    mirlbldy = mirlblbox[1,1] - mirlblbox[0,1]
    mirlbl1.translate(-mirlbldx-20,-mirlbldy*0.5)
    
    # fillet and offset text labels to ensure shapes overlap
    # (default font gives disconnected polygons)
    mirlbl1.fillet(0.2)
    mirlbl2 = gp.offset(mirlbl1,0.2,**lablspec) 
    devarr.add(mirlbl2)

for psi,ps in enumerate(phc_scl_all):
    # add scaling text labels above devices
    scllbl1 = gp.Text('S'+str(int(np.mod(ps*100,10))),15,
                      ((psi-1.5)*col_sp, 20+(nrows*nreps-1)*row_sp - yctr+row_offset),
                      **lablspec) 
    
    # fillet and offset text labels to ensure shapes overlap
    # (default font gives disconnected polygons)
    scllbl1.fillet(0.2)
    scllbl2 = gp.offset(scllbl1,0.2,**lablspec) 
    devarr.add(scllbl2)

# add markers at 4 corners (+/-225,+/-225)um
mks = sq4_markers('mks',gdslib,l_trbl=2,l_tlbr=1)
devarr.add(gp.CellArray(mks['cell'],2,2,(450,450),(-225,-225)))

# write to gds
gdslib.add(devarr,include_dependencies=True, overwrite_duplicate=True, 
           update_references=True)          # add main gds cell with sub-cells
gdslib.write_gds(gdsname)

# print gds info in console
print("GDS file saved in " + gdsname + "\n"
      "Top level cell is ")
for c in gdslib.top_level():
    print(c)
print("List of all cells: ")
allcells = gdslib.cells
for c,cv in allcells.items():
    print(c)

# preview gds
#gp.LayoutViewer()       # preview gds with gdspylayout viewer
os.startfile(gdsname)    # open gds with default program (KLayout recommended)