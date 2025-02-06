# https://github.com/jonasw247/TDA_applied_on_SMLM/blob/main/From_SMLM_to_machine_learning-Tutorial.ipynb

import os
#from tkinter.messagebox import NO
import numpy as np
import gudhi as gd

# normalize_counts = 1000, none
# If set to N, N blink events are selected for each cell randomly
# If nucleus has less counts it is not considered.
# This can be useful for the comparison of cells with a different count numbers.
# But be careful, different count numbers can be biologically important and should not always be normalized.

# cell_size_area = 100.0 None
# Normalize cells to size in um**2
# This can be useful for the comparison of cells with a different size. 
# But be careful, different sizes can be biologically important and should not always be normalized. 

## TODO: fix argument defaults
def run_SNAP_TDA(x,y, normalize_counts=1000, cell_size_area=None):

    # normalize counts
    if normalize_counts != None:

        input_data = np.vstack([np.array(x),np.array(y)])

        if input_data.shape[1] > normalize_counts + 1:
            mask = np.random.choice(np.arange( input_data.shape[1] ), normalize_counts, replace = False)
            input_data = input_data.T[mask].T

        
        if input_data.shape[1] < normalize_counts - 1:
            print('count number smaller than ', normalize_counts,'not used')
            output_data = 0
            return output_data

    if cell_size_area != None:

        ## TODO: FIX RELOAD OF DATA
        reloadData = np.genfromtxt(path, usecols= 6, delimiter=',').T
        area = reloadData[1]
        cell_size_area_nm = cell_size_area * 1000000
        input_data = input_data / area * cell_size_area_nm
                                 
    data_filtered = input_data.T

    # apply persistent homology with alpha complex
    acX = gd.AlphaComplex(points=data_filtered).create_simplex_tree()
    dgmX = acX.persistence()

    # get the hole / betti number dimension 1 values
    output_data = np.sqrt(acX.persistence_intervals_in_dimension(1))
    return output_data

output_data = run_SNAP_TDA(x,y)