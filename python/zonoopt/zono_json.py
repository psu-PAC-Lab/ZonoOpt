import json
from scipy import sparse
import numpy as np
from ._core import *

def sparse_to_dict(m):
    coo = m.tocoo()
    return {
        "row": coo.row.tolist(),
        "col": coo.col.tolist(),
        "data": coo.data.tolist(),
        "shape": coo.shape
    }

def dict_to_sparse(d):
    return sparse.coo_matrix((d["data"], (d["row"], d["col"])), shape=d["shape"]).tocsc()

def to_json(Z: HybZono, filename: str):
    """
    Write a HybZono object to a JSON file.

    Args:
        Z (HybZono): The HybZono object to serialize.
        filename (str): The name of the JSON file to write to.
    """
    if Z.is_hybzono():
        class_str = 'HybZono'
    elif Z.is_conzono():
        class_str = 'ConZono'
    elif Z.is_zono():
        class_str = 'Zono'
    elif Z.is_point():
        class_str = 'Point'
    elif Z.is_empty_set():
        class_str = 'EmptySet'
    else:
        raise ValueError('Unsupported HybZono subclass.')

    dict = {
        'class': class_str,
        'n': Z.get_n(),
        'Gc': sparse_to_dict(Z.get_Gc()),
        'Gb': sparse_to_dict(Z.get_Gb()),
        'c': Z.get_c().tolist(),
        'Ac': sparse_to_dict(Z.get_Ac()),
        'Ab': sparse_to_dict(Z.get_Ab()),
        'b': Z.get_b().tolist(),
        'zero_one_form': Z.is_0_1_form()
    }
    with open(filename, 'w') as f:
        json.dump(dict, f)

def from_json(filename: str) -> HybZono:
    """
    Read a HybZono object from a JSON file.

    Args:
        filename (str): The name of the JSON file to read from.

    Returns:
        HybZono: The deserialized HybZono object.
    """
    with open(filename, 'r') as f:
        dict = json.load(f)
    class_str = dict['class']
    n = dict['n']
    Gc = dict_to_sparse(dict['Gc'])
    Gb = dict_to_sparse(dict['Gb'])
    c = np.array(dict['c'])
    Ac = dict_to_sparse(dict['Ac'])
    Ab = dict_to_sparse(dict['Ab'])
    b = np.array(dict['b'])
    zero_one_form = dict['zero_one_form']

    if class_str == 'HybZono':
        return HybZono(Gc, Gb, c, Ac, Ab, b, zero_one_form)
    elif class_str == 'ConZono':
        return ConZono(Gc, c, Ac, b, zero_one_form)
    elif class_str == 'Zono':
        return Zono(Gc, c)
    elif class_str == 'Point':
        return Point(c)
    elif class_str == 'EmptySet':
        return EmptySet(n)
    else:
        raise ValueError('Unsupported HybZono subclass.')