import numpy as np


class ArrayWithCallback(np.ndarray):
    def __new__(cls, input_array):
        obj = np.asarray(input_array).view(cls)
        if hasattr(input_array, '_callback_fns'):
            obj._callback_fns = input_array._callback_fns
        else:
            obj._callback_fns = list()
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self._callback_fns = getattr(obj, '_callback_fns', None)

    def add_modified_callback(self, fn):
        self._callback_fns.append(fn)

    def _call_modified_callbacks(self):
        if hasattr(self, '_callback_fns') and self._callback_fns is not None:
            for fn in self._callback_fns:
                fn()

    def __setattr__(self, key, value):
        super().__setattr__(key, value)
        if key != '_callback_fns':
            self._call_modified_callbacks()

    def __setitem__(self, key, item):
        super().__setitem__(key, item)
        self._call_modified_callbacks()


class DictWithCallback(dict):
    def __init__(self, *args, **kwargs):
        self._callback_fns = list()
        super().__init__(*args, **kwargs)

    def add_modified_callback(self, fn):
        self._callback_fns.append(fn)

    def _call_modified_callbacks(self):
        if hasattr(self, '_callback_fns'):
            for fn in self._callback_fns:
                fn()

    def __setattr__(self, key, value):
        super().__setattr__(key, value)
        if key != '_callback_fns':
            self._call_modified_callbacks()

    def __setitem__(self, key, item):
        super().__setitem__(key, item)
        self._call_modified_callbacks()


class LowercaseKeyDictWithCallback(DictWithCallback):
    def __setitem__(self, key, item):
        super().__setitem__(key.lower(), item)

    def __getitem__(self, item):
        return super().__getitem__(item.lower())


class ListWithCallback(list):
    def __init__(self, *args, **kwargs):
        self._callback_fns = list()
        super().__init__(*args, **kwargs)

    def add_modified_callback(self, fn):
        self._callback_fns.append(fn)

    def _call_modified_callbacks(self):
        if hasattr(self, '_callback_fns'):
            for fn in self._callback_fns:
                fn()

    def __setattr__(self, key, value):
        super().__setattr__(key, value)
        if key != '_callback_fns':
            self._call_modified_callbacks()

    def __setitem__(self, key, item):
        super().__setitem__(key, item)
        self._call_modified_callbacks()


def to_iter(obj):
    """
    Helper function to convert scalar objects to single element iterable objects.
    Parameters
    ----------
    obj : any

    Returns
    -------
    Iterable
        If object is iterable, return back the input.  If the input is a scalar, return back a 1 element list containing
        the object
    """
    try:
        iter(obj)
        return obj
    except TypeError:
        return [obj]