from typing import Any
import os
import os.path
import sys
import sysconfig
import yaml

#------------------------------------------------------------------------------
#           SasktranifRegistry
#------------------------------------------------------------------------------

class SasktranifRegistry():

    #------------------------------------------------------------------------------
    #           __init__
    #------------------------------------------------------------------------------

    def __init__(self):
        self.registry = None
        self.load()

    #------------------------------------------------------------------------------
    #           load
    #------------------------------------------------------------------------------

    def load(self):
        if self.registry is None:
            if (os.path.exists(self.registry_filename())):
                with open(self.registry_filename(), 'rt') as f:
                    self.registry = yaml.safe_load(f)
            else:
                self.registry = {}


    #------------------------------------------------------------------------------
    #           save
    #------------------------------------------------------------------------------

    def save(self):
        if self.registry is not None:
            dirname = os.path.dirname(self.registry_filename())
            if (not os.path.exists(dirname)):
                os.makedirs(dirname)
            with open(self.registry_filename(), 'wt') as f:
                yaml.safe_dump( self.registry,f,default_flow_style=False)

    #------------------------------------------------------------------------------
    #           subkey
    #------------------------------------------------------------------------------

    def subkey(self, stringindex: str)->Any:
        stringindex = stringindex.lower()
        key = self.registry
        tokens = stringindex.split('/')
        for index in tokens:
            x = key.get(index)
            if x is None:
                key[index] = {}
                x = key[index]
            key = x
        return key

    #------------------------------------------------------------------------------
    #           setkey
    #------------------------------------------------------------------------------

    def setkey(self, keyname, value):
        keyname = keyname.lower()
        key = self.registry
        tokens = keyname.split('/')
        for index in tokens[:-1]:

            x = key.get(index)
            if x is None:
                key[index] = {}
                x = key[index]
            key = x
        key[tokens[-1]] = value
    #------------------------------------------------------------------------------
    #           registry_filename
    #------------------------------------------------------------------------------

    def registry_filename(self):

        filename = os.path.join(sysconfig.get_path('data'), 'share', 'usask-arg','registry','sasktranif','globalkey.yaml')
        return filename
