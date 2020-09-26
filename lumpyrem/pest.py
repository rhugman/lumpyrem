from lumpyrem import lumprem, run
import os

class Pst():
    """ A Pest setup class. Facilities generating PEST control, template and instruction files from an ennsemble of LUMPREM models.
    """
    def __init__(self, controlfile='temp.pst', workspace='.'):
        self.controlfile = os.path.join(controlfile)


    def write_pst(self, file=False):
        if file ==False:
            file = self.controlfile

        with open(file, 'w+') as f:
            f.write('pcf\n')