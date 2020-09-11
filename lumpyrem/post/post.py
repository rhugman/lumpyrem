def read_tsfile(file): 
    """Reads a MODFLOW6 timeseries file and return timeries as rec array.

    INPUT
    file: path to file.
   
    RETURNS
    a: rec array of timesteps and timeseries values"""

    import numpy as np

    #read names and methods
    start = 0
    textlist = []
    with open(filename) as f:
        for line in f:
            if 'END ATTRIBUTES' in line.upper():
                start = 0
            elif start:
                textlist.append([i for i in line.split()] )
            elif 'BEGIN ATTRIBUTES' in line.upper():
                start = 1

    tsnames = textlist[0][1:]
    methods = textlist[1][1:]
    tsnames.insert(0, 'time')
    tscount = len(tsnames)

    start = 0
    textlist = []
    with open(filename) as f:
        for line in f:
            if 'END TIMESERIES' in line.upper():
                start = 0
            elif start:
                #textlist.append([float(i) for i in line.split()] )
                textlist.append(tuple([float(i) for i in line.split()]))
            elif 'BEGIN TIMESERIES' in line.upper():
                start = 1

    a = np.array(textlist, dtype={'names':tsnames,
                                'formats':tscount*['f8']})

    return a


