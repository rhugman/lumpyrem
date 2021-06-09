import os
from lumpyrem import run
import numpy as np

class TimeSeries():
    """
    A class used to link a MODFLOW6 timeseries file to LUMPREM model outputs.

    Attributes
    ----------
    ts_file : str
        filename and path of MODFLOW6 timeseries file to write.
    lr_models : list
        list of lumpyrem Model objects.
    ts_names : list of str
        list of names for timeseries to include in the MODFLOW6 timeseries file.
    lumprem_ouput_cols : list of str
        list of LUMPREM output columns to import as timeseries. This list must match ts_names length and order.
    div_delta_t : bool
        True (Default) if LR2SERIES div_delta_t. False if LR2SERIES no_div_delta_t .
    workspace : path 
        Path to workspace folder. Default is current working directory.
    """

    def __init__(self,ts_file, lr_models, ts_names,
                      lumprem_output_cols,methods, 
                      div_delta_t=True, 
                      workspace=False, scales=None, timeoffset=' ', time_offset_method='next', sep='_', tssufix='modelname'):
        """Parameters
        ----------
        ts_file : str
            filename and path of MODFLOW6 timeseries file to write.
        lr_models : list
            list of lumpyrem Model objects.
        ts_names : list of str
            list of names for timeseries to include in the MODFLOW6 timeseries file.
        lumprem_ouput_cols : list of str
            list of LUMPREM output columns to import as timeseries. This list must match ts_names length and order.
        methods : list of str
            list of methdos to use in timeseries file. Must match sequence and length of ts_names.
        div_delta_t : bool or list
            True (Default) if LR2SERIES div_delta_t. False if LR2SERIES no_div_delta_t. Alternatively a list of str can be provided. It must match the length and sequence of ts_names.
        workspace : path 
            Path to workspace folder. Default is current working directory.
        scales : list of float
            list of floats to scale the lumprem outputs to. Must be in sequence and of same length as ts_names.
        tssufix: 'modelname', None or str
            adds sufix to the ts names. Use None to pass tsnames epxlicitly. 'modelname' appends the LUMPREM model name. TO DO: curently only works for a single model.
        """
        
        model_count = len(lr_models)
        col_count = len(ts_names)
        if col_count != len(lumprem_output_cols):
            print('ERROR! LUMPREM columns and timeseries names must be the same length.\n')
            return
        
        self.ts_file = ts_file
        if scales == None:
            self.scales = col_count*[1]
        else:
            self.scales = scales
        self.offsets = col_count*[0]
        self.methods = methods
        self.lr_models = lr_models

        if div_delta_t == True:
            self.div_delta = col_count*['div_delta_t']
        elif div_delta_t == False:
            self.div_delta = col_count*['no_div_delta_t']
        else:
            self.div_delta = div_delta_t

        self.lumprem_output_cols = lumprem_output_cols

        if tssufix=='modelname':
            self.ts_names = [i+sep+j.lumprem_model_name for i,j in zip(ts_names, col_count*lr_models)]
        elif tssufix==None:
            self.ts_names = ts_names
        else:
            self.ts_names = [i+sep+tssufix for i in ts_names]

        if workspace==False:
            self.workspace = os.getcwd()
        else:
            self.workspace = workspace
        self.timeoffset = timeoffset

        if timeoffset != ' ':
            self.time_offset_method = col_count*[time_offset_method]
        else:
            self.time_offset_method = col_count*['']
        
        self.sep = sep

    def write_ts(self):
        """Writes the MODFLOW6 timeseries file.

        Parameters
        ----------
        """
        #number of columns to include in the ts file
        count = len(self.ts_names)
        ts_file = os.path.join(self.workspace, self.ts_file+'.in')

        with open(ts_file, 'w') as f:
            for model in self.lr_models:
                model_name = model.lumprem_model_name
                f.write('READ_LUMPREM_OUTPUT_FILE lr_'+model_name+'.out '+str(count)+'\n')
                f.write('#  my_name     LUMPREM_name      divide_by_delta_t?\n\n')

                for col in range(count):
                    f.write("\t{0}\t\t{1}\t\t{2}".format(self.ts_names[col], self.lumprem_output_cols[col],self.div_delta[col]+'\n'))
                f.write('\n\n')

            f.write('WRITE_MF6_TIME_SERIES_FILE '+self.ts_file+' '+str(count*len(self.lr_models))+' '+str(self.timeoffset)+'\n')
            f.write("#\t{0}\t\t{1}\t\t{2}\t\t{3}\t\t{4}".format('ts_name','scale','offset','mf6method','time_offset_method\n\n'))
            for model in self.lr_models:
                model_name = model.lumprem_model_name
                for col in range(count):
                        f.write("\t{0}\t\t{1}\t\t{2}\t\t{3}\t{4}\t{5}".format(self.ts_names[col], self.scales[col],self.offsets[col],self.methods[col], self.time_offset_method[col], '#'+model_name+'\n'))

        f.close()
        print('MF6 timeseries file '+ts_file+' written to:\n'+ts_file)
        
        #write ts file
        filename = self.ts_file
        path = self.workspace
        run.run_process('lr2series', commands=[filename+'.in'],path=path)

def read_ts(filename):
    """Reads a modflow6 timeseries file and returns the timeseries as a rec array.

    Parameters
    ----------
    filename : str
        filename of ts file to read

    Returns
    -------
    a : numpy recarray 
        recarray with timesteps and timeseries value. Timeseries names are column names.
    tsnames : list of str
        list of timeseries names in the ts file
    """

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
    names = tsnames.copy()
    names.insert(0, 'time')
    methods = textlist[1][1:]
    count = len(names)

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


    a = np.array(textlist, dtype={'names':names,
                                'formats':count*['f8']})
    return a, tsnames, methods