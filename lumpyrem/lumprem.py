def set_uniquekeys(dictionary):
    unique_keys = {'start_date':'01/01/2000',
                'end_date':'31/12/2010',
                'nday_out':'monthly',
                'steps_per_day':5
                }
    unique_keys.update(dictionary)
    return unique_keys



def write_lumprep(infile, unique_keys, datasets):
    #write the unique keys to input file
    with open(infile, 'w+') as f:
        f.write('# File created using lumpy. \n\n')
        for key in unique_keys.keys():
            f.write("{0: <32}{1:}".format(key.upper(), str(unique_keys[key])+'\n'))
        f.write('\n')
    f.close()

    #write datasets to input file
    count=0

    for dset in datasets:
        count = count+1
        with open(infile, 'a') as f:
            f.write('# Lumprem dataset number '+str(count)+'\n')
            for key in dset.keys():
                if type(dset[key]) == tuple:
                    input1 = dset[key][0]
                    input2 = dset[key][1]
                else:
                    input1 = dset[key]
                    input2 = ''
                    if input1 == None:
                        continue
                f.write("{0: <32}{1:}{2:}".format(key.upper(), str(input1),'\t'+str(input2)+'\n'))
            f.write('\n')
    f.close()


def write_ts_dict(ts_name ,scale=1.0, offset=0.0, method='linearend'):
    tsdict={}
    tsdict = {ts_name: {'scale':scale,'offset':offset,'method':method}}
    return tsdict

def write_cols(listoftuples):
    text = ''
    count = len(listoftuples)
    print(count)
    for i in range(count):
        lumprem_name = listoftuples[i][1]
        ts_name = listoftuples[i][0]
        divdelta = listoftuples[i][2]
        text += "{0: <16}{1: <16}{2: <16}\n".format(lumprem_name,ts_name,divdelta)
        print(text)
    return text, count



def write_lr2series(filein,ts_dict, col_dict, modellist):
    #write the unique keys to input file
    with open(filein, 'w') as f:
        for lr in modellist:
            text, count = write_cols(col_dict[lr])
            f.write('READ_LUMPREM_OUTPUT_FILE '+lr+'.out '+str(count)+'\n')
            f.write('#  my_name     LUMPREM_name      divide_by_delta_t?\n\n')
            f.write(text+'\n')
        f.write('# Begin writting MF6 time series\n\n')
        f.write('\n')
        for ts in ts_dict.keys():
            tsnum = len(ts_dict[ts])
            f.write('WRITE_MF6_TIME_SERIES_FILE '+ts+' '+str(tsnum)+'\n')
            f.write('# ts_name     scale        ofset        mf6_method\n\n')
            
            for tsname in ts_dict[ts].keys():
                f.write("  {0: <15}{1: <10}{2: <10}{3:}\n".format(tsname,ts_dict[ts][tsname]['scale'],'\t'+str(ts_dict[ts][tsname]['offset']),str(ts_dict[ts][tsname]['method'])))
            f.write('\n')
    f.close()



    # runs a process 
def run_process(process, path=False, commands=[], print_output=True):
        """This calls a process and then executes a list of commands
        process - str name of process. Process needs to be added to environmental path
        path - directory where input and output files are stored
        commands - list of strings with input commands in sequence. """
        import os
        if path == False:
            path = os.getcwd()
        
        owd = os.getcwd()
        os.chdir(path)

        import subprocess
        p = subprocess.run([process], stdout=subprocess.PIPE,
                input='\n'.join(map(str, commands))+'\n', encoding='ascii')

        if print_output==True:
                print(p.stdout)
        os.chdir(owd)



def get_col_dict(lrout):
    with open(lrout, 'r') as file:
        a = file.readline().split()
    file.close()
    coldict = dict.fromkeys(a,None)
    return coldict

def read_cols(lrm,tuple_list):
    dictionary={}
    dictionary[lrm] = tuple_list
    return dictionary

def get_modellist(datasets):
    modellist = []
    for dset in datasets:
        modellist.append('lr_'+dset['lumprem_model_name'])

def read_ouput(model,filename, dropna=True):
    import pandas as pd
    df = pd.read_csv(filename, delim_whitespace=True)
    df.dropna(inplace=True)
    df = df.apply(pd.to_numeric)
    return df




class Simulation:
    def __init__(self, model_list):
        self.start_date = '01/01/2000'
        self.end_date = '31/12/2010'
        self.nday_out = 'monthly'
        self.steps_per_day = 5
        self.model_list = model_list
        self.model_names = []
        for m in model_list:
            self.model_names.append(m.lumprem_model_name)

    def write_simulation(self, infile='lumprep.in', lumprep=True):
            #write the unique keys to lumprep input file
            unique_keys = self.__dict__.copy()
            unique_keys.pop('model_list', None)
            unique_keys.pop('model_names', None)

            with open(infile, 'w') as f:
                f.write('# File created using lumpyrem. \n\n')
                for key in unique_keys.keys():
                    f.write("{0: <32}{1:}".format(key.upper(), str(unique_keys[key])+'\n'))
                f.write('\n')
            f.close()

            #write models to lumprep input file
            count=0
            for model in self.model_list:
                model_dict = model.__dict__
                count = count+1
                with open(infile, 'a') as f:
                    f.write('# Lumprem dataset number '+str(count)+'\n')
                    for key in model_dict:
                        if type(model_dict[key]) == tuple:
                            input1 = model_dict[key][0]
                            input2 = model_dict[key][1]
                        else:
                            input1 = model_dict[key]
                            input2 = ''
                            if input1 == None:
                                continue
                        f.write("{0: <32}{1:}{2:}".format(key.upper(), str(input1),'\t'+str(input2)+'\n'))
                    f.write('\n')
            f.close()

            if lumprep == True:
                import os
                cwd = os.getcwd()
                run_process('lumprep', commands=[infile])
            else:
                return


    def run_simulation(self):
        for model in self.model_list:
            model_name = model.lumprem_model_name

            run_process('lumprem', commands=['lr_'+model_name+'.in','lr_'+model_name+'.out'])

    def get_results(self):
        import pandas as pd
        results = pd.DataFrame()
        for m in self.model_names:
            abc = []
            filename = 'lr_'+str(m)+'.out'
            textlist = []

            with open(filename) as f:
                for line in f:
                    textlist.append([i for i in line.split()])

            floatlist=[]
            for i in textlist[1:-2]:
                floatlist.append([float(x) for x in i])
            df = pd.DataFrame(floatlist, columns=textlist[0])
            df['model_name'] = str(m)
            results = pd.concat([results,df])

        param2 = pd.DataFrame()
        for model in self.model_list:
            columns = model.__dict__.keys()
            param = pd.DataFrame(columns=columns)
            for k in columns:
                param.at[0,k] = model.__dict__[k]
                param[k] = pd.to_numeric(param[k],errors='ignore')
            param2 = pd.concat([param2,param])

        final = results.merge(param2,how='outer', left_on='model_name', right_on='lumprem_model_name')
)
        return final


class Model():
    def __init__(self, model_name, silofile):
        self.silofile = (silofile,'evap')
        self.rainfile = 'rain.dat'
        self.epotfile = 'epot.dat'
        self.vegfile = (0.2, 1.5)
        self.irrigfile = (1, 1.0)
        self.maxvol = 0.5
        self.irrigvolfrac = 0.5
        self.rdelay = None
        self.mdelay = None
        self.ks = None
        self.M = None
        self.L = None
        self.mflowmax = None
        self.offset = None
        self.factor1 = None
        self.factor2 = None
        self.power = None
        self.vol = None
        self.lumprem_model_name = model_name
        self.batch_file = None
        self.pest_control_file = None