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
    def __init__(self, model_list, start_date='01/01/2000',end_date = '31/12/2010', 
                    steps_per_day = 1,
                    nday_out ='monthly',
                    batch_file = 'run.bat',
                    pest_control_file = 'temp.pst',
                    silofile= ('silofile.txt','evap')):
        self.start_date = start_date
        self.end_date = end_date
        self.nday_out = nday_out
        self.steps_per_day = steps_per_day
        self.model_list = model_list
        self.model_names = []
        for m in model_list:
            self.model_names.append(m.lumprem_model_name)
        self.batch_file = batch_file
        self.pest_control_file = pest_control_file
        self.silofile = silofile
        print(self.silofile[0])
        print(self.silofile[1])

    def write_simulation(self, infile='lumprep.in'):
            #write the unique keys to lumprep input file
            unique_keys = ['start_date','end_date','nday_out','steps_per_day']
            selfdict = self.__dict__

            with open(infile, 'w') as f:
                f.write('# File created using lumpyrem. \n\n')
                for key in unique_keys:
                    f.write("{0: <32}{1:}".format(key.upper(), str(selfdict[key])+'\n'))
                f.write('\n')
            f.close()

            #write models to lumprep input file
            count=0
            for model in self.model_list:
                model_dict = model.__dict__
                count = count+1
        
                with open(infile, 'a') as f:
                    f.write('# Lumprem dataset number '+str(count)+'\n')
                    f.write("{0: <32}{1:}{2:}".format('SILOFILE', self.silofile[0],'\t'+str(self.silofile[1])+'\n'))

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
            with open(infile, 'a') as f:
                f.write("{0: <32}{1:}{2:}".format('BATCH_FILE_NAME',self.batch_file,'\n'))
                f.write("{0: <32}{1:}{2:}".format('PEST_CONTROL_FILE',self.pest_control_file,'\n'))
                f.write('\n')
            f.close()

            run_process('lumprep', commands=[infile])



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
        return final


class Model():
    def __init__(self, model_name, silofile=False):
        if silofile == True: #added to work with Instance class
            self.silofile = (silofile,'evap')
        self.rainfile = 'rain.dat'
        self.epotfile = 'epot.dat'
        self.vegfile = (0.2, 1.5)
        self.irrigfile = (1, 1.0)
        self.maxvol = 0.5
        self.irrigvolfrac = 0.5
        self.rdelay = 5.0
        self.mdelay = 1.0
        self.ks = 0.1
        self.M = 0.5
        self.L = 0.5
        self.mflowmax = 0.1
        self.offset = 0.0
        self.factor1 = 1.0
        self.factor2 = 1.0
        self.power = 0.0
        self.vol = self.maxvol/2
        self.lumprem_model_name = model_name

    
    def write_model(self, file=False, numdays=100, noutdays=None, nstep=1, mxiter=100, tol=1.0e-5, rbuf =[0.0], mbuf=[0.0]):
        """
        rbuf: the recharge delay buffer. List of volumes assumed to have left soil moisture on previous days. First lelement left soil moisture on preivous day. Second element two days previously, etc.
        mbuf: the macropore delay buffer. Same setup as rbuf.
        """
        if noutdays == None:
            noutdays = numdays

        if file == False:
            file = 'lr_'+self.lumprem_model_name+'.in'
        
        import numpy as np
        outdays = np.linspace(1,numdays,noutdays, dtype=int)

        with open(file, 'w') as f:
            #f.write('# File written using lumpyrem \n')
            f.write('* earth properties \n')
            f.write("{0: <4}{1:}{2:}".format(self.maxvol,self.irrigvolfrac,'\n'))
            f.write("{0: <4}{1:<4}{2:}".format(self.rdelay,self.mdelay,'\n'))
            f.write("{0: <4}{1:<4}{2:<4}{3:<4}{4:}".format(self.ks,self.M,self.L,self.mflowmax,'\n'))
            
            f.write('* volume to elevation\n')
            f.write("{0: <4}{1:<4}{2:<4}{3:<4}{4:}".format(self.offset,self.factor1,self.factor2,self.power,'\n'))
            
            f.write('* topographic surface\n')
            f.write("{0:}{1:}".format(self.offset,'\n'))

            f.write('* initial conditions\n')
            f.write("{0:}{1:}".format(self.vol,'\n'))
            f.write("{0: <4}{1:}{2:}".format(len(rbuf), len(mbuf),'\n'))
            for i in rbuf:
                f.write("{0:}{1:}".format(i,' '))
            f.write('\n')
            for i in mbuf:
                f.write("{0:}{1:}".format(i,' '))
            f.write('\n')
            
            f.write('* solution parameters\n')
            f.write("{0: <4}{1:<4}{2:<4}{3:}".format(nstep,mxiter,tol,'\n'))
            
            f.write('* timing information\n')
            f.write("{0: <4}{1:<4}{2:}".format(numdays,noutdays,'\n'))
            for i in outdays:
                f.write("{0:}{1:}".format(i,' '))
            f.write('\n')
            
            f.write('* data filenames\n')
            if type(self.vegfile) == tuple:
                f.write("{0: <4}{1:}{2:}".format(self.vegfile[0],self.vegfile[1],'\n'))
            else:
                f.write("{0:}{1:}".format(self.vegfile,'\n'))
            f.write("{0:}{1:}".format(self.rainfile,'\n'))
            f.write("{0:}{1:}".format(self.epotfile,'\n'))
            
            if type(self.irrigfile) == tuple:
                f.write("{0: <4}{1:}{2:}".format(self.irrigfile[0],self.irrigfile[1],'\n'))
            else:
                f.write("{0:}{1:}".format(self.irrigfile,'\n'))
    
    def run_model(self):
        model_name = self.lumprem_model_name
        run_process('lumprem', commands=['lr_'+model_name+'.in','lr_'+model_name+'.out'])

    def write_irigfile(self, irrigfile, numdays, irrig_start, fracyear, irrig_end=None, date_start='01/01/1900'):
        import datetime as dt
        from datetime import date
        def add_years(d, years):
            """Return a date that's `years` years after the date (or datetime)
            object `d`. Return the same calendar date (month and day) in the
            destination year, if it exists, otherwise use the following day
            (thus changing February 29 to March 1).

            """
            try:
                return d.replace(year = d.year + years)
            except ValueError:
                return d + (date(d.year + years, 1, 1) - date(d.year, 1, 1))
                

        date_start = dt.datetime.strptime(date_start, '%d/%m/%Y')
        irrig_start = date_start + dt.timedelta(days=int(irrig_start))
        irrig_end = irrig_start + dt.timedelta(days=int(365*fracyear))

        if irrig_start == 0:
            tsteps = [[1,1,0.5]]
        else:
            tsteps = [[1,0,0.0]]

        numdays = date_start+dt.timedelta(days=numdays)

        while irrig_start < numdays:
            tsteps.append([(irrig_start-date_start).days, 1,0.5])
            tsteps.append([(irrig_end-date_start).days, 0,0.0])
            irrig_start = add_years(irrig_start,years=1)
            irrig_end = add_years(irrig_end,years=1)

        
        # write the file
        with open(irrigfile,'w') as f:
            for row in tsteps:
                f.write("{0} {1} {2}{3}".format(row[0],row[1],row[2],'\n'))

        

