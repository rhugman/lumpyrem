import pandas as pd
import numpy as np
import os
from lumpyrem import run

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
                    f.write("{0: <32} {1:}".format(key.upper(), str(selfdict[key])+'\n'))
                f.write('\n')
            f.close()

            #write models to lumprep input file
            count=0
            for model in self.model_list:
                model_dict = model.__dict__
                count = count+1
        
                with open(infile, 'a') as f:
                    f.write('# Lumprem dataset number '+str(count)+'\n')
                    f.write("{0: <32} {1:}{2:}".format('SILOFILE', self.silofile[0],'\t'+str(self.silofile[1])+'\n'))

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
                f.write("{0: <32} {1:}{2:}".format('BATCH_FILE_NAME',self.batch_file,'\n'))
                f.write("{0: <32} {1:}{2:}".format('PEST_CONTROL_FILE',self.pest_control_file,'\n'))
                f.write('\n')
            f.close()

            run.run_process('lumprep', commands=[infile])

    def run_simulation(self):
        """Runs LUMPREM on models created using LUMPREP in the Simulation object.
        """
        for model in self.model_list:
            model_name = model.lumprem_model_name
            run.run_process('lumprem', commands=['lr_'+model_name+'.in','lr_'+model_name+'.out'])

    def get_results(self):
        """Returns a Pandas dataframe of model results and parameters from the Simulation object.
        """
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
