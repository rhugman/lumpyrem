import os
from lumpyrem import lumprem
class TimeSeries():
    def __init__(self,ts_file, lr_models, ts_names,
                      lumprem_ouput_cols, 
                      div_delta_t=True, 
                      workspace=False):
        
        model_count = len(lr_models)
        col_count = len(ts_names)
        if col_count != len(lumprem_ouput_cols):
            print('ERROR! LUMPREM columns and timesereis names must be the same length.\n')
            return
        
        self.ts_file = ts_file
        self.scales = model_count*[1]
        self.offsets = model_count*[0]
        self.methods = model_count*['linearend']
        self.lr_models = lr_models

        if div_delta_t == True:
            self.div_delta = model_count*['div_delta_t']
        else:
            elf.div_delta = model_count*['no_div_delta_t']

        self.lumprem_ouput_cols = lumprem_ouput_cols
        self.ts_names = ts_names

        if workspace==False:
            self.workspace = os.getcwd()
        else:
            self.workspace = workspace
       

    def write_ts(self):
        #number of columns to include in the ts file
        count = len(self.ts_names)
        ts_file = os.path.join(self.workspace, self.ts_file+'.in')

        with open(ts_file, 'w') as f:
            for model in self.lr_models:
                model_name = model.lumprem_model_name
                f.write('READ_LUMPREM_OUTPUT_FILE lr_'+model_name+'.out '+str(count)+'\n')
                f.write('#  my_name     LUMPREM_name      divide_by_delta_t?\n\n')

                for col in range(count):
                    f.write("\t{0}\t\t{1}\t\t{2}".format(self.ts_names[col]+'_'+model_name, self.lumprem_ouput_cols[col],self.div_delta[col]+'\n'))
                f.write('\n\n')

            f.write('WRITE_MF6_TIME_SERIES_FILE '+self.ts_file+' '+str(count*len(self.lr_models))+'\n')
            f.write("#\t{0}\t\t{1}\t\t{2}\t\t{3}".format('ts_name','scale','offset','mf6method\n\n'))
            for model in self.lr_models:
                model_name = model.lumprem_model_name
                for col in range(count):
                        f.write("\t{0}\t\t{1}\t\t{2}\t\t{3}\t{4}".format(self.ts_names[col]+'_'+model_name, self.scales[col],self.offsets[col],self.methods[col], '#'+model_name+'\n'))

        f.close()
        print('MF6 timeseries file '+ts_file+' written to:\n'+ts_file)
        
        #write ts file
        filename = self.ts_file
        path = self.workspace
        lumprem.run_process('lr2series', commands=[filename+'.in'],path=path)