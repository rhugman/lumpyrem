class TimeSeries():
    def __init__(self,ts_file):
        self.ts_file = ts_file
        self.ts_names = []
        self.scales = []
        self.offsets = []
        self.methods = []
        self.lr_models = [] # list of Model objects
        self.div_delta = 'div_by_delta'

    def write_ts(self, ts_file):
        #number of columns to include in the ts file
        count = len(self.ts_names)

        with open(ts_file, 'w') as f:
            print(self.lr_models)
            for model in self.lr_models:
                model_name = model.lumprem_model_name
                print(model_name)
                f.write('READ_LUMPREM_OUTPUT_FILE lr_'+model_name+'.out '+str(count)+'\n')
                f.write('#  my_name     LUMPREM_name      divide_by_delta_t?\n\n')
                #f.write(text+'\n')
            f.write('# Begin writting MF6 time series\n\n')
            f.write('\n')
        f.close()
