import pandas as pd
import numpy as np
import os
from lumpyrem import run

class Model():
    def __init__(self, model_name, silofile=False, workspace=False):
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

        if workspace==False:
            self.workspace = os.getcwd()
        else:
            self.workspace = workspace

    
    def write_model(self, file=False, numdays=100, noutdays=None, nstep=1, 
                          mxiter=100, tol=1.0e-5, rbuf =[0.0], mbuf=[0.0]):
        """
        rbuf: the recharge delay buffer. List of volumes assumed to have left soil moisture on previous days. First lelement left soil moisture on preivous day. Second element two days previously, etc.
        mbuf: the macropore delay buffer. Same setup as rbuf.
        """
        if noutdays == None:
            noutdays = numdays

        if file == False:
            file = 'lr_'+self.lumprem_model_name+'.in'
        
        file = os.path.join(self.workspace,file)

        if not os.path.exists(self.workspace):
            os.makedirs(self.workspace)

        import numpy as np
        outdays = np.linspace(1,numdays,noutdays, dtype=int)

        with open(file, 'w+') as f:
            #f.write('# File written using lumpyrem \n')
            f.write('* earth properties \n')
            f.write("{0: <4} {1:}{2:}".format(self.maxvol,self.irrigvolfrac,'\n'))
            f.write("{0: <4} {1:<4}{2:}".format(self.rdelay,self.mdelay,'\n'))
            f.write("{0: <4} {1:<4} {2:<4} {3:<4}{4:}".format(self.ks,self.M,self.L,self.mflowmax,'\n'))
            
            f.write('* volume to elevation\n')
            f.write("{0: <4} {1:<4} {2:<4} {3:<4}{4:}".format(self.offset,self.factor1,self.factor2,self.power,'\n'))
            
            f.write('* topographic surface\n')
            f.write("{0:}{1:}".format(self.offset,'\n'))

            f.write('* initial conditions\n')
            f.write("{0:}{1:}".format(self.vol,'\n'))
            f.write("{0: <4} {1:}{2:}".format(len(rbuf), len(mbuf),'\n'))
            for i in rbuf:
                f.write("{0:}{1:}".format(i,' '))
            f.write('\n')
            for i in mbuf:
                f.write("{0:}{1:}".format(i,' '))
            f.write('\n')
            
            f.write('* solution parameters\n')
            f.write("{0: <4} {1:<4} {2:<4}{3:}".format(nstep,mxiter,tol,'\n'))
            
            f.write('* timing information\n')
            f.write("{0: <4} {1:<4}{2:}".format(numdays,noutdays,'\n'))
            for i in outdays:
                f.write("{0:}{1:}".format(i,' '))
            f.write('\n')
            
            f.write('* data filenames\n')
            if type(self.vegfile) == tuple:
                f.write("{0: <4} {1:}{2:}".format(self.vegfile[0],self.vegfile[1],'\n'))
            else:
                f.write("{0:}{1:}".format(self.vegfile,'\n'))
            f.write("{0:}{1:}".format(self.rainfile,'\n'))
            f.write("{0:}{1:}".format(self.epotfile,'\n'))
            
            if type(self.irrigfile) == tuple:
                f.write("{0: <4} {1:}{2:}".format(self.irrigfile[0],self.irrigfile[1],'\n'))
            else:
                f.write("{0:}{1:}".format(self.irrigfile,'\n'))
        
        print('LUMPREM model input file written to: \n'+file+'\n')
    
    def run_model(self):
        model_name = self.lumprem_model_name
        path = self.workspace
        run.run_process('lumprem', commands=['lr_'+model_name+'.in','lr_'+model_name+'.out'],path=path)



    def write_irigfile(self, numdays, irrig_start, fracyear, irrigfile='irrig.in',irrig_end=None, date_start='01/01/1900'):
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
        irrigfile = os.path.join(self.workspace,irrigfile)

        with open(irrigfile,'w') as f:
            for row in tsteps:
                f.write("{0} {1} {2}{3}".format(row[0],row[1],row[2],'\n'))
        print('Irrigation input file written to: \n'+irrigfile)

        
    def get_results(self):
        """Returns a Pandas dataframe of model results and parameters from the Simulation object.
        """
        import pandas as pd

        filename = 'lr_'+str(self.lumprem_model_name)+'.out'
        filename = os.path.join(self.workspace, filename)
        
        textlist = []
        with open(filename) as f:
            for line in f:
                textlist.append([i for i in line.split()])

        floatlist=[]
        for i in textlist[1:-2]:
            floatlist.append([float(x) for x in i])
        df = pd.DataFrame(floatlist, columns=textlist[0])
        df['lumprem_model_name'] = str(self.lumprem_model_name)

        columns = self.__dict__.keys()
        param = pd.DataFrame(columns=columns)
        for k in columns:
            param.at[0,k] = self.__dict__[k]
            param[k] = pd.to_numeric(param[k],errors='ignore')

        final = df.merge(param,how='outer', left_on='lumprem_model_name', right_on='lumprem_model_name')
        return final
