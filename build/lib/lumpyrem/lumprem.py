import pandas as pd
import numpy as np
import os
from lumpyrem import run
import datetime as dt
from datetime import date

class Model():
    """
    A class used to represent a lumprem.in file. It facilitates the generating LUMPREM models and handling outputs.
    Default values are provded for all parameters however the user is advised to update those pertinnent to their case.

    Attributes
    ----------
    model_name : str
        name for LUMPREM model. Must be 3 characters or less.
    rainfile : str or tuple, optional
        name of rain file (default 'ran.dat').
    epotfile : 
        name of evaporation file (default 'epot.dat').
    vegfile : str or tuple, optional
        name of vegfile file or tuple with crop factor and gamma (default (0.5, 1)).
    irrigfile : str or tuple, optional
        name of irrigation file or tuple with irrigation code and groundwater extraction fraction (default (0, 0.0)).
    maxvol : float, optional
        volume of soil moisture (default 0.5)
    vol : float, optional
        volume at beggining of model run (default maxvol/2)
    irrigvolfrac: float, optional
        fraction of soil moisture irrigation must maintain, between 0.0 and 1.0. (default 0.5)
    rdelay : int, optional
        delay in days between water draining from soil moisture and being assigned to grounwater recharge (default 5)
    mdelay : int, optional
        delay in days between water draining from soil moisture and being assigned to macropore grounwater recharge (default 1)
    ks : float,optional
        saturated hydraulic conductivity (default 0.1)
    L : float, optional
        pore connectivity parameter (default 0.5)
    M : float, optional
        parameter determining the shape of drainage rate vs stored water relationship (default 0.5)
    mflowmax : float, optional
        maximum macropore recharegg allowed per day (default 0.1)
    offset : float, optional
        an offset, used to adjust volume to elevation (default 0.0)
    factor 1: float, optional
        factor1 in LUMPREM, used to adjust volume to elevation (default 1.0)
    factor 2: float, optional
        factor2 in LUMPREM, used to adjust volume to elevation (default 1.0)
    power: float, optional
        power in LUMPREM, used to adjust volume to elevation (default 1.0)

    elevmin : float, optional
        (default 0.0)
    elevmax : float, optional
        (default 0.0)
    surface : float, optional
        (default 0.0)

    silofile : bool, optional
        This is to facilitate interaction with the LUMPREP simulation class (default False).
    workspace : path 
        path to workspace folder. Default is current working directory
    """

    def __init__(self, model_name,rainfile='rain.dat',epotfile='epot.dat',vegfile=(0.5, 1),irrigfile=(0,0.0),maxvol=0.5,irrigvolfrac=0.5,
                rdelay=5,mdelay=1,ks=0.1,M=0.5,L=0.5,mflowmax=0.1,offset=0.0,factor1=2.0,factor2=3.0,power=0.5,elevmin=-9999.0, elevmax=10000.0,surface=0.0, vol=False, silofile=False, workspace=False):

        if silofile == True: #added to work with Instance class
            self.silofile = (silofile,'evap')
        self.rainfile = rainfile
        self.epotfile = epotfile
        self.vegfile = vegfile
        self.irrigfile = irrigfile
        self.maxvol = maxvol
        self.irrigvolfrac = irrigvolfrac
        self.rdelay = rdelay
        self.mdelay = mdelay
        self.ks = ks
        self.M = M
        self.L = L
        self.mflowmax = mflowmax
        self.offset = offset
        self.factor1 = factor1
        self.factor2 = factor2
        self.power = power
        self.elevmin = elevmin
        self.elevmax = elevmax
        self.surface = surface
        if vol ==False:
            self.vol = self.maxvol/2
        else:
            self.vol = vol
        self.lumprem_model_name = model_name

        if workspace==False:
            self.workspace = os.getcwd()
        else:
            self.workspace = workspace

    
    def write_model(self, file=False, numdays=100, noutdays=None, nstep=1, 
                          mxiter=100, tol=1.0e-5, rbuf =[0.0], mbuf=[0.0], start_date=None, end_date=None, print_output=True):
        """ Writes the LUMPREM model input files. 
        Default values are provded for all parameters however the user is advised to update those pertinnent to their case.

        Parameters
        ----------
        file : str, optional
            file name to write (default is the model_name with the sufix 'lr_' and file extension '.in'. For example: 'lr_abc.in')
        numdays : int, optional
            number of days of the model run (default 100)
        noutdays : int, optional
            number of days for which output is desired (default None, results in all days being recoreded)
        nstep : int, optonal
            number of steps into which each day is divided for iterative soil moisture computation (default 1)
        mxiter : int, optional
            max iterations per step 9default 100)
        tol : float, optional
            convergence tolerance
        rbuf : list, optional
            the recharge delay buffer for intial conditions. List of volumes assumed to have left soil moisture on previous days. First element left soil moisture on previous day. Second element two days previously, etc.
        mbuf : list, optional
            the macropore delay bufferfor intial conditions. Same setup as rbuf.
        start_date : str, optional
            date on which simualtion starts in 'dd/mm/yyyy' format. If provided noutdays can be assigned as 'monthly' to get ouputs on calendar months intervals.
        end_date : str, optional
            date str in 'dd/mm/yyyy' format. Requires start_date to be provided. If end_date is provided, numdays is calculated as difference between start_date and end_date. If numdays is provided this value is superceded by the calculted value.       
        """

        if noutdays == None:
            noutdays = numdays
        
        if start_date != None:
            start_date = dt.datetime.strptime(start_date, '%d/%m/%Y')

            if end_date == None:
                end_date = start_date + dt.timedelta(days=numdays)
            else:
                end_date = dt.datetime.strptime(end_date, '%d/%m/%Y')
                numdays = (end_date-start_date).days

            if noutdays == 'monthly':
                outdays=[1]
                date = start_date + dt.timedelta(days=1)
                while date <= end_date:
                    if date.day == 1:
                        timestep = (date-start_date).days
                        outdays.append(timestep)
                    date += dt.timedelta(days=1)
                noutdays = len(outdays)

            else:
                outdays = np.linspace(1,numdays,noutdays, dtype=int)
            
        else:
            outdays = np.linspace(1,numdays,noutdays, dtype=int)


        if file == False:
            file = 'lr_'+self.lumprem_model_name+'.in'
        
        file = os.path.join(self.workspace,file)

        if not os.path.exists(self.workspace):
            os.makedirs(self.workspace)
            
        

        with open(file, 'w+') as f:
            #f.write('# File written using lumpyrem \n')
            f.write('* earth properties \n')
            f.write("{0: <4} {1:}{2:}".format(self.maxvol,self.irrigvolfrac,'\n'))
            f.write("{0: <4} {1:<4}{2:}".format(self.rdelay,self.mdelay,'\n'))
            f.write("{0: <4} {1:<4} {2:<4} {3:<4}{4:}".format(self.ks,self.M,self.L,self.mflowmax,'\n'))
            
            f.write('* volume to elevation\n')
            f.write("{0: <4} {1:<4} {2:<4} {3:<4} {4} {5}\n".format(self.offset,self.factor1,self.factor2,self.power,self.elevmin, self.elevmax))
            
            f.write('* topographic surface\n')
            f.write("{0:}{1:}".format(self.surface,'\n'))

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
        
        if print_output==True:
            print('LUMPREM model input file written to: \n'+file+'\n')
    
    def run_model(self, print_output=True):
        """Runs the LUMPREM on model.
        """
        model_name = self.lumprem_model_name
        path = self.workspace
        run.run_process('lumprem', commands=['lr_'+model_name+'.in','lr_'+model_name+'.out'],path=path, print_output=print_output)



    def write_irigfile(self, numdays=None, irrig_start=0, fracyear=0.3, irrigfile='irrig.in',irrig_end=None, date_start='01/01/1900', date_end = None):
        """ Writes an irrigation input file to be used by LUMPREM.

        Parameters
        ----------
        numdays : int
            total number of days for the simulation
        irrig_start : str, int
            first date or day on which irrigation starts. Can be positive number of days or str in 'dd/mm/yyyy' format. If str, then it must be later than start_date. If date str is used then start_date must be set.
        fracyear : float, optional
            fraction of the year during which irrigation occurs. Must be between 0.0 and 1.0 (default 0.3)
        
        irrig_end : str, optional
            optionally the first date on which irrigation ends can be provided instead of fracyear. Must be str in forrmat 'dd/mm/yyyy'(default None). Requires start_date and irrig_start in str date format as well.
        date_start : str, optional
            date on which simulation starts (default '01/01/1900'). This aids in accoutning for month lengths and leap years.
        date_end : str, optional
            date on which simulation ends (default None). If provided, this supercedes numdays. This aids in accoutning for month lengths and leap years.
        irrigfile : str, optional
            name of irrigation file to write (default 'irrig.in')
        """

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

        if type(date_end) == str:
            numdays = dt.datetime.strptime(date_end, '%d/%m/%Y')
        else:
            numdays = date_start+dt.timedelta(days=numdays)


        if type(irrig_start) == str:
            irrig_start = dt.datetime.strptime(irrig_start, '%d/%m/%Y')
        else:
            irrig_start = date_start + dt.timedelta(days=int(irrig_start))

        if irrig_end == None:
            irrig_end = irrig_start + dt.timedelta(days=int(365*fracyear))
        else:
            irrig_end = dt.datetime.strptime(irrig_end, '%d/%m/%Y')


        if irrig_start == 0:
            tsteps = [[1,1,0.5]]
        else:
            tsteps = [[1,0,0.0]]

        

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
        """ Reads the results from the LUMPREM model and returns a Dataframe with parameters and results.

        Returns
        -------
        final : DataFrame
            Pandas dataframe of  model results and parameters.
        """

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
