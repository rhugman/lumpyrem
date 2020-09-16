import os
import subprocess

def run_process(process, path=False, commands=[], print_output=True):
        """This calls a process and then executes a list of commands.
        
        Parameters
        ----------
        process : str
            The name of the process to execute.
        path : str, optional
            path in which to execute commands. False (default) result sin commans being executed in current working directory.
        commands : list of str
            sequence of commands to pass to the process.
        print_output : bool, optional
                True, process output is printed. False, it is not.
            """

        if path == False:
            path = os.getcwd()
        
        owd = os.getcwd()
        os.chdir(path)

        p = subprocess.run([process], stdout=subprocess.PIPE,
                input='\n'.join(map(str, commands))+'\n', encoding='ascii')

        if print_output==True:
                print(p.stdout)
        os.chdir(owd)











