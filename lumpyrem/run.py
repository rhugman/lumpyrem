import os
import subprocess

def run_process(process, path=False, commands=[], print_output=True):
        """This calls a process and then executes a list of commands
        INPUT:
            Process: String name of process to execute. Process needs to be added to environmental path or exe in same folder.
            Path: directory where input and output files are stored.
            Commands: list of strings with input commands in sequence.
        """
        #import os
        if path == False:
            path = os.getcwd()
        
        owd = os.getcwd()
        os.chdir(path)

        #import subprocess
        p = subprocess.run([process], stdout=subprocess.PIPE,
                input='\n'.join(map(str, commands))+'\n', encoding='ascii')

        if print_output==True:
                print(p.stdout)
        os.chdir(owd)











