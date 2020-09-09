# runs a process 
def run_process(process, path, commands=[], print_output=True):
        """This calls a process and then executes a list of commands
        process - str name of process. Process needs to be added to environmental path
        path - directory where input and output files are stored
        commands - list of strings with input commands in sequence. """

        import os
        owd = os.getcwd()
        os.chdir(path)

        import subprocess
        p = subprocess.run([process], stdout=subprocess.PIPE,
                input='\n'.join(map(str, commands))+'\n', encoding='ascii')

        if print_output==True:
                print(p.stdout)
        os.chdir(owd)











