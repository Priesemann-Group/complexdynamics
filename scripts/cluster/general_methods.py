from argparse import ArgumentParser
import h5py
import yaml
import numpy as np
import os
from sympy import ShapeError

#from ESN import hyperparameter
#from ESN import ESN
#from ESN import high_dim_ESN


class progress_featback(object):
    """
    The class includes different featback functions which can be used to obtain featback for itterative aplications of similar steps.
    Some of which are reused from userbased function definitions.
    """
    def printProgressBar (iteration, total, prefix = '', suffix = '', decimals = 1, length = 100, fill = '█', printEnd = "\r",name=None):
        """
        Args:
            - iteration   - Required: current iteration (Int)
            - total       - Required: total iterations (Int)
            - prefix      - Optional: prefix string (Str)
            - suffix      - Optional: suffix string (Str)
            - decimals    - Optional: positive number of decimals in percent complete (Int)
            - length      - Optional: character length of bar (Int)
            - fill        - Optional: bar fill character (Str)
            - printEnd    - Optional: end character (e.g. "\r", "\r\n") (Str)

        return:
            - none
        notes:
            - The function returns in the last step similtaniously 100% and 100%-10^(-decimals)*1%.
            - The function is a modified version of the function from the surce:
            https://stackoverflow.com/questions/3173320/text-progress-bar-in-the-console (last visited 30.11.21).
        """
        percent = ("{0:." + str(decimals) + "f}").format(100 * (iteration / float(total)))
        filledLength = int(length * iteration // total)
        bar = fill * filledLength + '-' * (length - filledLength)
        if name is not None:
            print(f'\r{name}\t{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
        else:
            print(f'\r{prefix} |{bar}| {percent}% {suffix}', end = printEnd)
        # Print New Line on Complete
        if iteration == total-1: 
            percent = ("{0:." + str(decimals) + "f}").format(100)
            filledLength = length
            bar = fill * filledLength + '-' * (length - filledLength)
            if name is not None:
                print(f'\r{name}\t{prefix} |{bar}| {percent}% {suffix}')
            else:
                print(f'\r{prefix} |{bar}| {percent}% {suffix}')
        
    def print_percent_of_saved_frames(i, n):
        """
        Args:
            - iteration   - Required  : current iteration (Int)
            - total       - Required  : total iterations (Int)
        return:
            - none
        """
        if (100*i % n) == 0:
            print(f'Saved {100*i/n} % of the frames.', end = printEnd)

class yml:
    def read(PATH):
        """
        Opens a yaml file. 

        Args:   
            - PATH       -Required: Path of file to be loaded (str)

        Returns:
            - yaml_file if PATH variable is Path to yaml file, else raises Typeerror
        """
        try:
            with open(PATH, 'r') as file:
                yaml_file = yaml.safe_load(file)
        except FileNotFoundError:
            raise FileNotFoundError(
                PATH+" not found. Provide a path to an existing '.yml' File"
            )
            pass
        return yaml_file

    def create(PATH,dict=None):
        """
        Creates a yaml file from a dictionary

        Args: 
            - PATH       -Required: Path (with name) where the yaml file is saved(str)
            - dict       -Optional: Data to store in $PATH/$NAME.yml(dict)
        """
        try:
            with open(PATH, 'w') as file:
                yaml.dump(dict,file)
        except FileNotFoundError:
            raise FileNotFoundError(
                PATH+" not found. Provide an existing Path."
            )
    
    def change(PATH, action, variable=None, replacement = None, NewSavingPath=None):
        """
        changes an existing yaml file. 

        Args:   
            - PATH           -Required: Path of file to be loaded and edited

            - action         -Required:  'c' : change an existing variable
                                        'a' : add new variable
                                        'rm': remove existing variable

            - variable       - Optional: Ordered Array of names to the variable to work on. (str array)

            - replacement   - Optional: For add or change the entrie will be written in dict[variable]
                                        recommended type: dict or value (Int,float,str)

            - NewSavingPath - Optional: If not None, the changes are saved in a new file with the here saved path (+name)


        Returns:
            -yaml_file if PATH variable is Path to yaml file, else raises Typeerror

        Notes 
            -possebly way faster, files are overwritten with new format
        """
        try:
            with open(PATH, 'r+') as file:
                yaml_file = yaml.safe_load(file)
        
        except FileNotFoundError or PATH[len(PATH)-4:]=='.yml':
            raise FileNotFoundError(
                PATH+" not found. Provide a path to an existing '.yml' File"
            )
            pass
        if action=='c':
            if np.shape(variable)==(1,):
                yaml_file[variable[0]]=replacement
            elif np.shape(variable)==(2,):
                yaml_file[variable[0]][variable[1]]=replacement
            elif np.shape(variable)==(3,):
                yaml_file[variable[0]][variable[1]][variable[2]]=replacement
            elif np.shape(variable)==(4,):
                yaml_file[variable[0]][variable[1]][variable[2]][variable[3]]=replacement
            else:
                raise TypeError(
                    "Possible Errors follow:\n"
                    "'variable' and 'replacement' need to be defined for this action.\n"
                    "'variable' must be of type 'str array'\n"
                    "Only for levels of substructures have been added in the code. More can be included in yamlutils.py"
                )
            if NewSavingPath is None:
                yml.create(PATH,yaml_file)
            else:
                yml.create(NewSavingPath,yaml_file)
            return yaml_file

        elif action=='a':
            if np.shape(variable)==(1,):
                if variable[0] in yaml_file:
                    raise ValueError(
                        "The variable exitst already. Use 'c' to change it"
                    )
                else:
                    yaml_file[variable[0]]=replacement
            elif np.shape(variable)==(2,):
                if variable[0] in yaml_file:
                    if variable[1] in yaml_file[variable[0]]:
                        raise ValueError(
                            "The variable exitst already. Use 'c' to change it"
                        )
                    else:
                        yaml_file[variable[0]][variable[1]]=replacement
                else: 
                    yaml_file[variable[0]]={variable[1]:replacement}
            elif np.shape(variable)==(3,):
                if variable[0] in yaml_file:
                    if variable[1] in yaml_file[variable[0]]:
                        if variable[2] in yaml_file[variable[0]][variable[1]]:
                            raise ValueError(
                                "The variable exitst already. Use 'c' to change it"
                            )
                        else:
                            yaml_file[variable[0]][variable[1]][variable[2]]=replacement
                    else:
                        yaml_file[variable[0]][variable[1]]={variable[2]: replacement}
                else: 
                    yaml_file[variable[0]]={variable[1]:{variable[2]:replacement}}
            elif np.shape(variable)==(4,):
                if variable[0] in yaml_file:
                    if variable[1] in yaml_file[variable[0]]:
                        if variable[2] in yaml_file[variable[0]][variable[1]]:
                            if variable[3] in yaml_file[variable[0]][variable[1]][variable[2]]:
                                raise ValueError(
                                    "The variable exitst already. Use 'c' to change it"
                                )
                            else:
                                yaml_file[variable[0]][variable[1]][variable[2]][variable[3]] = replacement
                        else:
                            yaml_file[variable[0]][variable[1]][variable[2]]={variable[3]: replacement}
                    else:
                        yaml_file[variable[0]][variable[1]]={variable[2]: {variable[3]: replacement}}
                else: 
                    yaml_file[variable[0]]={variable[1]:{variable[2]: {variable[3]: replacement}}}
            else:
                raise TypeError(
                    "Possible Errors follow:\n"
                    "'variable' and 'replacement' need to be defined for this action.\n"
                    "'variable' must be of type 'str array'\n"
                    "Only for levels of substructures have been added in the code. More can be included in general_methods.py"
                )
            if NewSavingPath is None:
                yml.create(PATH,yaml_file)
            else:
                yml.create(NewSavingPath,yaml_file)
        
        elif action=='rm':
            if np.shape(variable)==(1,):
                del yaml_file[variable[0]]
            elif np.shape(variable)==(2,):
                del yaml_file[variable[0]][variable[1]]
            elif np.shape(variable)==(3,):
                del yaml_file[variable[0]][variable[1]][variable[2]]
            elif np.shape(variable)==(4,):
                del yaml_file[variable[0]][variable[1]][variable[2]][variable[3]]
            else:
                raise TypeError(
                    "Possible Errors follow:\n"
                    "'variable' and 'replacement' need to be defined for this action.\n"
                    "'variable' must be of type 'str array'\n"
                    "Only for levels of substructures have been added in the code. More can be included in yamlutils.py"
                )
            if NewSavingPath is None:
                yml.create(PATH,yaml_file)
            else:
                yml.create(NewSavingPath,yaml_file)

        else:
            raise TypeError(
                "Variable 'action' is '"+action+"' but needs to be 'c' to change a variable, 'a' to add a variable or 'rm' to rm a variable."
            )
        
        return yaml_file



def _get_path():
    """
    returns current Path
    """
    return os.path.dirname(os.path.abspath(__file__))+'/'