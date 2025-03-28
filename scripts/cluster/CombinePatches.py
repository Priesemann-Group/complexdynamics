from argparse import ArgumentParser
import numpy as np
from scipy.integrate import odeint
import functools as functool
import general_methods as gm
import time as t

import matplotlib.pyplot as plt

def _get_arguments():
    """Return script arguments as dictionary
    Args:
        - None
    Returns:
        - dictionary containing script arguments
    """
    parser = ArgumentParser()
    parser.add_argument(
        "ExperimentOutputFolder", help="Path to Output Folder"
    )
    parser.add_argument(
        "ExperimentInputFolder", help="Path to Input-Folder whith '.yml' File"
    )
    """ Potentially more arguments """
    return parser.parse_args()


def main():
    args = _get_arguments()
    # Path and file names
    #------------------------------
    PATH_IN = args.ExperimentInputFolder + '/'
    PATH_Out = args.ExperimentOutputFolder + '/'

    FileName_In = 'Parameter.yml'
    FileName_NPZ = 'InfectionData'

    Parameter = gm.yml.read(PATH_IN+FileName_In)
    PatchNr_tot = Parameter['PatchNumber'] #only in one dimension; the actual total Patchnumber is the square of it

    PATH_NPZ = PATH_IN+FileName_NPZ+'_Patch_'
    


    # Control Parameter
    # ----------------------------------------------
    Control = Parameter['Control']
    PatchNr = Parameter['PatchNumber']
    tau_min = Control['tau']['min'] 
    tau_max = Control['tau']['max']
    tau_res = Control['tau']['res']
    total_tau_vec = np.arange(tau_min, tau_max, tau_res)
    n_tau = len(total_tau_vec)

    M_max_min = Control['M_max']['min'] 
    M_max_max = Control['M_max']['max']
    M_max_res = Control['M_max']['res']
    total_M_max_vec = np.arange(M_max_min, M_max_max, M_max_res)
    n_M_max = len(total_M_max_vec)


    # Definition of complete matrices
    # ---------------------------------------------
    LLE = np.zeros((n_M_max, n_tau))
    Winding = np.zeros((n_M_max, n_tau))
    Imean = np.zeros((n_M_max, n_tau))
    
    AllCorrect=True

    for PatchNr in range(PatchNr_tot**2):
        
        n_xpatch = PatchNr%PatchNr_tot # Number of patches in x (tau) direction: itterator and n_xpatch are zero-based
        n_ypatch = PatchNr//PatchNr_tot # Number of patches in y (M_max) direction: itterator and n_xpatch are zero-based
        try: 
            Data = np.load(PATH_NPZ+str(PatchNr)+'.npz')


            #Im folgenen werden vier Faelle unterschieden da die Randstuecke zum die nicht teilbare menge an indizes enthalten, i.e. laenger sind
            if n_xpatch != PatchNr_tot-1 and n_ypatch != PatchNr_tot-1:
                # Testing if Patch is in corrrect domain
                if np.any(Data['tau_vec'] != total_tau_vec[(n_xpatch)*(n_tau//PatchNr_tot):(n_xpatch+1)*(n_tau//PatchNr_tot)]):
                    raise IndexError('Tau vector for Patch number: {} (x-patch number: {}) does not agree with expectations.'.format(PatchNr, n_xpatch))
                if np.any(Data['M_max_vec'] != total_M_max_vec[(n_ypatch)*(n_M_max//PatchNr_tot):(n_ypatch+1)*(n_M_max//PatchNr_tot)]):
                    raise IndexError('M_max vector for Patch number: {} (y-patch number: {}) does not agree with expectations.'.format(PatchNr, n_ypatch))
                #Saving Results
                LLE[(n_ypatch)*(n_M_max//PatchNr_tot):(n_ypatch+1)*(n_M_max//PatchNr_tot), (n_xpatch)*(n_tau//PatchNr_tot):(n_xpatch+1)*(n_tau//PatchNr_tot)] = Data['LLE']
                Winding[(n_ypatch)*(n_M_max//PatchNr_tot):(n_ypatch+1)*(n_M_max//PatchNr_tot), (n_xpatch)*(n_tau//PatchNr_tot):(n_xpatch+1)*(n_tau//PatchNr_tot)] = Data['Winding']
                Imean[(n_ypatch)*(n_M_max//PatchNr_tot):(n_ypatch+1)*(n_M_max//PatchNr_tot), (n_xpatch)*(n_tau//PatchNr_tot):(n_xpatch+1)*(n_tau//PatchNr_tot)] = Data['Imean']
                
            elif n_xpatch != PatchNr_tot-1 and n_ypatch == PatchNr_tot-1:
                # Testing if Patch is in corrrect domain
                if np.any(Data['tau_vec'] != total_tau_vec[(n_xpatch)*(n_tau//PatchNr_tot):(n_xpatch+1)*(n_tau//PatchNr_tot)]):
                    raise IndexError('Tau vector for Patch number: {} (x-patch number: {}) does not agree with expectations.'.format(PatchNr, n_xpatch))
                if np.any(Data['M_max_vec'] != total_M_max_vec[(n_ypatch)*(n_M_max//PatchNr_tot):]):
                    raise IndexError('M_max vector for Patch number: {} (y-patch number: {}) does not agree with expectations.'.format(PatchNr, n_ypatch))
                #Saving Results
                LLE[(n_ypatch)*(n_M_max//PatchNr_tot):, (n_xpatch)*(n_tau//PatchNr_tot):(n_xpatch+1)*(n_tau//PatchNr_tot)] = Data['LLE']
                Winding[(n_ypatch)*(n_M_max//PatchNr_tot):, (n_xpatch)*(n_tau//PatchNr_tot):(n_xpatch+1)*(n_tau//PatchNr_tot)] = Data['Winding']
                Imean[(n_ypatch)*(n_M_max//PatchNr_tot):, (n_xpatch)*(n_tau//PatchNr_tot):(n_xpatch+1)*(n_tau//PatchNr_tot)] = Data['Imean']

            elif n_xpatch == PatchNr_tot-1 and n_ypatch != PatchNr_tot-1:
                # Testing if Patch is in corrrect domain
                if np.any(Data['tau_vec'] != total_tau_vec[(n_xpatch)*(n_tau//PatchNr_tot):]):
                    raise IndexError('Tau vector for Patch number: {} (x-patch number: {}) does not agree with expectations.'.format(PatchNr, n_xpatch))
                if np.any(Data['M_max_vec'] != total_M_max_vec[(n_ypatch)*(n_M_max//PatchNr_tot):(n_ypatch+1)*(n_M_max//PatchNr_tot)]):
                    raise IndexError('M_max vector for Patch number: {} (y-patch number: {}) does not agree with expectations.'.format(PatchNr, n_ypatch))
                #Saving Results
                LLE[(n_ypatch)*(n_M_max//PatchNr_tot):(n_ypatch+1)*(n_M_max//PatchNr_tot), (n_xpatch)*(n_tau//PatchNr_tot):] = Data['LLE']
                Winding[(n_ypatch)*(n_M_max//PatchNr_tot):(n_ypatch+1)*(n_M_max//PatchNr_tot), (n_xpatch)*(n_tau//PatchNr_tot):] = Data['Winding']
                Imean[(n_ypatch)*(n_M_max//PatchNr_tot):(n_ypatch+1)*(n_M_max//PatchNr_tot), (n_xpatch)*(n_tau//PatchNr_tot):] = Data['Imean']

            elif n_xpatch == PatchNr_tot-1 and n_ypatch == PatchNr_tot-1:
                # Testing if Patch is in corrrect domain
                if np.any(Data['tau_vec'] != total_tau_vec[(n_xpatch)*(n_tau//PatchNr_tot):]):
                    raise IndexError('Tau vector for Patch number: {} (x-patch number: {}) does not agree with expectations.'.format(PatchNr, n_xpatch))
                if np.any(Data['M_max_vec'] != total_M_max_vec[(n_ypatch)*(n_M_max//PatchNr_tot):]):
                    raise IndexError('M_max vector for Patch number: {} (y-patch number: {}) does not agree with expectations.'.format(PatchNr, n_ypatch))
                #Saving Results
                LLE[(n_ypatch)*(n_M_max//PatchNr_tot):, (n_xpatch)*(n_tau//PatchNr_tot):] = Data['LLE']
                Winding[(n_ypatch)*(n_M_max//PatchNr_tot):, (n_xpatch)*(n_tau//PatchNr_tot):] = Data['Winding']
                Imean[(n_ypatch)*(n_M_max//PatchNr_tot):, (n_xpatch)*(n_tau//PatchNr_tot):] = Data['Imean']

            else:
                raise ValueError('No intented loop is taken.')
        except FileNotFoundError:
            AllCorrect=False
            pass
        
        gm.progress_featback.printProgressBar(iteration=PatchNr, total=PatchNr_tot**2, name='Combigning Patches:')

    if AllCorrect==False:
        print('Some Files are not found. Maybe calculations are not completed.')

    np.savez_compressed(PATH_Out+FileName_NPZ, 
                        tau_vec=total_tau_vec,
                        M_max_vec=total_M_max_vec,
                        LLE=LLE,
                        Winding=Winding,
                        Imean=Imean)
    print('Saved Data at: ', PATH_Out+FileName_NPZ+'.npz')

if __name__ =='__main__':
    main()