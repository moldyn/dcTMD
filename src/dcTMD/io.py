# -*- coding: utf-8 -*-
"""
"""
__all__ = []

from cgi import test
import numpy as np


def load_pullf(pullf_glob_pattern, pullf_files):
    if pullf_glob_pattern!=None:
        import glob
        files = glob.glob(pullf_glob_pattern)
    if pullf_files!=None:
        files = np.loadtxt(pullf_files, dtype=str)  # TODO reshape input?
    return files


def write_dG(o, N, x, av_forceintegral, wdiss):
    print('writing dG...')
    headers = ['#x',
               'force_integral',
               'wdiss',
               'corrected_force_integral',
               ]
    np.savetxt(o + str(N) + '_dG.dat',
               np.c_[x,
                     av_forceintegral,
                     wdiss,
                     av_forceintegral - wdiss,
                     ],
               fmt='%20.8f',
               header='    '.join(headers),
              )


"""
def write_dG(o, N, x, av_forceintegral, av_intcorr, wdiss):
    print('writing dG...')
    headers = ['#x',
               'force_integral',
               'frict_coeff',
               'wdiss',
               'corrected_force_integral',
               ]
    np.savetxt(o + str(N) + '_dG.dat',
               np.c_[x,
                     av_forceintegral,
                     av_intcorr,
                     wdiss,
                     av_forceintegral - wdiss,
                     ],
               fmt=['%15.8f', '%20.8f', '%20.8f', '%20.8f', '%20.8f'],
               header='    '.join(headers),
              )
"""

def write_output(o: str, N: int, **kwargs):
    print('writing dG...')
    results = np.asarray([(name, a) for name, a in kwargs.items() \
                            if name != 'errors'], dtype=object)
    header = results.T[0]
    arrays = np.vstack(results.T[1])
    if kwargs["errors"]:
        #header_err = ['s_' + name for name in header]
        #arrays_err = kwargs['errors']
        header = np.append(header, ['s_' + name for name in header if name != 'x'])
        arrays = np.vstack([arrays, kwargs['errors']])
    #header = [name for name in kwargs.keys() if name != 'errors']
    #arrays = np.asarray([a for a in kwargs.values()])

    np.savetxt(o + "_" + str(N) + '_dG.dat',
               arrays.T,
               fmt='%20.8f',
               header='    '.join(header),
              )
    

#def write_friction(o, N, x, autocorr_set, av_intcorr, blurred, runn_av, action,
#                   action_gaussmooth):
#    frictresult = open(o + str(N) + "_sig" + str(blurr) + "_frict.dat", "w")
#    frictresult.write(
#        "#x   ACF   frict_coeff   gauss_filtered_frict_coeff" +
#        "   av_window_frict_coeff   action   action_gausssmooth\n")
#    for i in range(length_data):
#        frictresult.write(
#            "{:15.8f} {:20.8f} {:20.8f} {:20.8f} {:20.8f} {:20.8f} {:20.8f}\n".format(
#                x[i], autocorr_set[i], av_intcorr[i], blurred[i],
#                runn_av[i], action[i], action_gaussmooth[i]))
#    frictresult.close()
    
    
"""
def write_friction(o, N, x, autocorr_set, av_intcorr, blurred, runn_av, action,
                   action_gaussmooth):
    print('Writing friction...')
    headers = ['#x',
               'ACF',
               'frict_coeff',
               'gauss_filtered_frict_coeff',
               'av_window_frict_coeff',
               'action',
               'action_gausssmooth',
               ]
    np.savetxt(o + str(N) + "_sig" + str(blurr) + "_frict.dat",
               np.c_[x,
                     autocorr_set,
                     av_intcorr,
                     blurred,
                     runn_av,
                     action,
                     action_gaussmooth,
                     ],
               fmt=['%15.8f', '%20.8f', '%20.8f', '%20.8f', '%20.8f', '%20.8f', 
                    '%20.8f'],
               header='    '.join(headers),
               )
"""