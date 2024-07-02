"""Create test data."""

import dcTMD
from dcTMD.storing import save, WorkSet, ForceSet
from dcTMD.dcTMD import WorkEstimator, ForceEstimator


# define variables
velocity = 0.001
res = 1
verbose = True
temperature = 300
n_resamples = 100
seed = 42
sigma = 0.1
mode = 'nearest'
bootstrapmode = 'std'

pullf_files = 'testdata/pullf_filenames.dat'
filenames = dcTMD.io.load_pullf(pullf_files)


# create ForceSet instance
forceset = ForceSet(velocity=velocity,
                    resolution=res,
                    verbose=verbose,
                    )
# fit/fill workset
forceset.fit(filenames)
# save workset
save('testdata/forceset', forceset)
# Instantiate a ForceEstimator instance and fit it with the ForceSet
# instance
forceestimator = ForceEstimator(temperature)
forceestimator.fit(forceset)
# smooth friction
forceestimator.smooth_friction(sigma, mode=mode)
save('testdata/forceestimator', forceestimator)

# create ForceSet instance
workset = WorkSet(velocity=velocity,
                  resolution=res,
                  verbose=verbose,
                  )
# fit/fill workset
workset.fit(filenames)
# save workset
save('testdata/workset', workset)
# Instantiate a ForceEstimator instance and fit it with the ForceSet
# instance
workeestimator = WorkEstimator(temperature)
workeestimator.fit(workset)
# smooth friction
workeestimator.smooth_friction(sigma, mode=mode)
# error estimation vis bootstrapping
workeestimator.estimate_free_energy_errors(n_resamples, bootstrapmode, seed)
save('testdata/workeestimator', workeestimator)
