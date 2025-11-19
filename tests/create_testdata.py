"""Create test data."""

import dcTMD
from dcTMD.storing import save, WorkSet, ForceSet
from dcTMD.dcTMD import WorkEstimator, ForceEstimator
from dcTMD.io import write_output
import numpy as np

# define variables
velocity = 0.001
res = 1
res10 = 10
verbose = True
temperature = 300
n_resamples = 100
seed = 42
sigma = 0.1
mode = 'nearest'
bootstrapmode = 'std'

pullf_files = 'testdata/pullf_filenames.dat'
filenames = dcTMD.io.load_pullf(pullf_files)


# define three vaid indices for memory kernel testing
def _valid_indices(forceestimator):
    n = len(forceestimator.force_set.time_)
    # Choose three well-spaced, guaranteed-valid indices.
    # Ensure they are >=1 and <= n-2 to avoid boundary weirdness.
    return np.array([1, max(2, n // 3), n - 2], dtype=int)


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
forceestimator.memory_kernel(_valid_indices(forceestimator))
save('testdata/forceestimator', forceestimator)

# create ForceSet instance with res10
forceset = ForceSet(velocity=velocity,
                    resolution=res10,
                    verbose=verbose,
                    )
# fit/fill workset
forceset.fit(filenames)
# save workset
save('testdata/forceset_res10', forceset)
# Instantiate a ForceEstimator instance and fit it with the ForceSet
# instance
forceestimator = ForceEstimator(temperature)
forceestimator.fit(forceset)
# smooth friction
forceestimator.smooth_friction(sigma, mode=mode)
forceestimator.memory_kernel(_valid_indices(forceestimator))
save('testdata/forceestimator_res10', forceestimator)

# create WorkSet instance
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
write_output(
    'testdata/workeestimator_dcTMDresults',
    workeestimator,
    filetype=('dat', 'npz')
)

# create WorkSet instance with res10
workset = WorkSet(velocity=velocity,
                  resolution=res10,
                  verbose=verbose,
                  )
# fit/fill workset
workset.fit(filenames)
# save workset
save('testdata/workset_res10', workset)
# Instantiate a ForceEstimator instance and fit it with the ForceSet
# instance
workeestimator = WorkEstimator(temperature)
workeestimator.fit(workset)
# smooth friction
workeestimator.smooth_friction(sigma, mode=mode)
# error estimation vis bootstrapping
workeestimator.estimate_free_energy_errors(n_resamples, bootstrapmode, seed)
save('testdata/workeestimator_res10', workeestimator)
