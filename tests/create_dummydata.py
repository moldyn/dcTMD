import numpy as np
import re

rng = np.random.default_rng(seed=100)

pullf_filenames = np.loadtxt('testdata/pullf_filenames.dat', dtype=str)


feature_array_2D = np.zeros((len(pullf_filenames), 100))
feature_array_3D = np.zeros((len(pullf_filenames), 100, 10))

names1D = []
names2D = []
for j, name_string in enumerate(pullf_filenames):
    # get running number from trajectory file
    rn_pattern = r'_(\d{2,4})'
    rn_match = re.search(rn_pattern, name_string)
    running_number = rn_match.group(1) if rn_match else None
    if j < 13:
        my = 1
    else:
        my = 10

    rnddist1D = np.random.normal(1, 1, 100)
    feature_array_2D[j] = rnddist1D
    np.savetxt(f'testdata/feature_1D_{running_number}.txt', rnddist1D)
    names1D.append(f'feature_1D_{running_number}.txt')
    rnddist2D = np.random.normal(1, 1, (100, 10))
    feature_array_3D[j] = rnddist2D
    np.savetxt(f'testdata/feature_2D_{running_number}.txt', rnddist2D)
    names2D.append(f'feature_2D_{running_number}.txt')

np.save('testdata/feature_array_2D.npy', feature_array_2D)
np.savetxt('testdata/feature_1D_filenamesname.txt', names1D, fmt='%s')
np.save('testdata/feature_array_3D.npy', feature_array_3D)
np.savetxt('testdata/feature_2D_filenamesname.txt', names2D, fmt='%s')
