# -*- coding: utf-8 -*-
# MIT License
# Copyright (c) 2022, Victor Tänzel, Miriam Jäger
# All rights reserved.
"""Classes that store constraint force data as work or force time traces.

The resulting force or work sets are needed for further analysis.
"""

__all__ = ['WorkSet', 'ForceSet', 'save', 'load']

import joblib

import numpy as np

from os.path import basename
from beartype import beartype
from beartype.typing import Optional, Any
from sklearn.base import BaseEstimator, TransformerMixin
from scipy.integrate import cumulative_trapezoid

from dcTMD._typing import (
    Int,
    Float,
    Str,
    ArrayLikeStr,
    Index1DArray,
    Float2DArray,
)


@beartype
def save(
    filename: Str,
    classobject,
) -> None:
    """
    Save a class object: a data handler or an estimator.

    Parameters
    ----------
    filename :
        File name to which classobject is saved.
    classobject :
        Instance of the data handler, i.e. a WorkSet or ForceSet instance, or
        of an estimator, i.e. a WorkEstimator or ForceEstimator instance.

    Examples
    --------
    >>> # Save Estimators and data handlers. Here: WorkSet.
    >>> # Save a WorkSet instance named work_set and load it again:
    >>> from dcTMD.storing import save, load
    >>> save(work_set, 'my_workset.joblib')  # noqa: F821
    >>> my_workset = load('my_workset.joblib')
    """
    joblib.dump(classobject, filename)


@beartype
def load(
    filename: Str,
) -> Any:
    """
    Load a data handler or an estimator.

    Parameters
    ----------
    filename :
        Name of the file containing the data handler.

    Returns
    -------
    handler:
        Loaded class object.

    Examples
    --------
    >>> # Loads estimators and data handlers. Here: WorkSet.
    >>> # Save a WorkSet instance named work_set and load it again:
    >>> from dcTMD.storing import save, load
    >>> save(work_set, 'my_workset.joblib')  # noqa: F821
    >>> my_workset = load('my_workset.joblib')
    """
    return joblib.load(filename)


@beartype
def _integrate_force(
    forceorworkset,
    force_data,
):
    """Integrate a force time trace and return the work time trace."""
    work_data = cumulative_trapezoid(
        force_data,
        forceorworkset.position_,
        initial=0,
    )
    return work_data[::forceorworkset.resolution]


@beartype
def _get_time_from_testfile(forceorworkset):
    """Read test force file to determine the time trace."""
    forceorworkset.time_ = np.loadtxt(
        forceorworkset.X[0],
        comments=('@', '#'),
        usecols=[0],
    )
    if forceorworkset.verbose:
        print(f'Using {forceorworkset.X[0]} to initialize arrays.')
        time_length = len(forceorworkset.time_)
        lenth_reduced = forceorworkset.time_[::forceorworkset.resolution]
        time_length_reduced = len(lenth_reduced)
        print(f'length of pullf file is {time_length}')
        print(f'reduced length is {time_length_reduced}')


class WorkSet(TransformerMixin, BaseEstimator):
    """
    Class for managing constraint work data.

    Parameters
    ----------
    velocity :
        Pulling velocity in nm/ps.
    resolution :
        Striding to reduce work time trace.
    verbose :
        Enables verbose mode.

    Attributes
    ----------
    work_:
        Constraint work time traces, in kJ/mol.
    names_ :
        Constraint force file names corresponding to work time traces.
    time_ :
        Time trace corresponding to the work, in ps.
    position_ :
        Positions time trace, product of time trace and velocity, in nm.

    Examples
    --------
    >>> # Load some file names listed in 'filenames'
    >>> import numpy as np
    >>> from dcTMD.storing import WorkSet
    >>> work_set = WorkSet(velocity=0.001, resolution=1)
    >>> work_set.fit(filenames)  # noqa: F821
    Loading & integrating force files: 100%|████| X/X [XX:XX<00:00,  X.XX/it]
    WorkSet(velocity=0.001)
    >>> work_set.work_.shape
    (N_trajectories_, len(time_))

    >>> # Reduce work by selecting some trajectories via their indices,
    >>> # for example the first three, and receive a new WorkSet instance
    >>> indices = np.array([0, 1, 2])
    >>> reduced_set = work_set.reduce(indices)
    >>> reduced_set.work_.shape
    (3, len(time_))
    """

    @beartype
    def __init__(
        self,
        velocity: Float,
        resolution: Int = 1,
        verbose: bool = False,
    ) -> None:
        """Initialize WorkSet class."""
        self.velocity = velocity
        self.resolution = resolution
        self.verbose = verbose

    @beartype
    def fit(
        self,
        X: ArrayLikeStr,  # noqa: WPS111 N803
        y: Optional[np.ndarray] = None,  # noqa: WPS111
    ):
        """
        Load constraint force files and calculate work time traces.

        Parameters
        ----------
        X :
            File names of constraint force files to be read in and integrated.
        y :
            Not used, present for scikit API consistency by convention.

        Returns
        -------
        self:
            Fitted estimator.
        """
        self.X = X  # noqa: WPS111 N803
        # read a test file for the time trace
        _get_time_from_testfile(self)
        # fill arrays with data
        self._fill_work()
        return self

    @beartype
    def transform(self, X, y=None) -> Float2DArray:  # noqa: WPS111 N803
        """Return work set."""
        return self.work_

    @beartype
    def reduce(self, indices: Index1DArray):
        """
        Reduce work set to a chosen subset and return new instance.

        Parameters
        ----------
        indices :
            Indices corresponding to the work trajectories that are kept in the
            work set.

        Returns
        -------
        self :
            Instance of WorkSet.
        """
        import copy
        reduced_work_set = copy.deepcopy(self)
        reduced_work_set.work_ = reduced_work_set.work_[indices]
        reduced_work_set.names_ = reduced_work_set.names_[indices]
        return reduced_work_set

    @beartype
    def _fill_work(self) -> None:
        """Help integrate the force files."""
        import tqdm
        self.work_ = np.zeros(
            (len(self.X), len(self.time_[::self.resolution])),
            dtype=float,
        )
        self.names_ = np.array([])
        self.position_ = self.time_ * self.velocity
        # read in data and fill work_array
        pbar = tqdm.tqdm
        for idx, file_name in pbar(
            enumerate(self.X),
            total=len(self.X),
            desc='Loading & integrating force files',
        ):
            if self.verbose:
                pbar.write(f'Reading file {file_name}')
            file_data = np.loadtxt(
                file_name,
                comments=('@', '#'),
                usecols=[1],
            )
            # test if file is not corrupted, else add it
            if file_data.shape == self.time_.shape:
                self.work_[idx, :] = _integrate_force(self, file_data)
                short_name = basename(file_name)
                self.names_ = np.append(self.names_, short_name)
            else:
                pbar.write(f'skip file {file_name}')
                pbar.write(
                    f'shape is {file_data.shape} and not {self.time_.shape}'
                )
        # removing rows with only zero, reduce positions resolution
        self.work_ = self.work_[~np.all(self.work_ == 0, axis=1)]
        self.position_ = self.position_[::self.resolution]


class ForceSet(TransformerMixin, BaseEstimator):
    """
    Class for managing constraint force data.

    Parameters
    ----------
    velocity :
        Pulling velocity in nm/ps.
    resolution :
        Striding to reduce work time trace. This parameter
        is only added for compatibility with WorkSet
    verbose :
        Enables verbose mode.

    Attributes
    ----------
    force_:
        Constraint force time traces, in kJ/mol.
    names_ :
        Constraint force file names corresponding to force time traces.
    time_ :
        Time trace corresponding to the force, in ps.
    position_ :
        Positions time trace, product of time trace and velocity, in nm.

    Examples
    --------
    >>> # Load some file names listed in 'filenames'
    >>> import numpy as np  # noqa: F401
    >>> from dcTMD.storing import ForceSet
    >>> force_set = ForceSet(velocity=0.001, resolution=1)
    >>> force_set.fit(filenames)
    Loading force files: 100%|████| X/X [XX:XX<00:00,  X.XX/it]
    ForceSet(velocity=0.001)
    >>> force_set.work_.shape
    (N_trajectories_, len(time_))
    """

    @beartype
    def __init__(
        self,
        velocity: Float,
        resolution: Int = 1,
        verbose: bool = False,
    ) -> None:
        """Initialize WorkSet class."""
        self.velocity = velocity
        self.resolution = resolution
        self.verbose = verbose

    @beartype
    def fit(
        self,
        X: ArrayLikeStr,  # noqa: WPS111 N803
        y: Optional[np.ndarray] = None,  # noqa: WPS111
    ):
        """
        Load constraint force files.

        Parameters
        ----------
        X :
            File names of constraint force files to be read in.
        y :
            Not used, present for scikit API consistency by convention.

        Returns
        -------
        self:
            Fitted estimator.
        """
        self.X = X  # noqa: WPS111 N803
        # read a test file for the time trace
        _get_time_from_testfile(self)
        # fill arrays with data
        self._fill_force()
        self.integrate()
        return self

    @beartype
    def transform(self, X, y=None) -> Float2DArray:  # noqa: WPS111 N803
        """Return force set."""
        return self.force_

    @beartype
    def integrate(self) -> None:
        """Integrate forces and return a WorkSet instance."""
        # (1) Instantiate a WorkSet with velocity, resolution and verbose
        # (2) Save names_, time_, position_ attributes manually
        # (3) Integrate the forces in force_ with _integrate_force
        # (4) Return WorkSet instance.
        # Be careful with the resolution so as not to reduce it twice.
        print('integrating forceset --> workset')
        self.work_ = _integrate_force(self, self.force_)[::self.resolution]
        print(self.work_.shape)

    @beartype
    def _fill_force(self) -> None:
        """Help load the force files."""
        # Load force files
        # Check if files are corrupt and build names_, position_ and work_
        import tqdm
        self.force_ = np.zeros(
            (len(self.X), len(self.time_)),
            dtype=float,
        )
        self.names_ = np.array([])
        self.position_ = self.time_ * self.velocity
        # read in data and fill force_array
        pbar = tqdm.tqdm
        for idx, file_name in pbar(
            enumerate(self.X),
            total=len(self.X),
            desc='Loading force files',
        ):
            if self.verbose:
                pbar.write(f'Reading file {file_name}')
            file_data = np.loadtxt(
                file_name,
                comments=('@', '#'),
                usecols=[1],
            )
            # test if file is not corrupted, else add it
            if file_data.shape == self.time_.shape:
                self.force_[idx, :] = file_data
                short_name = basename(file_name)
                self.names_ = np.append(self.names_, short_name)
            else:
                pbar.write(f'skip file {file_name}')
                pbar.write(
                    f'shape is {file_data.shape} and not {self.time_.shape}'
                )
        # removing rows with only zero
        self.force_ = self.force_[~np.all(self.force_ == 0, axis=1)]
