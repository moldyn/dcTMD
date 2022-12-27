# -*- coding: utf-8 -*-
"""
Classes storing constraint force data as work or force time traces.


MIT License
Copyright (c) 2022, Victor Tänzel, Miriam Jäger
All rights reserved.
"""

__all__ = ['WorkSet', 'ForceSet']

import numpy as np
from beartype import beartype
from beartype.typing import Optional
from sklearn.base import BaseEstimator, TransformerMixin

from dcTMD._typing import (
    Int,
    Float,
    ArrayLikeStr,
    Index1DArray,
    Float2DArray,
)


@beartype
def _get_time_from_testfile(manager):
    """Read test force file to determine the time trace for a data manager."""
    if manager.verbose:
        print(f'Using {manager.X[0]} to initialize arrays.')
        time_length = len(manager.time)
        time_length_reduced = len(manager.time[::manager.resolution])
        print(f'length of pullf file is {time_length}')
        print(f'reduced length is {time_length_reduced}')
    manager.time_ = np.loadtxt(
        manager.X[0],
        comments=('@', '#'),
        usecols=[0],
    )


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
    # Load some file names listed in 'filenames'
    >>> import numpy as np
    >>> from dcTMD.storing import WorkSet
    >>> work_set = WorkSet(velocity=0.001, resolution=1)
    >>> work_set.fit(filenames)
    Loading & integrating force files: 100%|████| X/X [XX:XX<00:00,  X.XX/it]
    WorkSet(velocity=0.001)
    >>> work_set.work_.shape
    (N_trajectories_, len(time_))

    # Reduce work by selecting some trajectories via their indices,
    # for example the first three, and receive a new WorkSet instance
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
    def fit(self, X: ArrayLikeStr, y: Optional[np.ndarray] = None):
        """
        Loads constraint force files and calculates work time traces.

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
        self.X = X
        # read a test file for the time trace
        _get_time_from_testfile(self)
        # fill arrays with data
        self._fill_work()
        return self

    @beartype
    def transform(self, X, y=None) -> Float2DArray:
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
        from scipy.integrate import cumulative_trapezoid
        self.work_ = np.zeros(
            (len(self.X), len(self.time_[::self.resolution])),
            dtype=float,
        )
        self.names_ = np.array([])
        self.position_ = self.time_ * self.velocity
        # read in data and fill work_array
        for idx, file_name in (pbar := tqdm.tqdm)(
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
            # test if file is corrupted, else add it
            if file_data.shape == self.time_.shape:
                self.work_[idx, :] = cumulative_trapezoid(
                    file_data,
                    self.position_,
                    initial=0,
                )[::self.resolution]
                self.names_ = np.append(self.names_, file_name)
            else:
                pbar.write(f'skip file {file_name}')
                pbar.write(
                    f'shape is {file_data.shape} and not {self.time_.shape}')
        # removing rows with only zero, reduce positions resolution
        self.work_ = self.work_[~np.all(self.work_ == 0, axis=1)]
        self.position_ = self.position_[::self.resolution]
