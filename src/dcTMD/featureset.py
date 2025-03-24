import glob
import re
import numpy as np
from tqdm import tqdm
from beartype import beartype
from beartype.typing import List, Union


class FeatureSet:
    """
    A class to create a array based on trajectory features.
    This array is later on used to calculate similarities between trajectories.

    The class loads trajectory-specific features and organizes them into arrays for further analysis.
    It assumes that all files have the same shape. The order in which they are loaded needs to be specified.
    """

    @beartype
    def __init__(
        self,
        filenames: Union[List[str], np.ndarray] = None,
        filenameprefix: Union[List[str], np.ndarray] = None,
        wildcard: str = '*{}*',
        verbose: bool = False
    ) -> None:
        """
        Initialize the FeatureSet.

        Parameters
        ----------
        filenames :
            A list or array of filenames. If provided, it will be used directly.
        filenameprefix :
            A list or array of names that contain the running number.
            Must be used together with wildcard.
        wildcard :
            A wildcard string pattern for generating filenames
            in combination with the running number from filenameprefix.
            Must be provided if filenameprefix is used, by default '*{}*'.
        verbose : bool, optional
            Whether to print logs during execution, by default False.

        Examples
        --------
        >>> filenames = ['traj_01_contacts.txt', 'traj_02_contacts.txt']
        >>> featureset = FeatureSet(filenames, verbose=True)
        >>> featureset.fill_array()
        >>> array = featureset.array
        >>> 100%|████████████| X/X [XX:XX<00:00, X.XXit/s]
        """
        self.verbose = verbose

        if filenames is not None:
            self.filenames = np.asarray(filenames)
        elif filenameprefix is not None and wildcard:
            self.filenames = self._get_filenames(filenameprefix, wildcard)
        else:
            raise ValueError(
                "Either `filenames` must be provided directly, "
                "or both `filenameprefix` and `wildcard` must be provided together."
            )
        if self.verbose:
            print(f"Loaded filenames: {self.filenames}")

    def _get_filenames(
        self,
        filename_prefix: Union[List[str], np.ndarray],
        filename_wildcard: str
    ) -> np.ndarray:
        """
        Generate filenames based on the provided wildcard.

        This method searches for the running number in the filename prefixes and
        inserts it into the wildcard pattern.

        Parameters
        ----------
        filename_prefix :
            A list or array of filename prefixes that contain the running number/index.
        filename_wildcard :
            A wildcard pattern (e.g., "wildcard_{}.xvg") for filenames.

        Returns
        -------
        np.ndarray
            An array of filenames matching the wildcard.

        Examples
        --------
        >>> filenames = ['prefix_01', 'prefix_02']
        >>> wc = 'wildcard_{}.xvg'
        >>> FeatureSet(filenames, wc)._get_filenames(filenames, wc)
        array(['wildcard_01.xvg', 'wildcard_02.xvg'])
        """
        filenames = []
        rn_pattern = r'_(\d{2,4})'

        for string in filename_prefix:
            rn_match = re.search(rn_pattern, string)
            running_number = rn_match.group(1) if rn_match else None
            if running_number:
                name = filename_wildcard.format(running_number)
                matched_name = glob.glob(name)
                if not matched_name:
                    raise FileNotFoundError(f'No file found matching pattern: {name}.')
                filenames.extend(matched_name)

        self.filenames = np.asarray(filenames).flatten()
        if len(filenames) != len(filename_prefix):
            raise ValueError(f"Number of files: {len(filenames)} does not match number of prefixes: {len(filename_prefix)}.")
        return self.filenames

    def _read_testfile(self) -> None:
        """
        Read the first file in the list of filenames to initialize the array shape.

        """
        testfile = np.loadtxt(self.filenames[0], comments=("@", "#"))
        self.fileshape = testfile.shape
        if self.verbose:
            print(f"Using {self.filenames[0]} to initialize arrays with shape {self.fileshape}")

    def fill_array(self) -> np.ndarray:
        """
        Load the data from trajectory files into a NumPy array.

        This method reads each file in self.filenames and fills the data into a pre-allocated NumPy array based on the file shape determined by `_read_testfile`.

        Files that cannot be loaded or its shape does not match the expected shape, are skipped.

        Returns
        -------
        np.ndarray
            A array where each entry corresponds to the data from a trajectory file.

        Examples
        --------
        >>> filenames = ['traj_01_contacts.txt', 'traj_02_contacts.txt']
        >>> featureset = FeatureSet(filenames, verbose=True)
        >>> featureset.fill_array()
        >>> array = featureset.array
        >>> 100%|████████████| X/X [XX:XX<00:00, X.XXit/s]
        """
        self._read_testfile()
        array = np.zeros(shape=(
            len(self.filenames),
            *self.fileshape,
        ))
        self.featureset_names = []

        with tqdm(
            total=len(self.filenames),
            desc='Loading files',
        ) as pbar:
            for i, current_file_name in enumerate(self.filenames):
                if self.verbose:
                    print(f"Reading file {current_file_name}")
                try:
                    input_data = np.loadtxt(current_file_name, comments=("@", "#"))
                except Exception:
                    raise FileNotFoundError(f"File {current_file_name} not found")

                if input_data.shape != self.fileshape:
                    print(f"Skipping file {current_file_name} due to shape mismatch")
                    continue

                array[i, :] = input_data
                self.featureset_names.append(current_file_name)
                pbar.update(1)

        self.array = array
        return self.array
