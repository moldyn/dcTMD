# -*- coding: utf-8 -*-
"""
Class for defining typing.

In this module all type hint types are defined. Adapted copy of the file with
the same name from MoSAIC (github.com/moldyn/MoSAIC).

MIT License
Copyright (c) 2021-2022, Daniel Nagel, Victor Tanzel
All rights reserved.
"""

import numpy as np
from beartype.typing import List, Union
from beartype.vale import Is, IsAttr, IsEqual

try:  # for python <= 3.8 use typing_extensions
    from beartype.typing import Annotated
except ImportError:
    from typing_extensions import Annotated


class NDim:
    """Class for creating Validators checking for desired dimensions."""

    def __class_getitem__(self, ndim):
        return IsAttr['ndim', IsEqual[ndim]]


class DType:
    """Class for creating Validators checking for desired dtype."""

    def __class_getitem__(self, dtype):
        return Is[lambda arr: np.issubdtype(arr.dtype, dtype)]


# Define Validators
# cast np.bool_ return type of np.all to bool to avoid tri-state boolean
# error, see beartype #153
IsDTypeLike = Is[lambda dtype: np.issubdtype(dtype, np.generic)]
IsPositive = Is[lambda arr: bool(np.all(arr >= 0))]
IsStd = Is[lambda val: val in ['std']]
IsLessThanOne = Is[lambda arr: bool(np.all(arr <= 1))]
IsStrictlyPositive = Is[lambda arr: bool(np.all(arr > 0))]

# Define Types
# # scalar datatypes
Int = Union[int, np.integer]
Float = Union[float, np.floating]
Str = Union[str, np.str_]
NumInRange0to1 = Annotated[Union[Int, Float], IsPositive & IsLessThanOne]

# beartype substitute for np.typing.DTypeLike
DTypeLike = Annotated[type, IsDTypeLike]

# array datatypes
IntNDArray = Annotated[np.ndarray, DType[np.integer]]
FloatNDArray = Annotated[np.ndarray, DType[np.floating]]
StrNDArray = Annotated[np.ndarray, DType[np.str_]]
StrStd = Annotated[Str, IsStd]

Index1DArray = Annotated[IntNDArray, NDim[1] & IsPositive]
Float1DArray = Annotated[FloatNDArray, NDim[1]]
Float2DArray = Annotated[FloatNDArray, NDim[2]]
Float3DArray = Annotated[FloatNDArray, NDim[3]]
Str1DArray = Annotated[StrNDArray, NDim[1]]

ArrayLikeFloat = Union[List[float], FloatNDArray]
ArrayLikeStr = Union[List[str], Str1DArray]
