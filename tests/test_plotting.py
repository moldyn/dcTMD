# -*- coding: utf-8 -*-
"""
Tests for the plotting module.

MIT License
Copyright (c) 2021-2022, Miriam Jäger
All rights reserved.
"""

import pytest
import numpy as np
import matplotlib.pyplot as plt
from os.path import dirname, join
from dcTMD.storing import load
from dcTMD.utils.plotting import (
    plot_dcTMD_results,
    plot_dG_Wdiss,
    plot_Gamma,
    plot_dG,
    plot_dG_werrors,
    plot_worklines,
    plot_histo_normaldist,
    plot_worknormalitychecks,
)


VELOCITY = 0.001
RESOLUTION = 1
VERBOSE = True
TEMPERATURE = 300
INDICES = np.array([1, 3, 5])
HERE = dirname(__file__)
TEST_FILE_DIR = join(HERE, 'testdata')


@pytest.fixture
def ref_workestimator(scope="session"):
    workeestimator_name = f'{TEST_FILE_DIR}/workeestimator'
    return load(filename=workeestimator_name)


@pytest.fixture
def ref_workset(scope="session"):
    workset_name = f'{TEST_FILE_DIR}/workset'
    return load(filename=workset_name)


def test_plot_dcTMD_results(ref_workestimator):
    x = ref_workestimator.position_
    friction = ref_workestimator.friction_
    fig, axs = plot_dcTMD_results(ref_workestimator, friction, x)
    assert isinstance(fig, plt.Figure)
    assert len(axs) == 2
    assert len(axs[0].lines) == 3
    assert len(axs[1].lines) == 1
    for ax in axs:
        assert isinstance(ax, plt.Axes)


def test_plot_dG_Wdiss(ref_workestimator):
    x = ref_workestimator.position_
    fig, ax = plt.subplots()
    plot_dG_Wdiss(ref_workestimator, ax, x)
    assert len(ax.lines) == 3
    assert ax.get_xlabel() == "position $x$ [nm]"
    assert ax.get_ylabel() == "[kJ/mol]"


def test_plot_Gamma():
    x = np.array([0, 1, 2, 3, 4])
    friction = np.array([0, 1, 2, 3, 4])
    fig, ax = plt.subplots()
    plot_Gamma(x, friction, ax)
    assert len(ax.lines) == 1
    assert ax.get_xlabel() == "position $x$ [nm]"
    assert ax.get_xlim() == (0, 4)


def test_plot_dG():
    fig, ax = plt.subplots()
    x = np.linspace(0, 10, 100)
    dG = np.sin(x)
    plot_dG(x, dG, ax, label='test')
    assert ax.get_xlabel() == r'position $x$ [nm]'
    assert ax.get_ylabel() == r'$\Delta G$ [kJ/mol]'
    assert ax.get_xlim() == (0, 10)


def test_plot_dG_werrors(ref_workestimator):
    fig, ax = plt.subplots()
    assert hasattr(ref_workestimator, 's_dG_')  # noqa: WPS421
    plot_dG_werrors(ref_workestimator, ax)
    # one for dG and one for the error bounds
    assert len(ax.lines) >= 1
    assert ax.get_ylabel() == r'$\Delta G$ [kJ/mol]'
    assert ax.get_xlabel() == r'position $x$ [nm]'


def test_plot_dG_werrors_no_errors(ref_workestimator, capsys):
    fig, ax = plt.subplots()
    # Remove the s_dG_ attribute
    if hasattr(ref_workestimator, 's_dG_'):
        del ref_workestimator.s_dG_
    plot_dG_werrors(ref_workestimator, ax)
    captured = capsys.readouterr()
    assert "no errors are determined" in captured.out
    # No lines should be plotted
    assert len(ax.lines) == 0


def test_plot_worklines(ref_workset):
    fig, ax = plt.subplots()
    plot_worklines(ref_workset, ax)
    assert ax.get_xlabel() == r'position $x$ [nm]'
    assert ax.get_ylabel() == r'work $W$ [kJ/mol]'


def test_plot_histo_normaldist():
    fig, ax = plt.subplots()
    data = np.random.normal(0, 1, 100)
    plot_histo_normaldist(data, ax, color='blue')
    assert ax.get_xlabel() == r'$P$'
    assert ax.get_ylabel() == r'$W$ [kJ/mol]'
    assert len(ax.get_lines()) == 1


def test_plot_worknormalitychecks(ref_workset):
    index = [2, 5, 8]
    fig, axs = plot_worknormalitychecks(
        ref_workset,
        index,
    )
    assert axs[0].get_xlabel() == 'position $x$ [nm]'
    assert axs[0].get_ylabel() == r'work $W$ [kJ/mol]'
    assert axs[1].get_xlabel() == r'$P$'
    assert axs[1].get_ylabel() == r'$W$ [kJ/mol]'
