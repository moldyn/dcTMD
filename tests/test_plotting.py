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
from dcTMD.dcTMD import WorkEstimator
from dcTMD.storing import load
from dcTMD.utils.plotting import (
    plot_dcTMD_results,
    plot_dG_Wdiss,
    plot_Gamma,
    plot_dG,
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
    workset_name = f'{TEST_FILE_DIR}/workset'
    workset = load(filename=workset_name)
    estimator = WorkEstimator(temperature=TEMPERATURE)
    return estimator.fit(workset)


def test_plot_dcTMD_results(ref_workestimator):
    x = ref_workestimator.position_
    friction = ref_workestimator.friction_
    fig, axs = plot_dcTMD_results(x, ref_workestimator, friction)
    assert isinstance(fig, plt.Figure)
    assert len(axs) == 2
    assert len(axs[0].lines) == 3
    assert len(axs[1].lines) == 1
    for ax in axs:
        assert isinstance(ax, plt.Axes)


def test_plot_dG_Wdiss(ref_workestimator):
    x = ref_workestimator.position_
    fig, ax = plt.subplots()
    plot_dG_Wdiss(x, ref_workestimator, ax)
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
    assert ax.get_ylabel() == "$\\Gamma$ [kJ/mol/(nm$^2$/ps)]"


def test_plot_dG():
    fig, ax = plt.subplots()
    x = np.linspace(0, 10, 100)
    dG = np.sin(x)
    plot_dG(x, dG, ax, label='test')
    assert ax.get_xlabel() == r'position $x$ [nm]'
    assert ax.get_ylabel() == r'$\Delta G$ [kJ/mol]'
    assert ax.get_xlim() == (0, 10)


def test_plot_worklines():
    fig, ax = plt.subplots()
    x = np.linspace(0, 10, 100)
    workset = np.random.rand(10, 100)
    plot_worklines(x, workset, ax)
    assert ax.get_xlabel() == r'position $x$ [nm]'
    assert ax.get_ylabel() == r'work $W$ [kJ/mol]'
    assert ax.get_xlim() == (0, 10)


def test_plot_histo_normaldist():
    fig, ax = plt.subplots()
    data = np.random.normal(0, 1, 100)
    plot_histo_normaldist(data, ax, color='blue')
    assert ax.get_xlabel() == r'$P$'
    assert ax.get_ylabel() == r'$W$ [kJ/mol]'
    assert len(ax.get_lines()) == 1


def test_plot_worknormalitychecks():
    x = np.linspace(0, 1, 10)
    workset = np.random.rand(5, 10)
    index = [2, 5, 8]
    fig, axs = plot_worknormalitychecks(x, workset, index, colors=None)
    assert axs[0].get_xlabel() == 'position $x$ [nm]'
    assert axs[0].get_ylabel() == r'work $W$ [kJ/mol]'
    assert axs[1].get_xlabel() == r'$P$'
    assert axs[1].get_ylabel() == r'$W$ [kJ/mol]'
