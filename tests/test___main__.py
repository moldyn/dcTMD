# -*- coding: utf-8 -*-
"""Tests for the cli script.
MIT License
Copyright (c) 2021-2022, Miriam Jäger
All rights reserved.
"""

from click.testing import CliRunner
from dcTMD.__main__ import main
from os.path import dirname, join

HERE = dirname(__file__)
TEST_FILE_DIR = join(HERE, 'testdata')


def test_main_work_mode(tmpdir):
    """Test main CLI in work mode."""
    # Create temporary directory for output files
    output = tmpdir.join('test_work')
    # Set up test command line arguments
    args = [
        '--mode', 'work',
        '--file', f'{TEST_FILE_DIR}/*.xvg',
        '--outname', str(output),
        '--temperature', '300',
        '--velocity', '0.01',
        '--res', '1',
        '--sigma', '0.1',
        '--verbose',
        '--plot',
        '--save_dataset',
    ]
    runner = CliRunner()
    clirunner_result = runner.invoke(main, args)
    assert clirunner_result.exit_code == 0


def test_main_force_mode(tmpdir):
    """Test main CLI in force mode."""
    # Create temporary directory for output files
    output = tmpdir.join('test_force')
    # Set up test command line arguments
    args = [
        '--mode', 'force',
        '--file', f'{TEST_FILE_DIR}/*.xvg',
        '--outname', str(output),
        '--temperature', '300',
        '--velocity', '0.01',
        '--res', '1',
        '--sigma', '0.1',
        '--verbose',
        '--plot',
        '--save_dataset',
    ]
    runner = CliRunner()
    clirunner_result = runner.invoke(main, args)
    assert clirunner_result.exit_code == 0


def test_main_no_plot(tmpdir):
    """Test main CLI without plot option."""
    # Create temporary directory for output files
    output = tmpdir.join('test_no_plot')
    # Set up test command line arguments
    args = [
        '--mode', 'work',
        '--file', f'{TEST_FILE_DIR}/*.xvg',
        '--outname', str(output),
        '--temperature', '300',
        '--velocity', '0.01',
        '--res', '1',
        '--sigma', '0.1',
        '--verbose',
        '--save_dataset',
    ]
    runner = CliRunner()
    clirunner_result = runner.invoke(main, args)
    assert clirunner_result.exit_code == 0


def test_main_no_save_dataset(tmpdir):
    """Test main CLI without save_dataset option."""
    # Create temporary directory for output files
    output = tmpdir.join('test_no_save')
    # Set up test command line arguments
    args = [
        '--mode', 'work',
        '--file', f'{TEST_FILE_DIR}/*.xvg',
        '--outname', str(output),
        '--temperature', '300',
        '--velocity', '0.01',
        '--res', '1',
        '--sigma', '0.1',
        '--verbose',
        '--plot',
    ]
    runner = CliRunner()
    clirunner_result = runner.invoke(main, args)
    assert clirunner_result.exit_code == 0


def test_main_minimal(tmpdir):
    """Test main CLI with minimal options."""
    # Create temporary directory for output files
    output = tmpdir.join('test_minimal')
    # Set up test command line arguments
    args = [
        '--file', f'{TEST_FILE_DIR}/*.xvg',
        '--outname', str(output),
        '--temperature', '300',
        '--velocity', '0.01',
    ]
    runner = CliRunner()
    clirunner_result = runner.invoke(main, args)
    assert clirunner_result.exit_code == 0