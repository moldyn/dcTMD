# -*- coding: utf-8 -*-
"""Tests for the cli script.
MIT License
Copyright (c) 2021-2022, Miriam Jäger
All rights reserved.
"""
from click.testing import CliRunner
import tempfile
import os
import shutil
from dcTMD.__main__ import main


def test_main():
    # Create temporary directory for output files
    with tempfile.TemporaryDirectory() as tmpdirname:
        # Set up test command line arguments
        args = [
            '--mode', 'work',
            '--file', 'testdata/*.xvg',
            '--outname', os.path.join(tmpdirname, 'test'),
            '--temperature', '300',
            '--velocity', '0.01',
            '--res', '1',
            '--sigma', '0.1',
            '--resamples', '10',
            '--verbose',
            '--plot',
            '--save_dataset'
        ]
        runner = CliRunner()
        result = runner.invoke(main, args)
        assert result.exit_code == 0

        # Check that output files were created
        # assert os.path.exists(os.path.join(tmpdirname, 'test.npz'))
        # assert os.path.exists(os.path.join(tmpdirname, 'test.dat'))

        # Clean up
        shutil.rmtree(tmpdirname)
