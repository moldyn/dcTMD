# -*- coding: utf-8 -*-
"""Tests for the cli script.
MIT License
Copyright (c) 2021-2022, Miriam Jäger
All rights reserved.
"""
from click.testing import CliRunner
from dcTMD.__main__ import main


def test_main(tmpdir):
    # Create temporary directory for output files
    output = tmpdir.join('test')
    # Set up test command line arguments
    args = [
        '--mode', 'work',
        '--file', 'testdata/*.xvg',
        '--outname', output,
        '--temperature', '300',
        '--velocity', '0.01',
        '--res', '1',
        '--sigma', '0.1',
        '--verbose',
        '--plot',
        '--save_dataset'
    ]
    runner = CliRunner()
    result = runner.invoke(main, args)
    assert result.exit_code == 0
