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


def test_main(tmpdir):
    """Test main CLI."""
    # Create temporary directory for output files
    output = tmpdir.join('test')
    # Set up test command line arguments
    args = [
        '--mode', 'work',
        '--file', f'{TEST_FILE_DIR}/*.xvg',
        '--outname', f'{output}',
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
