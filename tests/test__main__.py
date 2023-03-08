# -*- coding: utf-8 -*-
"""Tests for the cli script.
MIT License
Copyright (c) 2021-2022, Miriam Jäger
All rights reserved.
"""
from click.testing import CliRunner

from dcTMD.__main__ import main


def test_main():
    runner = CliRunner()
    result = runner.invoke(main)
    assert result.exit_code == 0
    assert '--help' in result.output
