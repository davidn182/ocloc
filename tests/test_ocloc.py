import os
import sys
import pytest

import obspy


from ocloc.ocloc import read_correlation_file


class TestReadCorrelationFile:

    # Test a failure.
    def test_file_not_found(self):
        fname = 'tests/data/badsac.sac'
        with pytest.raises(FileNotFoundError, match='file does not exist'):
            read_correlation_file(fname)

    # 
    def test_file_not_found(self):
        fname = 'tests/data/KEF_O01_1413547247_100.sac'
        tr = read_correlation_file(fname)

        date = obspy.UTCDateTime(2014, 10, 17, 12, 0, 47)
        assert tr.stats.average_date == date
        assert tr.stats.number_of_days == 100
        assert tr.stats.station_pair == 'KEF_O01'


# pytest -- runs all tests
# pytest tests/test_bla.py  -- runs this file
# pytest tests/test_bla.py::test_a  -- runs this function
# pytest tests/test_bla.py::TestA::test_a  -- runs this function in this class
# flags: -svv
# https://docs.pytest.org
#
# capsys - captures system output
# tmpdir - gives temporary directory
