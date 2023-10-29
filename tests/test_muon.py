import sys
from unittest.mock import MagicMock
sys.modules["drift_chamber.widget"] = MagicMock()

import numpy as np
import pytest
from pytest_mock import MockerFixture

from drift_chamber import Muon


def test_init(mocker):
    mock_log_normal = mocker.patch.object(np.random, "lognormal")
    mock_cos = mocker.patch.object(np, "cos")
    mock_cos.return_value = 1
    mock_random = mocker.patch.object(np.random, "random")
    zen = 1
    r = 0
    azi = 1
    x_coord = 1
    y_coord = 1
    charge_random_var = 0
    mock_random.side_effect = [
        charge_random_var,
        (zen - 0.5 * np.pi) / np.pi,
        r,
        azi / (2 * np.pi),
        x_coord,
        y_coord
    ]
    muon = Muon()
    mock_log_normal.assert_called_once_with(6.55, 1.8)
    mock_cos.assert_called_once_with(1.0)
    assert mock_random.call_count == 6
    assert muon.energy == mock_log_normal.return_value
    assert muon.zen == zen
    assert muon.azi == azi
    assert muon.x_coord == x_coord * 100 - 50 
    assert muon.y_coord == y_coord * 80 - 15
