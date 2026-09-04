# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""pytest configuration for desimodel tests.

Some tests need $DESI_SURVEYOPS to point at a checkout of the surveyops
tiles files.  Rather than downloading a full snapshot from a remote server
during CI (which is prone to intermittent network failures), fall back to a
minimal fixture bundled with the package.  A real $DESI_SURVEYOPS snapshot,
if already present in the environment, always takes precedence over the
bundled fixture data.

This module is imported by pytest before any test modules in this directory
are collected, so setting $DESI_SURVEYOPS here (at import time, not just in
a fixture) ensures it is already set by the time test modules check for it
at import time, e.g. in test_io.py.
"""
import os

_bundled_surveyops = os.path.join(os.path.dirname(__file__), 'data', 'surveyops')


def _has_real_surveyops():
    """Return ``True`` if $DESI_SURVEYOPS already points at real tile data."""
    surveyops = os.environ.get('DESI_SURVEYOPS')
    if surveyops is None:
        return False
    for opsdir in (os.path.join(surveyops, 'ops'), os.path.join(surveyops, 'trunk', 'ops')):
        if os.path.isfile(os.path.join(opsdir, 'tiles-main.ecsv')):
            return True
    return False


if not _has_real_surveyops():
    os.environ['DESI_SURVEYOPS'] = _bundled_surveyops
