# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
desimodel.install
=================

Install data files not handled by pip install.
"""

import re

def default_install_dir():
    """Return the default install directory.

    Returns
    -------
    :class:`str`
        The path to the install directory.
    """
    import importlib.resources
    return str(importlib.resources.files('desimodel'))


def assert_svn_exists():
    """Assert svn command exists and raise an informative error if not"""

    from subprocess import check_output, CalledProcessError
    try:
        r = check_output(['svn', '--version'])
    except OSError as e:
        raise AssertionError("svn command is not executable. Install svn to use the install script. Original Error is: '{0}'.".format(e.strerror))
    except CalledProcessError as e:
        raise AssertionError("The svn command ({0}) on this system does not work. Output is: '{1}'.".format(e.cmd, e.output))

def get_svn_version(desimodel_version=None):
    """
    Return which svn version should be checked out for given desimodel_version

    Parameters
    ----------
    desimodel_version : :class:`str`, optional
        The version X.Y.Z to download, trunk, or something of the
        form branches/... Defaults to package version if x.y.z tagged,
        otherwise trunk.

    Returns
    -------
    :class:`str`
        svn version string
    """

    if desimodel_version is None:
        from . import __version__ as this_version
        # check for tag-like versions x.y.z or x.y
        if (re.match(r'^\d+\.\d+\.\d+$', this_version) is not None or
            re.match(r'\d+\.\d+$', this_version) is not None):
            svn_version = 'tags/' + this_version
        # otherwise use trunk if version wasn't specified
        else:
            svn_version = 'trunk'
    elif (desimodel_version == 'trunk') or ('branches/' in desimodel_version):
        svn_version = desimodel_version
    else:
        svn_version = 'tags/' + desimodel_version

    return svn_version

def svn_export(desimodel_version=None, svn_checkout=False):
    """Create a :command:`svn export` command suitable for downloading a
    particular desimodel version.

    Parameters
    ----------
    desimodel_version : :class:`str`, optional
        The version X.Y.Z to download, trunk, or something of the
        form branches/... Defaults to package version if x.y.z tagged,
        otherwise trunk.
    svn_checkout : bool, default False
        If True, svn checkout instead of svn export

    Returns
    -------
    :class:`list`
        A :command:`svn` command in list form, suitable for passing to
        :class:`subprocess.Popen`.
    """
    svn_version = get_svn_version(desimodel_version)
    if svn_checkout:
        svn_subcommand = 'checkout'
    else:
        svn_subcommand = 'export'

    return ["svn", svn_subcommand,
            f"https://desi.lbl.gov/svn/code/desimodel/{svn_version}/data"]


def install(desimodel=None, version=None, svn_checkout=False, dry_run=False):
    """Primary workhorse function.

    Parameters
    ----------
    desimodel : :class:`str`, optional
        Allows the install directory to be explicitly set.
    version : :class:`str`, optional
        Allows the desimodel version to be explicitly set.
    svn_checkout : bool, default False
        If True, svn checkout instead of svn export
    dry_run : bool, default False
        If True, print commands but don't actually get the data

    Raises
    ------
    :class:`RuntimeError`
        Standard error output from svn export command when status is non-zero.
    """
    from os import chdir, environ
    from os.path import exists, join
    from subprocess import Popen, PIPE
    try:
        install_dir = environ['DESIMODEL']
    except KeyError:
        if desimodel is not None:
            install_dir = desimodel
        else:
            install_dir = default_install_dir()
    if exists(join(install_dir, 'data')):
        raise ValueError("{0} already exists!".format(join(install_dir,
                                                           'data')))
    assert_svn_exists()

    chdir(install_dir)

    svn_version = get_svn_version(version)
    print(f'Installing desimodel data {svn_version} to {install_dir}')

    command = svn_export(version, svn_checkout)
    if dry_run:
        cmdstr = ' '.join(command)
        print(f'Dry run, would have run "{cmdstr}"')
    else:
        proc = Popen(command, stdout=PIPE, stderr=PIPE)
        out, err = proc.communicate()
        status = proc.returncode
        if status != 0:
            raise RuntimeError(err.rstrip())


def main():
    """Entry point for the :command:`install_desimodel_data` script.

    Returns
    -------
    :class:`int`
        An integer suitable for passing to :func:`sys.exit`.
    """
    from sys import argv
    from argparse import ArgumentParser
    desc = """Install desimodel data.

This script will attempt to download and install the desimodel data/ directory.
The script will attempt to attempt to install the data in the following
locations, in order of preference:

1. :envvar:`DESIMODEL`, that is, the directory specified by the
   environment variable.
2. The value set with the -d option on the command line.
3. A directory relative to the file containing this script.  This directory
   is currently {0}.

If the data directory already exists, this script will not do anything.
""".format(default_install_dir())
    parser = ArgumentParser(description=desc, prog=argv[0])
    parser.add_argument('-d', '--desimodel', action='store', dest='desimodel',
                        metavar='DESIMODEL',
                        help=('Place the data/ directory in this directory. ' +
                              'In other words, the environment variable ' +
                              'DESIMODEL should be set to this directory.'))
    parser.add_argument('-D', '--desimodel-version', action='store',
                        dest='desimodel_version', metavar='VERSION',
                        help='Explicitly set the version to download.')
    parser.add_argument('--checkout', action='store_true',
                        help='svn checkout instead of svn export data')
    parser.add_argument('--dry-run', action='store_true',
                        help="Print actions but don't actually install the data")
    options = parser.parse_args()
    try:
        install(options.desimodel, options.desimodel_version,
                svn_checkout=options.checkout, dry_run=options.dry_run)
    except (ValueError, RuntimeError) as e:
        print(e.message)
        return 1
    return 0
