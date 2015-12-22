# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
desimodel.install
=================

Install data files not handled by pip install.
"""
#
from __future__ import absolute_import, division, print_function, unicode_literals
#
#
#
def main():
    """Entry point for the install_desimodel_data script.
    """
    from . import __version__ as desimodel_version
    from sys import argv
    from subprocess import Popen, PIPE
    from argparse import ArgumentParser
    from os import chdir, environ
    from os.path import basename, dirname, exists, join
    default_install_dir = dirname(dirname(dirname(dirname(dirname(__file__)))))
    desc = """Install desimodel data.

This script will attempt to download and install the desimodel data/ directory.
The script will attempt to attempt to install the data in the following
locations, in order of preference:

1. $DESIMODEL, that is, the directory specified by the environment variable.
2. The value set with the -d option on the command line.
3. A directory relative to the file containing this script.  This directory
   is currently {0}.

If the data directory already exists, this script will not do anything.
""".format(default_install_dir)
    parser = ArgumentParser(description=desc, prog=argv[0])
    parser.add_argument('-d', '--desimodel', action='store', dest='desimodel',
        metavar='DESIMODEL',
        help=('Place the data/ directory in this directory.  '+
        'In other words, the environment variable DESIMODEL should be set to this directory.'))
    options = parser.parse_args()
    try:
        install_dir = environ['DESIMODEL']
    except KeyError:
        if options.desimodel is not None:
            install_dir = options.desimodel
        else:
            install_dir = default_install_dir
    if exists(join(install_dir, 'data')):
        print("{0} already exists!".format(join(install_dir, 'data')))
        return 1
    chdir(install_dir)
    command = ["svn", "export", "https://desi.lbl.gov/svn/code/desimodel/tags/{0}/data".format(desimodel_version)]
    # print(' '.join(command))
    proc = Popen(command, stdout=PIPE, stderr=PIPE)
    out, err = proc.communicate()
    # if basename(install_dir) == desimodel_version:
    #     print(install_dir)
    # print(options.datahome)
    return 0
