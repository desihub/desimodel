#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst
#
# Standard imports
#
import glob
import os
import sys
#
# setuptools' sdist command ignores MANIFEST.in
#
from distutils.command.sdist import sdist as DistutilsSdist
from setuptools import setup, find_packages
#
# DESI support code.
#
import desiutil.setup as ds
#
# First provide helpful messages if contributors try and run legacy commands
# for tests or docs.
#
API_HELP = """
Note: Generating api.rst files is no longer done using 'python setup.py api'. Instead
you will need to run:

    desi_api_file

which is part of the desiutil package. If you don't already have desiutil installed, you can install it with:

    pip install desiutil
"""

MODULE_HELP = """
Note: Generating Module files is no longer done using 'python setup.py api'. Instead
you will need to run:

    desiInstall

or

    desi_module_file

depending on your exact situation.  desiInstall is preferred.  Both commands are
part of the desiutil package. If you don't already have desiutil installed, you can install it with:

    pip install desiutil
"""

VERSION_HELP = """
Note: Generating version strings is no longer done using 'python setup.py version'. Instead
you will need to run:

    desi_update_version [-t TAG] desiutil

which is part of the desiutil package. If you don't already have desiutil installed, you can install it with:

    pip install desiutil
"""

TEST_HELP = """
Note: running tests is no longer done using 'python setup.py test'. Instead
you will need to run:

    pytest

If you don't already have pytest installed, you can install it with:

    pip install pytest
"""

DOCS_HELP = """
Note: building the documentation is no longer done using
'python setup.py {0}'. Instead you will need to run:

    sphinx-build -W --keep-going -b html doc doc/_build/html

If you don't already have Sphinx installed, you can install it with:

    pip install Sphinx
"""

message = {'api': API_HELP,
           'module_file': MODULE_HELP,
           'test': TEST_HELP,
           'version': VERSION_HELP,
           'build_docs': DOCS_HELP.format('build_docs'),
           'build_sphinx': DOCS_HELP.format('build_sphinx'), }

for m in message:
    if m in sys.argv:
        print(message[m])
        sys.exit(1)

#
# Begin setup
#
setup_keywords = dict()
#
# THESE SETTINGS NEED TO BE CHANGED FOR EVERY PRODUCT.
#
setup_keywords['name'] = 'desimodel'
setup_keywords['description'] = 'DESI hardware information'
setup_keywords['author'] = 'DESI Collaboration'
setup_keywords['author_email'] = 'desi-data@desi.lbl.gov'
setup_keywords['license'] = 'BSD'
setup_keywords['url'] = 'https://github.com/desihub/desimodel'
#
# END OF SETTINGS THAT NEED TO BE CHANGED.
#
setup_keywords['version'] = ds.get_version(setup_keywords['name'])
#
# Use README.rst as long_description.
#
setup_keywords['long_description'] = ''
if os.path.exists('README.rst'):
    with open('README.rst') as readme:
        setup_keywords['long_description'] = readme.read()
#
# Set other keywords for the setup function.  These are automated, & should
# be left alone unless you are an expert.
#
# Treat everything in bin/ except *.rst as a script to be installed.
#
if os.path.isdir('bin'):
    setup_keywords['scripts'] = [fname for fname in glob.glob(os.path.join('bin', '*'))
        if not os.path.basename(fname).endswith('.rst')]
setup_keywords['provides'] = [setup_keywords['name']]
setup_keywords['python_requires'] = '>=3.5'
setup_keywords['zip_safe'] = False
# setup_keywords['use_2to3'] = False
setup_keywords['packages'] = find_packages('py')
setup_keywords['package_dir'] = {'': 'py'}
setup_keywords['cmdclass'] = {'sdist': DistutilsSdist}
setup_keywords['test_suite']='{name}.test.{name}_test_suite'.format(**setup_keywords)
#
# Autogenerate command-line scripts.
#
setup_keywords['entry_points'] = {'console_scripts':['install_desimodel_data = desimodel.install:main']}
#
# Add internal data directories.
#
# setup_keywords['package_data'] = {'desimodel.test': ['t/*']}
#
# Run setup command.
#
setup(**setup_keywords)
