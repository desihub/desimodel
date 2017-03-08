===============================
How to update inputs from DocDB
===============================

Basic Setup
===========

Make branches of both desimodel github code and svn data::

    git checkout -b update_inputs
    base=https://desi.lbl.gov/svn/code/desimodel
    svn copy $base/trunk $base/branches/update_inputs
    svn checkout $base/branches/update_inputs/data
    export DESIMODEL=`pwd`
    export docdb=https://desi.lbl.gov/DocDB/cgi-bin/private/RetrieveFile

Add entries to $HOME/.netrc to enable downloading DocDB files without
having to enter a password every time::

    machine desi.lbl.gov
    login StephenBailey
    password NotMyRealPassword

Inputs to update
================

::

    data/inputs/
        fiberpos/
            pos_on_z1.txt               - DEPRECATED
            cassette_order.txt          - from DESI-2721
            DESI-0530-posloc.txt        - from DESI-0530
        throughput/
            DESI-0334-spectro/
                blue-thru.txt
                nir-thru-500.txt
                red-thru.txt
            DESI-0347-throughput.txt


Overall DESI parameters
-----------------------

Get the overall DESI parameters file from DESI-0347::

    DOC347VER=11
    cd $DESIMODEL/data/
    wget -O desi.yaml "$docdb?docid=347;filename=desi.yaml;version=$DOC347VER"

Mapping of fiber positions to focal plane locations
===================================================

In data/inputs/fiberpos:

  * Update cassette_order.txt from DESI-2721
  * Update DESI-0530-posloc.txt from DESI-0530 tab
    PositionerAndFiducialLocations cells B40-J589
  * Combine these into the output data/focalplane/fiberpos* files

::

    python $DESIMODEL/bin/randomize_fibers.py \
        --input $DESIMODEL/data/inputs/fiberpos/DESI-0530-posloc.txt \
        --cassettes $DESIMODEL/data/inputs/fiberpos/cassette_order.txt \
        --outdir $DESIMODEL/data/focalplane

See fiberpos.rst for more details.

Throughput
==========

CCD Quantum Efficiency
----------------------

In data/inputs/throughput/DESI-0334-spectro/, download the following files
from https://desi.lbl.gov/DocDB/cgi-bin/private/ShowDocument?docid=334:

  * blue-thru.txt
  * red-thru.txt
  * nir-thru-250.txt

See $DESIMODEL/data/inputs/throughput/DESI-0334-spectro/README for cut-and-paste
wget instructions.

Overall DESI throughput
-----------------------

In data/inputs/throughput, download https://desi.lbl.gov/DocDB/cgi-bin/private/RetrieveFile?docid=347;filename=Inst_Throughput_forDESIModel11.txt;version=11
as DESI-0347-throughput.txt, e.g.::

    DOC347VER=11
    cd $DESIMODEL/data/inputs/throughput
    wget -O DESI-0347-throughput.txt "$docdb?docid=347;filename=Inst_Throughput_forDESIModel11.txt;version=11"

Combine throughput files
------------------------

::

  python $DESIMODEL/bin/combine_throughputs.py -o $DESIMODEL/data/throughput

Other Inputs
============

Other input files not covered here, added as part of
`desimodel PR #29 <https://github.com/desihub/desimodel/pull/29>`_:

  * Data for the figures on the "Geometric Blur" tab of DESI-0347,
    though not yet in DocDB as a data file itself:

    * inputs/throughput/raytracing.txt

  * Reformatted by hand from DESI-0347 excel spreadsheet:

    * throughput/DESI-0347_blur.ecsv
    * throughput/DESI-0347_offset.ecsv
