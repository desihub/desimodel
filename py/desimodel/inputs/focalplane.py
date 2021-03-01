# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
desimodel.inputs.focalplane
==============================

Utilities for constructing a focalplane model.
"""
import os
import datetime
import re
import csv
import yaml
import glob

import configobj
import numpy as np
from astropy.table import Table, Column

from desiutil.log import get_logger

from . import docdb

from ..io import datadir, findfile, load_focalplane

from .focalplane_utils import (
    compute_theta_phi_range,
    rotate_petals,
    create_nominal,
    device_compare,
    device_printdiff,
    update_exclusions,
    valid_states,
    create_tables,
    load_petal_fiber_map,
    propagate_state,
)


def devices_from_fiberpos(fp):
    """Populate focalplane properties from a fiberpos file.

    This is only used for testing consistency with the previous fiberpos
    files.  It should not be used for work with the real instrument.  The
    focalplane properties are modified in place.

    Args:
        fp (dict):  The focalplane dictionary.

    Returns:
        None

    """
    # Get the mapping AND device info from the fiberpos instead.
    # We CANNOT use the desimodel "load_fiberpos()" function here, since
    # it seems to strip out the ETC devices.  Manually read the file.
    fiberpos_file = findfile("focalplane/fiberpos-all.fits")
    fpos = Table.read(fiberpos_file)

    # There is no "PETAL_ID" for the fake fiberpos, so we use PETAL
    for pet, dev, devtyp, blk, fib, xoff, yoff in zip(
        fpos["PETAL"],
        fpos["DEVICE"],
        fpos["DEVICE_TYPE"],
        fpos["SLITBLOCK"],
        fpos["BLOCKFIBER"],
        fpos["X"],
        fpos["Y"],
    ):
        fp[pet][dev]["DEVICE_ID"] = "FAKE"
        fp[pet][dev]["DEVICE_TYPE"] = devtyp
        fp[pet][dev]["CABLE"] = 0
        fp[pet][dev]["CONDUIT"] = "NA"
        fp[pet][dev]["FWHM"] = 0.0
        fp[pet][dev]["FRD"] = 0.0
        fp[pet][dev]["ABS"] = 0.0
        fp[pet][dev]["SLITBLOCK"] = blk
        fp[pet][dev]["BLOCKFIBER"] = fib
        fp[pet][dev]["OFFSET_X"] = xoff
        fp[pet][dev]["OFFSET_Y"] = yoff
        fp[pet][dev]["OFFSET_T"] = 0.0
        fp[pet][dev]["OFFSET_P"] = 0.0
        fp[pet][dev]["MIN_T"] = 0.0
        fp[pet][dev]["MIN_P"] = 0.0
        fp[pet][dev]["LENGTH_R1"] = 0.0
        fp[pet][dev]["LENGTH_R2"] = 0.0
        if (devtyp == "POS") or (devtyp == "ETC"):
            # This is a positioner.
            fp[pet][dev]["MAX_T"] = 380.0
            fp[pet][dev]["MAX_P"] = 200.0
            fp[pet][dev]["LENGTH_R1"] = 3.0
            fp[pet][dev]["LENGTH_R2"] = 3.0
    return


def devices_from_files(
    fp, posdir=None, fillfake=False, fakeoffset=False, fibermaps=None
):
    """Populate focalplane properties from information in files.

    This populates the focalplane with device information gathered from
    the "pos_settings" files in svn and from the "Petal verification"
    files on DocDB.

    The focalplane dictionary is modified in place.

    Args:
        fp (dict):  The focalplane dictionary.
        posdir (str):  Directory containing the many positioner conf files.
        fillfake (bool):  If true, fill missing POS and ETC locations with
            a fake nominal positioner.
        fakeoffset (bool):  If true, use theta / phi offsets that matched very
            old versions of fiberassign.
        fibermaps (list):  (optional) Override list of tuples (DocDB number,
            DocDB version, DocDB csv file) of where to find the petal mapping
            files.

    Returns:
        None

    """
    log = get_logger()

    fp = load_petal_fiber_map(existing=fp, fibermaps=fibermaps)

    # Parse all the positioner files.
    pos = dict()
    pospat = re.compile(r"unit_(.*).conf")
    if posdir is not None:
        for root, dirs, files in os.walk(posdir):
            for f in files:
                posmat = pospat.match(f)
                if posmat is not None:
                    file_dev = posmat.group(1)
                    pfile = os.path.join(root, f)
                    print("parsing {}".format(pfile), flush=True)
                    cnf = configobj.ConfigObj(pfile, unrepr=True)
                    # Is this device used?
                    if ("DEVICE_LOC" not in cnf) or (int(cnf["DEVICE_LOC"]) < 0):
                        continue
                    if "PETAL_ID" not in cnf:
                        continue
                    pet = int(cnf["PETAL_ID"])
                    if (pet < 0) or (pet not in fp):
                        continue
                    # Check that the positioner ID in the file name matches
                    # the file contents.
                    if file_dev != cnf["POS_ID"]:
                        msg = (
                            "positioner file {} has device {} in its name "
                            "but contains POS_ID={}".format(f, file_dev, cnf["POS_ID"])
                        )
                        raise RuntimeError(msg)
                    # Add properties to dictionary
                    pos[file_dev] = cnf
            break
    for devid, props in pos.items():
        pet = int(props["PETAL_ID"])
        dev = int(props["DEVICE_LOC"])
        if dev not in fp[pet]:
            # This should never happen- all possible device locations
            # should have been pre-populated before calling this function.
            msg = "Device location {} on petal ID {} does not exist".format(dev, pet)
            raise RuntimeError(msg)
        fp[pet][dev]["DEVICE_ID"] = devid
        t_min, t_max, p_min, p_max = compute_theta_phi_range(
            props["PHYSICAL_RANGE_T"], props["PHYSICAL_RANGE_P"]
        )
        fp[pet][dev]["OFFSET_T"] = props["OFFSET_T"]
        fp[pet][dev]["OFFSET_P"] = props["OFFSET_P"]
        fp[pet][dev]["LENGTH_R1"] = props["LENGTH_R1"]
        fp[pet][dev]["LENGTH_R2"] = props["LENGTH_R2"]
        fp[pet][dev]["MIN_T"] = t_min
        fp[pet][dev]["MAX_T"] = t_max
        fp[pet][dev]["MIN_P"] = p_min
        fp[pet][dev]["MAX_P"] = p_max

    if fillfake:
        t_min, t_max, p_min, p_max = compute_theta_phi_range(380.0, 200.0)
        for petal in list(sorted(fp.keys())):
            devlist = list(sorted(fp[petal].keys()))
            for dev in devlist:
                if (petal not in fp) or (dev not in fp[petal]):
                    print(fp[petal][dev], flush=True)
                devtyp = fp[petal][dev]["DEVICE_TYPE"]
                if (devtyp != "POS") and (devtyp != "ETC"):
                    continue
                if fp[petal][dev]["DEVICE_ID"] == "NONE":
                    fp[petal][dev]["LENGTH_R1"] = 3.0
                    fp[petal][dev]["LENGTH_R2"] = 3.0
                    if fakeoffset:
                        fp[petal][dev]["OFFSET_T"] = 0.0
                        fp[petal][dev]["OFFSET_P"] = 0.0
                        fp[petal][dev]["MIN_T"] = 0.0
                        fp[petal][dev]["MAX_T"] = 380.0
                        fp[petal][dev]["MIN_P"] = 0.0
                        fp[petal][dev]["MAX_P"] = 200.0
                    else:
                        fp[petal][dev]["OFFSET_T"] = -170.0
                        fp[petal][dev]["OFFSET_P"] = -5.0
                        fp[petal][dev]["MIN_T"] = t_min
                        fp[petal][dev]["MAX_T"] = t_max
                        fp[petal][dev]["MIN_P"] = p_min
                        fp[petal][dev]["MAX_P"] = p_max
    return


def create(
    testdir=None,
    posdir=None,
    fibermaps=None,
    petalloc=None,
    startvalid=None,
    fillfake=False,
    fakeoffset=False,
    fakefiberpos=False,
    reset=False,
):
    """Construct DESI focalplane and state files.

    This function gathers information from the following sources:
        - Petal verification files on DocDB
        - Positioner device configuration files (e.g. from svn).
        - DESI-0530, to get the mapping from device ID to device type as well
          as the nominal device X/Y offsets on petal 0 (for fillfake option).
        - Exclusion configobj files in $DESIMODEL/data/focalplane.

    Args:
        testdir (str):  Override the output directory for testing.
        posdir (str):  Directory containing the many positioner conf files.
            If None, simulate identical, nominal positioners.  A None value
            will force fillfake=True.
        fibermaps (list):  Override list of tuples (DocDB number,
            DocDB version, DocDB csv file) of where to find the petal mapping
            files.
        petalloc (dict):  Mapping of petal ID to petal location.
        startvalid (str):  The first time when this focalplane model is valid.
            ISO 8601 format string.
        fillfake (bool):  If True, fill missing device locations with fake
            positioners with nominal values for use in simulations.
        fakeoffset (bool):  If True, artificially sets the theta / phi angle
            offsets to zero.  This replicates the behavior of legacy
            fiberassign and should only be used for testing.
        fakefiberpos (bool):  If True, ignore the real fibermaps and load the
            old fiberpos file to get the mapping.  Only useful for testing.
        reset (bool):  If True, ignore all previous focalplane models and
            start with all positioners "good".  Default propagates the state
            of most recent model, after verifying that the positioners are the
            same.

    Returns:
        None

    """
    log = get_logger()

    outdir = testdir
    if outdir is None:
        outdir = os.path.join(datadir(), "focalplane")

    if posdir is None:
        fillfake = True

    if fakefiberpos and (posdir is not None):
        raise RuntimeError(
            "Cannot specify both fake positioners from fiberpos and real"
            " devices from posdir"
        )

    if startvalid is None:
        startvalid = datetime.datetime.utcnow()
    else:
        startvalid = datetime.datetime.strptime(startvalid, "%Y-%m-%dT%H:%M:%S")
    file_date = startvalid.isoformat(timespec="seconds")

    if (petalloc is None) and (posdir is not None):
        raise RuntimeError(
            "If specifying posdir, must also specify the petal locations."
        )

    # The mapping of petal IDs to locations
    if petalloc is None:
        petalloc = {x: x for x in range(10)}

    # Create a focalplane containing all possible petal locations
    # and devices.
    fp = create_nominal(petalloc)

    if fakefiberpos:
        devices_from_fiberpos(fp)
    else:
        devices_from_files(
            fp,
            posdir=posdir,
            fillfake=fillfake,
            fakeoffset=fakeoffset,
            fibermaps=fibermaps,
        )

    # Now rotate the X / Y offsets based on the petal location.
    rotate_petals(fp)

    # Construct the focaplane and state tables

    nrows = 0
    allpetals = list(sorted(fp.keys()))
    for petal in allpetals:
        devlist = list(sorted(fp[petal].keys()))
        nrows += len(devlist)

    out_fp, out_state = create_tables(nrows)

    # Populate the table
    dev_cols = [
        "OFFSET_X",
        "OFFSET_Y",
        "OFFSET_T",
        "OFFSET_P",
        "LENGTH_R1",
        "LENGTH_R2",
        "MIN_T",
        "MAX_T",
        "MIN_P",
        "MAX_P",
    ]

    row = 0
    for petal in allpetals:
        devlist = list(sorted(fp[petal].keys()))
        for dev in devlist:
            out_fp[row]["PETAL"] = fp[petal][dev]["PETAL"]
            out_fp[row]["DEVICE"] = dev
            out_fp[row]["LOCATION"] = fp[petal][dev]["PETAL"] * 1000 + dev
            out_fp[row]["PETAL_ID"] = petal
            out_fp[row]["DEVICE_ID"] = fp[petal][dev]["DEVICE_ID"]
            out_fp[row]["DEVICE_TYPE"] = fp[petal][dev]["DEVICE_TYPE"]
            out_fp[row]["SLITBLOCK"] = fp[petal][dev]["SLITBLOCK"]
            out_fp[row]["BLOCKFIBER"] = fp[petal][dev]["BLOCKFIBER"]
            out_fp[row]["CABLE"] = fp[petal][dev]["CABLE"]
            out_fp[row]["CONDUIT"] = fp[petal][dev]["CONDUIT"]
            out_fp[row]["FWHM"] = fp[petal][dev]["FWHM"]
            out_fp[row]["FRD"] = fp[petal][dev]["FRD"]
            out_fp[row]["ABS"] = fp[petal][dev]["ABS"]
            if fp[petal][dev]["SLITBLOCK"] < 0:
                # This must not be a POS device
                out_fp[row]["FIBER"] = -1
            else:
                out_fp[row]["FIBER"] = (
                    fp[petal][dev]["PETAL"] * 500
                    + fp[petal][dev]["SLITBLOCK"] * 25
                    + fp[petal][dev]["BLOCKFIBER"]
                )
            if (fp[petal][dev]["DEVICE_ID"] != "NONE") or fillfake:
                for col in dev_cols:
                    if col in fp[petal][dev]:
                        out_fp[row][col] = fp[petal][dev][col]
            out_state[row]["LOCATION"] = out_fp[row]["LOCATION"]
            out_state[row]["STATE"] = valid_states["OK"]
            out_state[row]["EXCLUSION"] = "default"
            out_state[row]["TIME"] = file_date
            row += 1

    # Now load the file(s) with the exclusion polygons
    # Add the legacy polygons to the dictionary for reference.
    # Also add an "unknown" polygon set which includes a large circle for the
    # theta arm that is the size of the patrol radius.
    excl = dict()

    # First the THETA arm.
    circs = [[[0.0 + 3.0, 0.0], 2.095]]
    seg = [
        [2.095 + 3.0, -0.474],
        [1.358 + 3.0, -2.5],
        [-0.229 + 3.0, -2.5],
        [-1.241 + 3.0, -2.792],
        [-2.095 + 3.0, -0.356],
    ]
    segs = [seg]
    shp_theta = dict()
    shp_theta["circles"] = circs
    shp_theta["segments"] = segs

    # Now the PHI arm
    circs = [[[0.0 + 3.0, 0.0], 0.967]]
    seg_upper = [[-3.0 + 3.0, 0.990], [0.0 + 3.0, 0.990]]
    seg_lower = [
        [-2.944 + 3.0, -1.339],
        [-2.944 + 3.0, -2.015],
        [-1.981 + 3.0, -1.757],
        [-1.844 + 3.0, -0.990],
        [0.0 + 3.0, -0.990],
    ]
    segs = [seg_upper, seg_lower]
    shp_phi = dict()
    shp_phi["circles"] = circs
    shp_phi["segments"] = segs

    excl["legacy"] = dict()
    excl["legacy"]["theta"] = shp_theta
    excl["legacy"]["phi"] = shp_phi

    excl["unknown"] = dict()
    excl["unknown"]["theta"] = dict()
    excl["unknown"]["theta"]["circles"] = [[[0.0, 0.0], 6.0]]
    excl["unknown"]["theta"]["segments"] = list()
    excl["unknown"]["phi"] = dict()
    excl["unknown"]["phi"]["circles"] = list()
    excl["unknown"]["phi"]["segments"] = list()

    # Get all available exclusion polygons from the desimodel data directory.

    fpdir = os.path.join(datadir(), "focalplane")
    excl_match = os.path.join(fpdir, "exclusions_*.conf")
    excl_files = glob.glob(excl_match)
    update_exclusions(excl, excl_files)

    # Propagate the device state.  Unless we have the reset option, we need
    # to load the current focalplane model and try to use the state from that.
    # We will choose a time that is one second before the currently selected
    # time.

    if not reset:
        dtime = datetime.timedelta(seconds=1)
        oldtime = startvalid - dtime

        oldfp, oldexcl, oldstate, oldtmstr = load_focalplane(oldtime)
        checkrows = np.where(oldstate["LOCATION"] == 7000)[0]
        print(oldstate[checkrows])

        log.info("Comparing generated focalplane to one from %s", oldtmstr)

        # Compare the old and new.  These are the device properties we care
        # about when propagating state.  In particular if positioner
        # calibration has modified arm lengths and angle ranges we don't
        # care.
        checkcols = [
            "PETAL",
            "DEVICE",
            "PETAL_ID",
            "DEVICE_ID",
            "DEVICE_TYPE",
            "SLITBLOCK",
            "BLOCKFIBER",
        ]
        diff = device_compare(oldfp, out_fp, checkcols)

        device_printdiff(diff)

        if len(diff) > 0:
            msg = (
                "Existing focalplane device properties have changed."
                "  Refusing to propagate the device state.  Use the 'reset'"
                " option to start with a new device state."
            )
            raise RuntimeError(msg)

        propagate_state(out_state, excl, oldstate, oldexcl)

    # Ensure that the default polygon has been defined.
    if "default" not in excl.keys():
        raise RuntimeError("No default exclusion polygon found in available files")

    # Now write out all of this collected information.  Also write out an
    # initial "state" log as a starting point.  Note that by having log
    # files (which contain datestamps) also have a "starting" date, it means
    # that we don't need a single log for the entire survey.

    out_fp_file = os.path.join(outdir, "desi-focalplane_{}.ecsv".format(file_date))
    out_excl_file = os.path.join(outdir, "desi-exclusion_{}.yaml".format(file_date))
    out_state_file = os.path.join(outdir, "desi-state_{}.ecsv".format(file_date))

    out_fp.write(out_fp_file, format="ascii.ecsv", overwrite=True)
    del out_fp

    out_state.write(out_state_file, format="ascii.ecsv", overwrite=True)
    del out_state

    # Now write out the exclusion polygons.  Since these are not tabular, we
    # write to a YAML file.

    with open(out_excl_file, "w") as pf:
        yaml.dump(excl, pf, default_flow_style=False)

    return
