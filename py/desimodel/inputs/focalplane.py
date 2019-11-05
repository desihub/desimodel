# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
desimodel.inputs.focalplane
==============================

Utilities for constructing a focalplane model.
"""
import os
from datetime import datetime
import re
import csv
import yaml

import configobj
import numpy as np
from astropy.table import Table, Column

from desiutil.log import get_logger

from . import docdb

from ..io import datadir, findfile


def _compute_theta_phi_range(phys_t, phys_p):
    """Compute the min/max range about the initial offset.

    Based on the "full_range" defined in plate_control/petal/posmodel.py

    Args:
        phys_t (float):  PHYSICAL_RANGE_T in degrees.
        phys_p (float):  PHYSICAL_RANGE_P in degrees.

    Returns:
        (tuple):  The (theta_min, theta_max, phi_min, phi_max) angles.

    """
    t_min = -0.5 * phys_t
    t_max = 0.5 * phys_t
    p_min = 185.0 - phys_p
    p_max = 185.0
    return (t_min, t_max, p_min, p_max)


def _create_device():
    """Create an empty device property dictionary.
    """
    props = dict()
    props["PETAL"] = -1
    props["DEVICE"] = -1
    props["PETAL_ID"] = -1
    props["DEVICE_ID"] = "NONE"
    props["DEVICE_TYPE"] = "NONE"
    props["CABLE"] = -1
    props["CONDUIT"] = "NA"
    props["FWHM"] = 0.0
    props["FRD"] = 0.0
    props["ABS"] = 0.0
    props["OFFSET_X"] = 0.0
    props["OFFSET_Y"] = 0.0
    props["SLITBLOCK"] = -1
    props["BLOCKFIBER"] = -1
    props["OFFSET_T"] = 0.0
    props["OFFSET_P"] = 0.0
    props["MIN_T"] = 0.0
    props["MAX_T"] = 0.0
    props["MIN_P"] = 0.0
    props["MAX_P"] = 0.0
    props["LENGTH_R1"] = 0.0
    props["LENGTH_R2"] = 0.0
    return props


def _rotate_petals(fp):
    """Rotate the X/Y offsets according to petal location.

    The X / Y offsets of each device are rotated to the petal
    location for that device.  The focalplane dictionary is
    modified in place.

    Args:
        fp (dict):  The focalplane dictionary.

    Returns:
        None

    """
    # Now rotate the X / Y offsets based on the petal location.
    petals = list(sorted(fp.keys()))
    for petal in petals:
        devlist = list(sorted(fp[petal].keys()))
        for dev in devlist:
            # The petal location of this petal ID
            petal_loc = fp[petal][dev]["PETAL"]
            # Petal 0 is at the "bottom"; See DESI-0530.  The X/Y and
            # positioner theta offset are defined with the petal in location 3
            # We need to rotate from petal location 3 to desired location.
            petalrot_deg = (float(7 + petal_loc) * 36.0) % 360.0
            petalrot_rad = np.radians(petalrot_deg)
            x = fp[petal][dev]["OFFSET_X"]
            y = fp[petal][dev]["OFFSET_Y"]
            fp[petal][dev]["OFFSET_X"] = \
                np.cos(petalrot_rad) * x - np.sin(petalrot_rad) * y
            fp[petal][dev]["OFFSET_Y"] = \
                np.sin(petalrot_rad) * x + np.cos(petalrot_rad) * y
            fp[petal][dev]["OFFSET_T"] += petalrot_deg
    return


def _create_nominal(petal_loc):
    """Create a nominal focalplane layout.

    This uses DocDB 0530 to construct a nominal focalplane.  All
    positioner devices are assigned to their nominal X/Y locations,
    nominal theta / phi offsets and ranges, and nominal arm lengths.

    Quantities not specified in 0530, such as physical petal and device IDs,
    are set to -1.

    The input petal_loc dictionary is required, and only petal IDs in this
    dictionary will be created.

    Note:  The X/Y offsets used are those relative to the petal
    when placed at location 3.  After modifying these with data
    from other sources, the final offsets are rotated into place.

    Args:
        petal_loc (dict):  Dictionary of petal ID to petal location.

    Returns:
        (dict):  Dictionary of petal location properties, containing
            dictionaries of device properties.

    """
    log = get_logger()
    fp = dict()

    xls_fp_layout = docdb.download(
        530, 14, "DESI-0530-v14 (Focal Plane Layout).xlsx")
    xls_sheet = "PositionerAndFiducialLocations"
    rowmin, rowmax = 49, 591
    headers = docdb.xls_read_row(xls_fp_layout, xls_sheet, rowmin-1, "B", "S")
    assert headers[0] == "device_location_id"
    assert headers[1] == "device_type"
    xls_devloc = docdb.xls_read_col(
        xls_fp_layout, xls_sheet, "B", rowmin, rowmax, dtype=np.int32)
    xls_devtype = docdb.xls_read_col(
        xls_fp_layout, xls_sheet, "C", rowmin, rowmax, dtype=str)
    xls_dev_nominal_x = docdb.xls_read_col(
        xls_fp_layout, xls_sheet, "D", rowmin, rowmax, dtype=np.float64)
    xls_dev_nominal_y = docdb.xls_read_col(
        xls_fp_layout, xls_sheet, "E", rowmin, rowmax, dtype=np.float64)
    devtype = dict()
    dev_nominal_xy = dict()
    for loc, typ in zip(xls_devloc, xls_devtype):
        devtype[int(loc)] = typ
    for loc, x, y in zip(xls_devloc, xls_dev_nominal_x, xls_dev_nominal_y):
        dev_nominal_xy[int(loc)] = (x, y)

    petals = list(sorted(petal_loc.keys()))
    device_locs = list(sorted(devtype.keys()))
    for petal in petals:
        pt = dict()
        for loc in device_locs:
            # Create an empty device
            props = _create_device()
            props["PETAL"] = petal_loc[petal]
            props["PETAL_ID"] = petal
            props["DEVICE_TYPE"] = devtype[loc]
            x, y = dev_nominal_xy[loc]
            props["OFFSET_X"] = x
            props["OFFSET_Y"] = y
            pt[loc] = props
        fp[petal] = pt
    return fp


def _devices_from_fiberpos(fp):
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
            fpos["PETAL"], fpos["DEVICE"], fpos["DEVICE_TYPE"],
            fpos["SLITBLOCK"], fpos["BLOCKFIBER"], fpos["X"], fpos["Y"]):
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
        if (devtyp != "POS") and (devtyp != "ETC"):
            # This is a positioner.
            fp[pet][dev]["MAX_T"] = 380.0
            fp[pet][dev]["MAX_P"] = 200.0
            fp[pet][dev]["LENGTH_R1"] = 3.0
            fp[pet][dev]["LENGTH_R2"] = 3.0
    return


def _devices_from_files(fp, posdir=None, fillfake=False, fakeoffset=False,
                        fibermaps=None):
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
    if fibermaps is None:
        fibermaps = [
            (4042, 5, "Petal_2_final_verification.csv"),
            (4043, 7, "Petal_3_final_verification.csv"),
            (4807, 2, "Petal_4_final_verification.csv"),
            (4808, 3, "Petal_5_final_verification.csv"),
            (4809, 2, "Petal_6_final_verification.csv"),
            (4190, 6, "Petal_7_final_verification.csv"),
            (4806, 4, "Petal_8_final_verification.csv"),
            (4810, 3, "Petal_9_final_verification.csv"),
            (4868, 5, "Petal_10_final_verification.csv"),
            (4883, 4, "Petal_11_final_verification.csv")
        ]
    for docnum, docver, docname in fibermaps:
        fmfile = None
        try:
            fmfile = docdb.download(docnum, docver, docname)
        except IOError:
            msg = "Could not download {}".format(docname)
            log.error(msg)
        fmslitcheck = dict()
        firstline = True
        with open(fmfile, newline="") as csvfile:
            reader = csv.reader(csvfile, delimiter=",")
            cols = dict()
            for row in reader:
                if firstline:
                    for cnum, elem in enumerate(row):
                        nm = elem.strip().rstrip()
                        cols[nm] = cnum
                    firstline = False
                else:
                    pet = int(row[cols["PETAL_ID"]])
                    dev = int(row[cols["DEVICE_LOC"]])
                    cable = int(row[cols["Cable_ID"]])
                    conduit = row[cols["Conduit"]]
                    fwhm = float(row[cols["FWHM@f/3.9"]])
                    fthrough = float(row[cols["FRD_Throughput"]])
                    athrough = float(row[cols["Abs_Throuhgput"]])
                    blkfib = row[cols["slit_position"]].split(":")
                    blk = int(blkfib[0])
                    fib = int(blkfib[1])
                    if blk not in fmslitcheck:
                        fmslitcheck[blk] = dict()
                    if fib in fmslitcheck[blk]:
                        msg = "Petal ID {}, slitblock {}, blockfiber {}" \
                            " already assigned to device {}.  " \
                            "Reassigning to device {}".format(
                                pet, blk, fib, fmslitcheck[blk][fib], dev
                            )
                        log.warning(msg)
                    fmslitcheck[blk][fib] = dev
                    if (pet not in fp) or (dev not in fp[pet]):
                        print("FAIL:  petal {}, dev {} not in fp".format(pet, dev), flush=True)
                    fp[pet][dev]["SLITBLOCK"] = blk
                    fp[pet][dev]["BLOCKFIBER"] = fib
                    fp[pet][dev]["CABLE"] = cable
                    fp[pet][dev]["CONDUIT"] = conduit
                    fp[pet][dev]["FWHM"] = fwhm
                    fp[pet][dev]["FRD"] = fthrough
                    fp[pet][dev]["ABS"] = athrough
    # HARD-CODED modifications.  These changes are to work around features
    # in the files on DocDB.  Remove these as they are fixed upstream.
    # Note that once we get information from the database, then these may
    # no longer be needed.
    # ---------------------------
    # DESI-4807v2-Petal_4_final_verification.csv
    # Petal ID 04 has a typo.  Device location 357 should be slitblock
    # 19 and blockfiber 23 (it is marked as 24)
    log.info("Correcting petal ID 4, location 357")
    fp[4][357]["SLITBLOCK"] = 19
    fp[4][357]["BLOCKFIBER"] = 23
    # ---------------------------
    # DESI-4809v2-Petal_6_final_verification.csv
    # Petal ID 06 is missing an entry for device location 261.  This device
    # location is assigned to positioner M03120 in the pos_settings files.
    # Assign it to the one missing fiber location.
    log.info("Correcting petal ID 6, location 261")
    fp[6][261]["DEVICE_TYPE"] = "POS"
    fp[6][261]["DEVICE_ID"] = "NONE" # Populated below from pos_settings
    fp[6][261]["SLITBLOCK"] = 19
    fp[6][261]["BLOCKFIBER"] = 22
    fp[6][261]["CABLE"] = 6
    fp[6][261]["CONDUIT"] = "E0"     # This conduit has one fewer than F3
    fp[6][261]["FWHM"] = 0.0         # No information from file
    fp[6][261]["FRD"] = 0.0          # No information from file
    fp[6][261]["ABS"] = 0.0          # No information from file
    # ---------------------------
    # DESI-4883v4-Petal_11_final_verification.csv
    # Petal ID 11 is missing an entry for device location 484.  This
    # device location is assigned to positioner M06847 in the pos_settings
    # files.  Assign it to the one missing fiber location.
    log.info("Correcting petal ID 11, location 484")
    fp[11][484]["DEVICE_TYPE"] = "POS"
    fp[11][484]["DEVICE_ID"] = "NONE"  # Populated below from pos_settings
    fp[11][484]["SLITBLOCK"] = 3
    fp[11][484]["BLOCKFIBER"] = 3
    fp[11][484]["CABLE"] = 4
    fp[11][484]["CONDUIT"] = "G0"     # This conduit has one fewer than G1
    fp[11][484]["FWHM"] = 0.0         # No information from file
    fp[11][484]["FRD"] = 0.0          # No information from file
    fp[11][484]["ABS"] = 0.0          # No information from file

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
                    if ("DEVICE_LOC" not in cnf) \
                            or (int(cnf["DEVICE_LOC"]) < 0):
                        continue
                    if ("PETAL_ID" not in cnf):
                        continue
                    pet = int(cnf["PETAL_ID"])
                    if (pet < 0) or (pet not in fp):
                        continue
                    # Check that the positioner ID in the file name matches
                    # the file contents.
                    if file_dev != cnf["POS_ID"]:
                        msg = "positioner file {} has device {} in its name "\
                            "but contains POS_ID={}"\
                            .format(f, file_dev, cnf["POS_ID"])
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
            msg = "Device location {} on petal ID {} does not exist".format(
                dev, pet
            )
            raise RuntimeError(msg)
        fp[pet][dev]["DEVICE_ID"] = devid
        t_min, t_max, p_min, p_max = _compute_theta_phi_range(
            props["PHYSICAL_RANGE_T"], props["PHYSICAL_RANGE_P"])
        fp[pet][dev]["OFFSET_T"] = props["OFFSET_T"]
        fp[pet][dev]["OFFSET_P"] = props["OFFSET_P"]
        fp[pet][dev]["LENGTH_R1"] = props["LENGTH_R1"]
        fp[pet][dev]["LENGTH_R2"] = props["LENGTH_R2"]
        fp[pet][dev]["MIN_T"] = t_min
        fp[pet][dev]["MAX_T"] = t_max
        fp[pet][dev]["MIN_P"] = p_min
        fp[pet][dev]["MAX_P"] = p_max

    if fillfake:
        t_min, t_max, p_min, p_max = _compute_theta_phi_range(380.0, 200.0)
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


def _collision_to_segments(raw):
    rx = raw[:, 0]
    ry = raw[:, 1]
    sg = [[float(x), float(y)] for x, y in zip(rx, ry)]
    start = list(sg[0])
    sg.append(start)
    return [sg]


def update_exclusions(excl, paths=list()):
    """Update exclusion polygons in a focalplane model.

    Args:
        excl (dict):  Dictionary of exclusion polygons, modified in place.
        paths (list):  List of file paths to append to the exclusions.

    Returns:
        None

    """
    log = get_logger()
    # NOTE:  The GFA and Petal exclusion polygons in these files are for
    # the petal in the default position (location 3).  They will be
    # rotated by downstream codes like fiberassign.  If the petal locations
    # in focalplane coordinates are very different from nominal, we may
    # want to read and store explicit polygons for each petal.  TBD.

    for pf in paths:
        # Add shapes from other files.
        log.info("Loading exclusion polygons from {}".format(pf))
        exprops = configobj.ConfigObj(pf, unrepr=True)
        if "NAME" not in exprops:
            msg = "exclusion file {} does not contain a NAME parameter"\
                .format(pf)
            raise RuntimeError(msg)
        nm = exprops["NAME"]
        props = dict()
        ktheta_raw = np.transpose(np.array(exprops["KEEPOUT_THETA"]))
        props["theta"] = dict()
        props["theta"]["segments"] = \
            _collision_to_segments(ktheta_raw)
        props["theta"]["circles"] = list()
        kphi_raw = np.transpose(np.array(exprops["KEEPOUT_PHI"]))
        props["phi"] = dict()
        props["phi"]["segments"] = _collision_to_segments(kphi_raw)
        props["phi"]["circles"] = list()
        kpetal_raw = np.transpose(np.array(exprops["KEEPOUT_PTL"]))
        props["petal"] = dict()
        props["petal"]["segments"] = \
            _collision_to_segments(kpetal_raw)
        props["petal"]["circles"] = list()
        kgfa_raw = np.transpose(np.array(exprops["KEEPOUT_GFA"]))
        props["gfa"] = dict()
        props["gfa"]["segments"] = \
            _collision_to_segments(kgfa_raw)
        props["gfa"]["circles"] = list()
        excl[nm] = props
    return


def create(testdir=None, posdir=None, polyfile=None, fibermaps=None,
           petalloc=None, startvalid=None, fillfake=False,
           exclusion="legacy", fakeoffset=False, fakefiberpos=False):
    """Construct DESI focalplane and state files.

    This function gathers information from the following sources:
        - Petal verification files on DocDB
        - Positioner device configuration files (e.g. from svn).
        - DESI-0530, to get the mapping from device ID to device type as well
          as the nominal device X/Y offsets on petal 0 (for fillfake option).
        - A "collision" file containing the exclusion polygons to use.

    Args:
        testdir (str):  Override the output directory for testing.
        posdir (str):  Directory containing the many positioner conf files.
            If None, simulate identical, nominal positioners.  A None value
            will force fillfake=True.
        polyfile (str):  File containing the exclusion polygons.  If None,
            Use the "legacy" polygons historically included in fiberassign.
        fibermaps (list):  Override list of tuples (DocDB number,
            DocDB version, DocDB csv file) of where to find the petal mapping
            files.
        petalloc (dict):  Mapping of petal ID to petal location.
        startvalid (str):  The first time when this focalplane model is valid.
            ISO 8601 format string.
        fillfake (bool):  If True, fill missing device locations with fake
            positioners with nominal values for use in simulations.
        exclusion (str):  The name of the default exclusion polygons.
        fakeoffset (bool):  If True, artificially sets the theta / phi angle
            offsets to zero.  This replicates the behavior of legacy
            fiberassign and should only be used for testing.
        fakefiberpos (bool):  If True, ignore the real fibermaps and load the
            old fiberpos file to get the mapping.  Only useful for testing.

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
            " devices from posdir")

    if startvalid is None:
        startvalid = datetime.utcnow()
    else:
        startvalid = datetime.strptime(startvalid, "%Y-%m-%dT%H:%M:%S")

    if (petalloc is None) and (posdir is not None):
        raise RuntimeError(
            "If specifying posdir, must also specify the petal locations."
        )

    # The mapping of petal IDs to locations
    if petalloc is None:
        petalloc = {x: x for x in range(10)}

    # Create a focalplane containing all possible petal locations
    # and devices.
    fp = _create_nominal(petalloc)

    # FIXME:  Add another option here to "use the database".

    if fakefiberpos:
        _devices_from_fiberpos(fp)
    else:
        _devices_from_files(
            fp, posdir=posdir, fillfake=fillfake, fakeoffset=fakeoffset,
            fibermaps=fibermaps
        )

    # Now rotate the X / Y offsets based on the petal location.
    _rotate_petals(fp)

    # Now load the file(s) with the exclusion polygons
    # Add the legacy polygons to the dictionary for reference.
    # Also add an "unknown" polygon set which includes a large circle for the
    # theta arm that is the size of the patrol radius.
    poly = dict()

    # First the THETA arm.
    circs = [
        [[0.0+3.0, 0.0], 2.095]
    ]
    seg = [
        [2.095+3.0, -0.474],
        [1.358+3.0, -2.5],
        [-0.229+3.0, -2.5],
        [-1.241+3.0, -2.792],
        [-2.095+3.0, -0.356]
    ]
    segs = [seg]
    shp_theta = dict()
    shp_theta["circles"] = circs
    shp_theta["segments"] = segs

    # Now the PHI arm
    circs = [
        [[0.0+3.0, 0.0], 0.967]
    ]
    seg_upper = [
        [-3.0+3.0, 0.990],
        [0.0+3.0, 0.990]
    ]
    seg_lower = [
        [-2.944+3.0, -1.339],
        [-2.944+3.0, -2.015],
        [-1.981+3.0, -1.757],
        [-1.844+3.0, -0.990],
        [0.0+3.0, -0.990]
    ]
    segs = [seg_upper, seg_lower]
    shp_phi = dict()
    shp_phi["circles"] = circs
    shp_phi["segments"] = segs

    poly["legacy"] = dict()
    poly["legacy"]["theta"] = shp_theta
    poly["legacy"]["phi"] = shp_phi

    poly["unknown"] = dict()
    poly["unknown"]["theta"] = dict()
    poly["unknown"]["theta"]["circles"] = [
        [[0.0, 0.0], 6.0]
    ]
    poly["unknown"]["theta"]["segments"] = list()
    poly["unknown"]["phi"] = dict()
    poly["unknown"]["phi"]["circles"] = list()
    poly["unknown"]["phi"]["segments"] = list()

    if polyfile is not None:
        update_exclusions(poly, polyfile)

    # Now write out all of this collected information.  Also write out an
    # initial "state" log as a starting point.  Note that by having log
    # files (which contain datestamps) also have a "starting" date, it means
    # that we don't need a single log for the entire survey.

    file_date = startvalid.isoformat(timespec="seconds")
    out_fp_file = os.path.join(
        outdir, "desi-focalplane_{}.ecsv".format(file_date))
    out_poly_file = os.path.join(
        outdir, "desi-exclusion_{}.yaml".format(file_date))
    out_state_file = os.path.join(
        outdir, "desi-state_{}.ecsv".format(file_date))

    # First the focalplane file

    nrows = 0
    allpetals = list(sorted(fp.keys()))
    for petal in allpetals:
        devlist = list(sorted(fp[petal].keys()))
        nrows += len(devlist)

    out_cols = [
        Column(name="PETAL", length=nrows, dtype=np.int32,
               description="Petal location [0-9]"),
        Column(name="DEVICE", length=nrows, dtype=np.int32,
               description="Device location on the petal"),
        Column(name="LOCATION", length=nrows, dtype=np.int32,
               description="PETAL * 1000 + DEVICE"),
        Column(name="PETAL_ID", length=nrows, dtype=np.int32,
               description="The physical petal ID"),
        Column(name="DEVICE_ID", length=nrows, dtype=np.dtype("a9"),
               description="The physical device ID string"),
        Column(name="DEVICE_TYPE", length=nrows, dtype=np.dtype("a3"),
               description="The device type (POS, ETC, FIF)"),
        Column(name="SLITBLOCK", length=nrows, dtype=np.int32,
               description="The slit block where this fiber goes"),
        Column(name="BLOCKFIBER", length=nrows, dtype=np.int32,
               description="The fiber index within the slit block"),
        Column(name="CABLE", length=nrows, dtype=np.int32,
               description="The cable ID"),
        Column(name="CONDUIT", length=nrows, dtype=np.dtype("a3"),
               description="The conduit"),
        Column(name="FIBER", length=nrows, dtype=np.int32,
               description="PETAL * 500 + SLITBLOCK * 25 + BLOCKFIBER"),
        Column(name="FWHM", length=nrows, dtype=np.float64,
               description="FWHM at f/3.9"),
        Column(name="FRD", length=nrows, dtype=np.float64,
               description="FRD Throughput"),
        Column(name="ABS", length=nrows, dtype=np.float64,
               description="ABS Throughput"),
        Column(name="OFFSET_X", length=nrows, dtype=np.float64,
               description="X location of positioner center", unit="mm"),
        Column(name="OFFSET_Y", length=nrows, dtype=np.float64,
               description="Y location of positioner center", unit="mm"),
        Column(name="OFFSET_T", length=nrows, dtype=np.float64,
               description="THETA zero point angle", unit="degrees"),
        Column(name="OFFSET_P", length=nrows, dtype=np.float64,
               description="PHI zero point angle", unit="degrees"),
        Column(name="LENGTH_R1", length=nrows, dtype=np.float64,
               description="Length of THETA arm", unit="mm"),
        Column(name="LENGTH_R2", length=nrows, dtype=np.float64,
               description="Length of PHI arm", unit="mm"),
        Column(name="MAX_T", length=nrows, dtype=np.float64,
               description="Maximum THETA angle relative to OFFSET_T",
               unit="degrees"),
        Column(name="MIN_T", length=nrows, dtype=np.float64,
               description="Minimum THETA angle relative to OFFSET_T",
               unit="degrees"),
        Column(name="MAX_P", length=nrows, dtype=np.float64,
               description="Maximum PHI angle relative to OFFSET_P",
               unit="degrees"),
        Column(name="MIN_P", length=nrows, dtype=np.float64,
               description="Minimum PHI angle relative to OFFSET_P",
               unit="degrees"),
    ]

    out_fp = Table()
    out_fp.add_columns(out_cols)

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
        "MAX_P"
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
                out_fp[row]["FIBER"] = \
                    fp[petal][dev]["PETAL"] * 500 \
                    + fp[petal][dev]["SLITBLOCK"] * 25 \
                    + fp[petal][dev]["BLOCKFIBER"]
            if (fp[petal][dev]["DEVICE_ID"] != "NONE") or fillfake:
                for col in dev_cols:
                    if col in fp[petal][dev]:
                        out_fp[row][col] = fp[petal][dev][col]
            row += 1

    out_fp.write(out_fp_file, format="ascii.ecsv")
    del out_fp

    # Now the state file.

    out_cols = [
        Column(name="TIME", length=nrows, dtype=np.dtype("a20"),
               description="The timestamp of the event (UTC, ISO format)"),
        Column(name="PETAL", length=nrows, dtype=np.int32,
               description="Petal location [0-9]"),
        Column(name="DEVICE", length=nrows, dtype=np.int32,
               description="Device location on the petal"),
        Column(name="LOCATION", length=nrows, dtype=np.int32,
               description="Global device location (PETAL * 1000 + DEVICE)"),
        Column(name="STATE", length=nrows, dtype=np.uint32,
               description="State bit field (good == 0)"),
        Column(name="EXCLUSION", length=nrows, dtype=np.dtype("a9"),
               description="The exclusion polygon for this device"),
    ]

    out_state = Table()
    out_state.add_columns(out_cols)

    row = 0
    for petal in allpetals:
        devlist = list(sorted(fp[petal].keys()))
        for dev in devlist:
            out_state[row]["TIME"] = file_date
            out_state[row]["PETAL"] = fp[petal][dev]["PETAL"]
            out_state[row]["LOCATION"] = fp[petal][dev]["PETAL"] * 1000 + dev
            out_state[row]["DEVICE"] = dev
            out_state[row]["STATE"] = 0
            out_state[row]["EXCLUSION"] = exclusion
            row += 1

    out_state.write(out_state_file, format="ascii.ecsv")

    # Now write out the exclusion polygons.  Since these are not tabular, we
    # write to a YAML file.

    with open(out_poly_file, "w") as pf:
        yaml.dump(poly, pf, default_flow_style=False)

    return
