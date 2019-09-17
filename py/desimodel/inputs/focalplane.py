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


def create(testdir=None, posdir=None, polyfile=None, fibermaps=None,
           petalloc=None, startvalid=None, fillfake=False, exclusion="legacy",
           fakeoffset=False, fakefiberpos=False):
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

    # Process the fibermap files from DocDB.  We start with this step so that
    # we can build the list of all physical PETAL_IDs that have a mapping.
    # Later, we will get positioner data only for devices that exist in the
    # mapping.
    allpetals = set()
    fp = dict()

    if fakefiberpos:
        # Get the mapping AND device info from the fiberpos instead.
        # We CANNOT use the desimodel "load_fiberpos()" function here, since
        # it seems to strip out the ETC devices.  Manually read the file.
        fiberpos_file = findfile('focalplane/fiberpos-all.fits')
        fpos = Table.read(fiberpos_file)

        # There is no "PETAL_ID" for the fake fiberpos, so we use PETAL
        for pet, dev, devtyp, blk, fib, xoff, yoff in zip(
                fpos["PETAL"], fpos["DEVICE"], fpos["DEVICE_TYPE"],
                fpos["SLITBLOCK"], fpos["BLOCKFIBER"], fpos["X"], fpos["Y"]):
            if (devtyp != "POS") and (devtyp != "ETC"):
                # only consider science positioners and sky monitors.
                continue
            allpetals.add(pet)
            if pet not in fp:
                fp[pet] = dict()
            if dev not in fp[pet]:
                empty = dict()
                empty["DEVICE_ID"] = "FAKE"
                fp[pet][dev] = empty
            fp[pet][dev]["SLITBLOCK"] = blk
            fp[pet][dev]["BLOCKFIBER"] = fib
            fp[pet][dev]["CABLE"] = 0
            fp[pet][dev]["CONDUIT"] = "NA"
            fp[pet][dev]["FWHM"] = 0.0
            fp[pet][dev]["FRD"] = 0.0
            fp[pet][dev]["ABS"] = 0.0
            fp[pet][dev]["DEVICE_ID"] = "FAKE"
            fp[pet][dev]["OFFSET_X"] = xoff
            fp[pet][dev]["OFFSET_Y"] = yoff
            fp[pet][dev]["OFFSET_T"] = 0.0
            fp[pet][dev]["OFFSET_P"] = 0.0
            fp[pet][dev]["MIN_T"] = 0.0
            fp[pet][dev]["MAX_T"] = 380.0
            fp[pet][dev]["MIN_P"] = 0.0
            fp[pet][dev]["MAX_P"] = 200.0
            fp[pet][dev]["LENGTH_R1"] = 3.0
            fp[pet][dev]["LENGTH_R2"] = 3.0
    else:
        for docnum, docver, docname in fibermaps:
            fmfile = None
            try:
                fmfile = docdb.download(docnum, docver, docname)
            except IOError:
                msg = "Could not download {}".format(docname)
                log.error(msg)
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
                        allpetals.add(pet)
                        dev = int(row[cols["DEVICE_LOC"]])
                        cable = int(row[cols["Cable_ID"]])
                        conduit = row[cols["Conduit"]]
                        fwhm = float(row[cols["FWHM@f/3.9"]])
                        fthrough = float(row[cols["FRD_Throughput"]])
                        athrough = float(row[cols["Abs_Throuhgput"]])
                        blkfib = row[cols["slit_position"]].split(":")
                        blk = int(blkfib[0])
                        fib = int(blkfib[1])
                        if pet not in fp:
                            fp[pet] = dict()
                        if dev not in fp[pet]:
                            empty = dict()
                            empty["DEVICE_ID"] = "NONE"
                            fp[pet][dev] = empty
                        fp[pet][dev]["SLITBLOCK"] = blk
                        fp[pet][dev]["BLOCKFIBER"] = fib
                        fp[pet][dev]["CABLE"] = cable
                        fp[pet][dev]["CONDUIT"] = conduit
                        fp[pet][dev]["FWHM"] = fwhm
                        fp[pet][dev]["FRD"] = fthrough
                        fp[pet][dev]["ABS"] = athrough

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

    # Extract and organize the information we are tracking

    for devid, props in pos.items():
        pet = int(props["PETAL_ID"])
        dev = int(props["DEVICE_LOC"])
        if dev not in fp[pet]:
            # This device location must not be a POS type.  We still want
            # its information though.
            fp[pet][dev] = dict()
            fp[pet][dev]["SLITBLOCK"] = -1
            fp[pet][dev]["BLOCKFIBER"] = -1
            fp[pet][dev]["CABLE"] = -1
            fp[pet][dev]["CONDUIT"] = "NA"
            fp[pet][dev]["FWHM"] = 0.0
            fp[pet][dev]["FRD"] = 0.0
            fp[pet][dev]["ABS"] = 0.0
        fp[pet][dev]["DEVICE_ID"] = devid
        t_min, t_max, p_min, p_max = _compute_theta_phi_range(
            props["PHYSICAL_RANGE_T"], props["PHYSICAL_RANGE_P"])
        # These values are incorrect in many cases.  Use nominal values
        # instead (see below).
        # fp[pet][dev]["OFFSET_X"] = props["OFFSET_X"]
        # fp[pet][dev]["OFFSET_Y"] = props["OFFSET_Y"]
        fp[pet][dev]["OFFSET_X"] = 0.0
        fp[pet][dev]["OFFSET_Y"] = 0.0
        fp[pet][dev]["OFFSET_T"] = props["OFFSET_T"]
        fp[pet][dev]["OFFSET_P"] = props["OFFSET_P"]
        fp[pet][dev]["LENGTH_R1"] = props["LENGTH_R1"]
        fp[pet][dev]["LENGTH_R2"] = props["LENGTH_R2"]
        fp[pet][dev]["MIN_T"] = t_min
        fp[pet][dev]["MAX_T"] = t_max
        fp[pet][dev]["MIN_P"] = p_min
        fp[pet][dev]["MAX_P"] = p_max

    # Use DocDB 0530 to get the mapping of device location to type.  Also use
    # this for X/Y offsets in case we are simulating positioners.
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

    for petal in sorted(allpetals):
        devlist = list(sorted(fp[petal].keys()))
        for dev in devlist:
            fp[petal][dev]["DEVICE_TYPE"] = devtype[dev]

    # If we have a map of petal ID to petal location, then use it.  Otherwise
    # just assign them in petal ID order.
    if petalloc is None:
        loc = 0
        for petal in sorted(allpetals):
            devlist = list(sorted(fp[petal].keys()))
            for dev in devlist:
                fp[petal][dev]["PETAL"] = loc
            loc += 1
    else:
        for petal in allpetals:
            devlist = list(sorted(fp[petal].keys()))
            for dev in devlist:
                fp[petal][dev]["PETAL"] = petalloc[petal]

    # If the fillfake option is specified, fill all device locations of POS and
    # ETC type with a fake positioner if there is none there.

    if fillfake:
        t_min, t_max, p_min, p_max = _compute_theta_phi_range(380.0, 200.0)
        for petal in sorted(allpetals):
            devlist = list(sorted(fp[petal].keys()))
            for dev in devlist:
                # The petal location of this petal ID
                petal_loc = fp[petal][dev]["PETAL"]
                if fp[petal][dev]["DEVICE_ID"] == "NONE":
                    fp[petal][dev]["OFFSET_X"] = 0.0
                    fp[petal][dev]["OFFSET_Y"] = 0.0
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
                    fp[petal][dev]["LENGTH_R1"] = 3.0
                    fp[petal][dev]["LENGTH_R2"] = 3.0

    # The pos_settings files seem to have incorrect / unreliable
    # X / Y offsets.  Override these with the nominal values for now
    # and investigate this more when moving to the DB.
    for petal in sorted(allpetals):
        devlist = list(sorted(fp[petal].keys()))
        for dev in devlist:
            if dev in dev_nominal_xy:
                xnom, ynom = dev_nominal_xy[dev]
                fp[petal][dev]["OFFSET_X"] = xnom
                fp[petal][dev]["OFFSET_Y"] = ynom

    # Now rotate the X / Y offsets based on the petal location.
    for petal in sorted(allpetals):
        devlist = list(sorted(fp[petal].keys()))
        for dev in devlist:
            # The petal location of this petal ID
            petal_loc = fp[petal][dev]["PETAL"]
            # Petal 0 is at the "bottom"; See DESI-0530.  Here we
            # rotate the nominal X/Y positions on petal 0 to the
            # petal being considered.
            phi = np.radians((float(7 + petal_loc) * 36.0) % 360.0)
            x = fp[petal][dev]["OFFSET_X"]
            y = fp[petal][dev]["OFFSET_Y"]
            fp[petal][dev]["OFFSET_X"] = np.cos(phi) * x - np.sin(phi) * y
            fp[petal][dev]["OFFSET_Y"] = np.sin(phi) * x + np.cos(phi) * y
            # print(
            #     "Petal {}, location {}, phi {}, ({}, {}) --> ({}, {})"
            #     .format(
            #         petal, petal_loc, phi, x, y, fp[petal][dev]["OFFSET_X"],
            #         fp[petal][dev]["OFFSET_Y"]
            #     ),
            #     flush=True
            # )

    # Now load the file(s) with the exclusion polygons
    # Add the legacy polygons to the dictionary for reference.
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

    def collision_to_segments(raw):
        rx = raw[:, 0]
        ry = raw[:, 1]
        sg = [[float(x), float(y)] for x, y in zip(rx, ry)]
        start = list(sg[0])
        sg.append(start)
        return [sg]

    if polyfile is not None:
        # Add shapes from other files.  The convention used here for
        exprops = configobj.ConfigObj(polyfile, unrepr=True)
        poly["default"] = dict()
        ktheta_raw = np.transpose(np.array(exprops["KEEPOUT_THETA"]))
        poly["default"]["theta"] = dict()
        poly["default"]["theta"]["segments"] = \
            collision_to_segments(ktheta_raw)
        poly["default"]["theta"]["circles"] = list()
        kphi_raw = np.transpose(np.array(exprops["KEEPOUT_PHI"]))
        poly["default"]["phi"] = dict()
        poly["default"]["phi"]["segments"] = collision_to_segments(kphi_raw)
        poly["default"]["phi"]["circles"] = list()
        kpetal_raw = np.transpose(np.array(exprops["KEEPOUT_PTL"]))
        poly["default"]["petal"] = dict()
        poly["default"]["petal"]["segments"] = \
            collision_to_segments(kpetal_raw)
        poly["default"]["petal"]["circles"] = list()
        kgfa_raw = np.transpose(np.array(exprops["KEEPOUT_GFA"]))
        poly["default"]["gfa"] = dict()
        poly["default"]["gfa"]["segments"] = \
            collision_to_segments(kgfa_raw)
        poly["default"]["gfa"]["circles"] = list()

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
    for petal in sorted(allpetals):
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
    for petal in sorted(allpetals):
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

    out_fp.write(out_fp_file, format='ascii.ecsv')
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
    for petal in sorted(allpetals):
        devlist = list(sorted(fp[petal].keys()))
        for dev in devlist:
            out_state[row]["TIME"] = file_date
            out_state[row]["PETAL"] = fp[petal][dev]["PETAL"]
            out_state[row]["LOCATION"] = fp[petal][dev]["PETAL"] * 1000 + dev
            out_state[row]["DEVICE"] = dev
            out_state[row]["STATE"] = 0
            out_state[row]["EXCLUSION"] = exclusion
            row += 1

    out_state.write(out_state_file, format='ascii.ecsv')

    # Now write out the exclusion polygons.  Since these are not tabular, we
    # write to a YAML file.

    with open(out_poly_file, "w") as pf:
        yaml.dump(poly, pf, default_flow_style=False)

    return
