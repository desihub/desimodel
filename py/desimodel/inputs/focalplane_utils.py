# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
desimodel.inputs.focalplane_utils
=====================================

Helpers for constructing a focalplane model.
"""
import hashlib
import configobj
import csv
import numpy as np

from astropy.table import Table, Column

from desiutil.log import get_logger

from . import docdb


valid_states = {
    "OK": 0,
    "STUCK": 2,
    "BROKEN": 4,
    "RESTRICT": 8,
}


def device_loc_to_type(loc):
    """Get the fixed, hardcoded device type for a device location.
    """
    if loc in [461, 501]:
        return "ETC"
    elif loc in [541, 542]:
        return "GIF"
    elif loc in [11, 75, 150, 239, 321, 439, 482, 496, 517, 534]:
        return "FIF"
    elif loc in [
        38, 331, 438, 460, 478, 479, 480, 481, 497, 498, 499, 500, 513, 514,
        515, 516, 527, 528, 529, 530, 531, 535, 536, 537, 538, 539, 540
    ]:
        return "NON"
    else:
        return "POS"


def compute_theta_phi_range(phys_t, phys_p):
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


def restricted_positioner_phi(radius, theta_arm, phi_arm, offset_p, min_p, max_p):
    """Compute the MIN_P angle needed to restrict positioner reach.

    Given the positioner arm lengths, desired maximum reach, and PHI offset
    and min / max, compute the new minimum phi angle needed to keep the
    positioner within the retracted radius.

    Args:
        radius (float):  The restricted radius in mm.
        theta_arm (float):  The theta arm length in mm.
        phi_arm (float):  The phi arm length in mm.
        offset_p (float):  The OFFSET_P phi zero point.
        min_p (float):  The MIN_P minimum phi angle.
        max_p (float):  The MAX_P maximum phi angle.

    Returns:
        (float):  The restricted MIN_P value.

    """
    phi_zero = np.radians(offset_p)
    phi_min = np.radians(min_p)
    phi_max = np.radians(max_p)
    # Use law of cosines to find max opening angle
    opening = np.degrees(
        np.arccos(
            (radius**2 - theta_arm**2 - phi_arm**2) /
            (-2.0 * theta_arm * phi_arm)
        )
    )
    # Phi min is relative to the offset
    return (180.0 - opening) - offset_p


def create_device():
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


def rotate_petals(fp):
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


def load_petal_fiber_map(existing=None, fibermaps=None):
    """Loads info from petal verification files.

    This loads the petal verification files from DocDB and populates a dictionary
    of properties.

    Args:
        existing (dict):  An existing dictionary to populate.  If None, a new one is
            created and returned.
        fibermaps (list):  (optional) Override list of tuples (DocDB number,
            DocDB version, DocDB csv file) of where to find the petal mapping
            files.

    Returns:
        (dict):  A dictionary of dictionaries with the device location info
            for every petal ID.

    """
    fp = existing
    if existing is None:
        fp = dict()
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
                    if existing is None:
                        if (pet not in fp):
                            fp[pet] = dict()
                        if (dev not in fp[pet]):
                            fp[pet][dev] = dict()
                    else:
                        if (pet not in fp) or (dev not in fp[pet]):
                            print("FAIL:  petal {}, dev {} not in fp".format(
                                pet, dev
                            ), flush=True)
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
    if 261 not in fp[6]:
        fp[6][261] = dict()
    fp[6][261]["DEVICE_TYPE"] = "POS"
    fp[6][261]["DEVICE_ID"] = "NONE"  # Populated later from pos_settings
    fp[6][261]["SLITBLOCK"] = 19
    fp[6][261]["BLOCKFIBER"] = 22
    fp[6][261]["CABLE"] = 6
    fp[6][261]["CONDUIT"] = "E0"      # This conduit has one fewer than F3
    fp[6][261]["FWHM"] = 0.0          # No information from file
    fp[6][261]["FRD"] = 0.0           # No information from file
    fp[6][261]["ABS"] = 0.0           # No information from file
    # ---------------------------
    # DESI-4883v4-Petal_11_final_verification.csv
    # Petal ID 11 is missing an entry for device location 484.  This
    # device location is assigned to positioner M06847 in the pos_settings
    # files.  Assign it to the one missing fiber location.
    log.info("Correcting petal ID 11, location 484")
    if 484 not in fp[11]:
        fp[11][484] = dict()
    fp[11][484]["DEVICE_TYPE"] = "POS"
    fp[11][484]["DEVICE_ID"] = "NONE"  # Populated below from pos_settings
    fp[11][484]["SLITBLOCK"] = 3
    fp[11][484]["BLOCKFIBER"] = 3
    fp[11][484]["CABLE"] = 4
    fp[11][484]["CONDUIT"] = "G0"     # This conduit has one fewer than G1
    fp[11][484]["FWHM"] = 0.0         # No information from file
    fp[11][484]["FRD"] = 0.0          # No information from file
    fp[11][484]["ABS"] = 0.0          # No information from file
    return fp


def create_nominal(petal_loc):
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
            props = create_device()
            props["PETAL"] = petal_loc[petal]
            props["PETAL_ID"] = petal
            props["DEVICE_TYPE"] = devtype[loc]
            x, y = dev_nominal_xy[loc]
            props["OFFSET_X"] = x
            props["OFFSET_Y"] = y
            pt[loc] = props
        fp[petal] = pt
    return fp


def device_compare(fpold, fpnew, check):
    """Compare two sets of focalplane device properties.

    Args:
        fpold (Table):  The original device properties.
        fpnew (Table):  The new device properties.
        check (list):  The column names to check for equality.

    Returns:
        (dict):  A dictionary containing the differences.  The keys are
            the LOCATION value, and the value is a dict with "old" and "new"
            keys that contain the table rows that differ.

    """
    out = dict()
    olddiff = np.setdiff1d(fpold[:]["LOCATION"], fpnew[:]["LOCATION"])
    rows = np.arange(len(fpold), dtype=np.int)[fpold[:]["LOCATION"] in olddiff]
    for r in rows:
        loc = fpold[r]["LOCATION"]
        out[loc] = dict()
        out[loc]["old"] = fpold[r]
        out[loc]["new"] = None
    totdiff = set(olddiff)
    newdiff = np.setdiff1d(fpnew[:]["LOCATION"], fpold[:]["LOCATION"])
    rows = np.arange(len(fpnew), dtype=np.int)[fpnew[:]["LOCATION"] in newdiff]
    for r in rows:
        loc = fnew[r]["LOCATION"]
        out[loc] = dict()
        out[loc]["new"] = fpnew[r]
        out[loc]["old"] = None
    totdiff.update(newdiff)

    # Go through locations found in both tables and look for differences in
    # the column values.

    old_loc_to_row = dict()
    for indx, row in enumerate(fpold):
        old_loc_to_row[row["LOCATION"]] = indx

    for row in fpnew:
        loc = row["LOCATION"]
        if loc in totdiff:
            continue
        oldrow = fpold[old_loc_to_row[loc]]
        for col in check:
            if row[col] != oldrow[col]:
                out[loc] = dict()
                out[loc]["old"] = np.copy(oldrow)
                out[loc]["new"] = np.copy(row)
                break
    return out


def device_printdiff(diff):
    """Print a diff dictionary created with device_compare().
    """
    for loc, df in diff.items():
        print("Location {:04d}:".format(loc))
        ol = df["old"]
        nw = df["new"]
        if ol is None:
            print("  OLD:  None")
        else:
            print("  OLD:")
            for col in ol.dtype.names:
                print("    {}:  {}".format(col, ol[col]))
        if nw is None:
            print("  NEW:  None")
        else:
            print("  NEW:")
            for col in nw.dtype.names:
                print("    {}:  {}".format(col, nw[col]))
    print("", flush=True)
    return


def collision_to_segments(raw):
    rx = raw[:, 0]
    ry = raw[:, 1]
    sg = [[float(x), float(y)] for x, y in zip(rx, ry)]
    start = list(sg[0])
    if not np.allclose(sg[0], sg[-1]):
        sg.append(start)
    return [sg]


def exclusions_equal(ex1, ex2):
    """Return True if the two polygons are equal, else False"""
    if len(ex1["segments"]) != len(ex2["segments"]):
        return False
    if len(ex1["circles"]) != len(ex2["circles"]):
        return False
    for slist1, slist2 in zip(ex1["segments"], ex2["segments"]):
        if len(slist1) != len(slist2):
            return False
        for s1, s2 in zip(slist1, slist2):
            if not np.allclose(s1, s2):
                return False
    for c1, c2 in zip(ex1["circles"], ex2["circles"]):
        if not np.allclose(c1, c2):
            return False
    return True


def hash_exclusion(excl):
    exhash = hashlib.md5()
    polynames = list(sorted(excl.keys()))
    for nm in polynames:
        for seglist in excl[nm]["segments"]:
            for seg in seglist:
                segstr = "{:0.4f}{:0.4f}".format(seg[0], seg[1])
                exhash.update(segstr.encode("utf-8"))
        for cir in excl[nm]["circles"]:
            cent = cir[0]
            rad = cir[1]
            cirstr = "{:0.4f}{:0.4f}{:0.4f}".format(cent[0], cent[1], rad)
            exhash.update(cirstr.encode("utf-8"))
    return exhash.hexdigest()


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
        log.info("Loading exclusion polygons from %s", pf)
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
            collision_to_segments(ktheta_raw)
        props["theta"]["circles"] = list()
        kphi_raw = np.transpose(np.array(exprops["KEEPOUT_PHI"]))
        props["phi"] = dict()
        props["phi"]["segments"] = collision_to_segments(kphi_raw)
        props["phi"]["circles"] = list()
        kpetal_raw = np.transpose(np.array(exprops["KEEPOUT_PTL"]))
        props["petal"] = dict()
        props["petal"]["segments"] = \
            collision_to_segments(kpetal_raw)
        props["petal"]["circles"] = list()
        kgfa_raw = np.transpose(np.array(exprops["KEEPOUT_GFA"]))
        props["gfa"] = dict()
        props["gfa"]["segments"] = \
            collision_to_segments(kgfa_raw)
        props["gfa"]["circles"] = list()
        excl[nm] = props
    return


def create_tables(n_fp_rows, n_state_rows=None):
    """Create empty focalplane and state tables.

    This function keeps the construction of the table schema in a single
    place that can be used across the code.

    Args:
        n_fp_rows (int):  The number of rows in the focalplane table.
        n_state_rows (int):  The number of rows in the state table.  If None, use
            the same number of rows as the focalplane table.

    Returns:
        (tuple):  The (focalplane, state) tables.

    """
    if n_fp_rows is None or n_fp_rows < 1:
        raise ValueError("number of focaplane table rows must be an integer > 0")
    if n_state_rows is None:
        n_state_rows = n_fp_rows

    fp_cols = [
        Column(name="PETAL", length=n_fp_rows, dtype=np.int32,
               description="Petal location [0-9]"),
        Column(name="DEVICE", length=n_fp_rows, dtype=np.int32,
               description="Device location on the petal"),
        Column(name="LOCATION", length=n_fp_rows, dtype=np.int32,
               description="PETAL * 1000 + DEVICE"),
        Column(name="PETAL_ID", length=n_fp_rows, dtype=np.int32,
               description="The physical petal ID"),
        Column(name="DEVICE_ID", length=n_fp_rows, dtype=np.dtype("a9"),
               description="The physical device ID string"),
        Column(name="DEVICE_TYPE", length=n_fp_rows, dtype=np.dtype("a3"),
               description="The device type (POS, ETC, FIF)"),
        Column(name="SLITBLOCK", length=n_fp_rows, dtype=np.int32,
               description="The slit block where this fiber goes"),
        Column(name="BLOCKFIBER", length=n_fp_rows, dtype=np.int32,
               description="The fiber index within the slit block"),
        Column(name="CABLE", length=n_fp_rows, dtype=np.int32,
               description="The cable ID"),
        Column(name="CONDUIT", length=n_fp_rows, dtype=np.dtype("a3"),
               description="The conduit"),
        Column(name="FIBER", length=n_fp_rows, dtype=np.int32,
               description="PETAL * 500 + SLITBLOCK * 25 + BLOCKFIBER"),
        Column(name="FWHM", length=n_fp_rows, dtype=np.float64,
               description="FWHM at f/3.9"),
        Column(name="FRD", length=n_fp_rows, dtype=np.float64,
               description="FRD Throughput"),
        Column(name="ABS", length=n_fp_rows, dtype=np.float64,
               description="ABS Throughput"),
        Column(name="OFFSET_X", length=n_fp_rows, dtype=np.float64,
               description="X location of positioner center", unit="mm"),
        Column(name="OFFSET_Y", length=n_fp_rows, dtype=np.float64,
               description="Y location of positioner center", unit="mm"),
        Column(name="OFFSET_T", length=n_fp_rows, dtype=np.float64,
               description="THETA zero point angle", unit="degrees"),
        Column(name="OFFSET_P", length=n_fp_rows, dtype=np.float64,
               description="PHI zero point angle", unit="degrees"),
        Column(name="LENGTH_R1", length=n_fp_rows, dtype=np.float64,
               description="Length of THETA arm", unit="mm"),
        Column(name="LENGTH_R2", length=n_fp_rows, dtype=np.float64,
               description="Length of PHI arm", unit="mm"),
        Column(name="MAX_T", length=n_fp_rows, dtype=np.float64,
               description="Maximum THETA angle relative to OFFSET_T",
               unit="degrees"),
        Column(name="MIN_T", length=n_fp_rows, dtype=np.float64,
               description="Minimum THETA angle relative to OFFSET_T",
               unit="degrees"),
        Column(name="MAX_P", length=n_fp_rows, dtype=np.float64,
               description="Maximum PHI angle relative to OFFSET_P",
               unit="degrees"),
        Column(name="MIN_P", length=n_fp_rows, dtype=np.float64,
               description="Minimum PHI angle relative to OFFSET_P",
               unit="degrees"),
    ]

    fp = Table()
    fp.add_columns(fp_cols)

    state_cols = [
        Column(name="TIME", length=n_state_rows, dtype=np.dtype("a20"),
               description="The timestamp of the event (UTC, ISO format)"),
        Column(name="LOCATION", length=n_state_rows, dtype=np.int32,
               description="Global device location (PETAL * 1000 + DEVICE)"),
        Column(name="STATE", length=n_state_rows, dtype=np.uint32,
               description="State bit field (good == 0)"),
        Column(name="POS_T", length=n_state_rows, dtype=np.float32,
               description="Current estimate of Theta arm angle"),
        Column(name="POS_P", length=n_state_rows, dtype=np.float32,
               description="Current estimate of Phi arm angle"),
        Column(name="MIN_P", length=n_state_rows, dtype=np.float32,
               description="Current minimum Phi angle (restricted reach)"),
        Column(name="EXCLUSION", length=n_state_rows, dtype=np.dtype("a16"),
               description="The exclusion polygon for this device"),
    ]

    state = Table()
    state.add_columns(state_cols)
    return (fp, state)


def propagate_state(state, excl, oldstate, oldexcl):
    """Propagate state to a new focalplane model.

    This takes a new state and exclusions and sets this new state to be the
    same as the old one for all locations that are specified in the oldstate.
    Any exclusions not defined in the new dictionary are copied from the old.

    This function assumes that the old and new focalplane models have
    already been verified to be identical...

    Args:
        state (Table):  The new state table, modified in place.
        excl (dict):  The new exclusions, modified in place.
        oldstate (Table):  The old state table.
        oldexcl (dict):  The old exclusions.

    Returns:
        None

    """
    old_row = {
        y: x for x, y in enumerate(oldstate["LOCATION"])
    }
    copyexcl = set()
    for r in range(len(state)):
        loc = state[r]["LOCATION"]
        if loc in old_row:
            for col in ["STATE", "EXCLUSION", "MIN_P", "POS_P", "POS_T"]:
                state[r][col] = oldstate[old_row[loc]][col]
            copyexcl.add(state[r]["EXCLUSION"])
    for xcopy in copyexcl:
        if xcopy not in excl:
            excl[xcopy] = oldexcl[xcopy]
    return
