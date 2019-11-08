# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
desimodel.inputs.focalplane_utils
=====================================

Helpers for constructing a focalplane model.
"""
import configobj
import numpy as np

from desiutil.log import get_logger

from . import docdb


valid_states = {
    "OK": 0,
    "STUCK": 2,
    "BROKEN": 4,
}


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


def _compare_direction(out, first, fkey, second, skey):
    """Do one direction of the comparison.

    This updates "out" with the things in "second" that do not exist in
    "first".

    """
    newpet = list(sorted(second.keys()))
    for petal in newpet:
        newdev = list(sorted(second[petal].keys()))
        if petal not in first:
            # This is a whole new petal...
            if petal not in out:
                out[petal] = dict()
            for dev in newdev:
                props = list(sorted(second[petal][dev].keys()))
                if dev not in out[petal]:
                    out[petal][dev] = dict()
                for p in props:
                    if p not in out[petal][dev]:
                        out[petal][dev][p] = dict()
                    out[petal][dev][p][fkey] = None
                    out[petal][dev][p][skey] = second[petal][dev][p]
        else:
            for dev in newdev:
                props = list(sorted(second[petal][dev].keys()))
                if dev not in first[petal]:
                    # This device is missing
                    if petal not in out:
                        out[petal] = dict()
                    if dev not in out[petal]:
                        out[petal][dev] = dict()
                    for p in props:
                        if p not in out[petal][dev]:
                            out[petal][dev][p] = dict()
                        out[petal][dev][p][fkey] = None
                        out[petal][dev][p][skey] = second[petal][dev][p]
                else:
                    for p in props:
                        if (p not in first[petal][dev]) or (
                            first[petal][dev][p] != second[petal][dev][p]
                        ):
                            # This property is missing or mismatched
                            if petal not in out:
                                out[petal] = dict()
                            if dev not in out[petal]:
                                out[petal][dev] = dict()
                            if p not in out[petal][dev]:
                                out[petal][dev][p] = dict()
                            if p in first[petal][dev]:
                                out[petal][dev][p][fkey] = first[petal][dev]
                                out[petal][dev][p][skey] = \
                                    second[petal][dev][p]
                            else:
                                out[petal][dev][p][fkey] = None
                                out[petal][dev][p][skey] = \
                                    second[petal][dev][p]
    return


def compare(fpold, fpnew):
    """Compare two sets of focalplane device properties.

    Args:
        fpold (dict):  The original device properties.
        fpnew (dict):  The new device properties.

    Returns:
        (dict):  A dictionary with the same structure as the inputs, but
            only differing values are included.  One additional nested
            level is added to the "leaves" of the hierarchy and this
            contains the keys "old" and "new".

    """
    out = dict()
    _compare_direction(out, fpold, "old", fpnew, "new")
    _compare_direction(out, fpnew, "new", fpold, "old")
    return out


def collision_to_segments(raw):
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
