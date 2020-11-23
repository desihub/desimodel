# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
desimodel.inputs.focalplane_sync
===================================

Tools for checking and synchronizing focalplane state to online system.
"""
import os
import datetime
import shutil

import ast

import yaml

import numpy as np

from astropy.table import Table, Column

from desiutil.log import get_logger

from ..io import datadir, findfile, load_focalplane

from .focalplane_utils import (
    create_tables,
    device_loc_to_type,
    load_petal_fiber_map,
    collision_to_segments,
    exclusions_equal,
    valid_states,
    device_compare,
    device_printdiff,
    update_exclusions,
    restricted_positioner_phi,
    hash_exclusion,
)


def load_fp_calibs(path):
    """Load data dumped from the online system."""
    fpcal = Table.read(path)
    return fpcal


def convert_fp_calibs(fpcal):
    """Convert the online system information.

    This returns a tuple containing the focalplane, the current state, and the set of
    unique exclusion polygons found in the file.

    Args:
        fpcal (Table):  The table loaded from a dump from the online system.

    Returns:
        (tuple):  The (focalplane, state, exclusions, time string) loaded from
            the cal file.

    """
    # Parse the calibration time
    cal_time_str = fpcal.meta["DATE_RETRIEVED"]
    cal_time_raw = datetime.datetime.strptime(cal_time_str, "%Y-%m-%dT%H:%M:%S%z")
    cal_time = cal_time_raw.replace(tzinfo=None)

    state_time_str = cal_time.isoformat(timespec="seconds")

    # Get the default exclusion polygons for theta, phi, GFA, and petal
    # boundaries

    excl = dict()

    ktheta_str = fpcal.meta["general_keepout_T"]

    ktheta_raw = np.transpose(
        np.array(ast.literal_eval(fpcal.meta["general_keepout_T"]))
    )
    kphi_raw = np.transpose(np.array(ast.literal_eval(fpcal.meta["general_keepout_P"])))
    kpetal_raw = np.transpose(np.array(ast.literal_eval(fpcal.meta["keepout_PTL"])))
    kgfa_raw = np.transpose(np.array(ast.literal_eval(fpcal.meta["keepout_GFA"])))

    kp = dict()
    kp["theta"] = dict()
    kp["theta"]["segments"] = collision_to_segments(ktheta_raw)
    kp["theta"]["circles"] = list()

    kp["phi"] = dict()
    kp["phi"]["segments"] = collision_to_segments(kphi_raw)
    kp["phi"]["circles"] = list()

    kp["petal"] = dict()
    kp["petal"]["segments"] = collision_to_segments(kpetal_raw)
    kp["petal"]["circles"] = list()

    kp["gfa"] = dict()
    kp["gfa"]["segments"] = collision_to_segments(kgfa_raw)
    kp["gfa"]["circles"] = list()

    excl["default"] = kp

    # Also make a set of exclusions for "retracted" positioners.  In this case,
    # the PHI polygon is set to a circle of 2.1mm, which will extend to the
    # edge of the nominal THETA keepout polygon.  This edge is the "outer clear
    # rotation envelope".

    outer_clear_rotation = 2.1  # mm

    retrct = dict(excl["default"])
    retrct["phi"]["segments"] = list()
    retrct["phi"]["circles"] = [[[0.0, 0.0], outer_clear_rotation]]
    excl["retracted"] = retrct

    # Get the fiber map from device location to spectrographs.
    fmap = load_petal_fiber_map()

    n_rows = len(fpcal)

    fp, state = create_tables(n_rows)

    kindx = 0

    for r in range(n_rows):
        d = fpcal[r]
        # First set the focalplane properties
        fp["PETAL"][r] = d["PETAL_LOC"]
        fp["PETAL_ID"][r] = d["PETAL_ID"]
        fp["DEVICE"][r] = d["DEVICE_LOC"]
        fp["DEVICE_ID"][r] = d["POS_ID"]
        fp["LOCATION"][r] = d["PETAL_LOC"] * 1000 + d["DEVICE_LOC"]
        fp["DEVICE_TYPE"][r] = device_loc_to_type(fp["DEVICE"][r])
        fp["SLITBLOCK"][r] = -1
        fp["BLOCKFIBER"][r] = -1
        fp["CABLE"][r] = -1
        fp["CONDUIT"][r] = "NA"
        fp["FWHM"][r] = 0.0
        fp["FRD"][r] = 0.0
        fp["ABS"][r] = 0.0
        if d["PETAL_ID"] in fmap:
            # We have some information about the device to fiber mapping.
            fp["SLITBLOCK"][r] = fmap[d["PETAL_ID"]][d["DEVICE_LOC"]]["SLITBLOCK"]
            fp["BLOCKFIBER"][r] = fmap[d["PETAL_ID"]][d["DEVICE_LOC"]]["BLOCKFIBER"]
            fp["CABLE"][r] = fmap[d["PETAL_ID"]][d["DEVICE_LOC"]]["CABLE"]
            fp["CONDUIT"][r] = fmap[d["PETAL_ID"]][d["DEVICE_LOC"]]["CONDUIT"]
            fp["FWHM"][r] = fmap[d["PETAL_ID"]][d["DEVICE_LOC"]]["FWHM"]
            fp["FRD"][r] = fmap[d["PETAL_ID"]][d["DEVICE_LOC"]]["FRD"]
            fp["ABS"][r] = fmap[d["PETAL_ID"]][d["DEVICE_LOC"]]["ABS"]
        fp["LENGTH_R1"][r] = d["LENGTH_R1"]
        fp["LENGTH_R2"][r] = d["LENGTH_R2"]
        fp["OFFSET_X"][r] = d["OFFSET_X_CS5"]
        fp["OFFSET_Y"][r] = d["OFFSET_Y_CS5"]
        fp["MIN_P"][r] = d["MIN_P"]
        fp["MAX_P"][r] = d["MAX_P"]
        fp["MIN_T"][r] = d["MIN_T"]
        fp["MAX_T"][r] = d["MAX_T"]
        fp["OFFSET_P"][r] = d["OFFSET_P"]

        # FIXME:  This is in focal surface coordinates, but is petal-local.  Need
        # to rotate to correct petal location.  We should use desimeter for this
        # conversion instead of the ideal case here.
        petalrot_deg = (float(7 + fp["PETAL"][r]) * 36.0) % 360.0
        fp["OFFSET_T"][r] = d["OFFSET_T"] + petalrot_deg

        # Set the FIBER
        fp["FIBER"][r] = -1
        if fp["SLITBLOCK"][r] >= 0:
            fp["FIBER"][r] = (
                fp["PETAL"][r] * 500 + fp["SLITBLOCK"][r] * 25 + fp["BLOCKFIBER"][r]
            )

        # Handle the keepout polygons for this device.  Every row specifies a keepout
        # even if many are close to identical.  Here we build up a dictionary of unique
        # keepouts and find the one to use for this positioner.

        dktheta = np.transpose(np.array(ast.literal_eval(d["KEEPOUT_T"])))
        dpolytheta = dict()
        dpolytheta["circles"] = list()
        dpolytheta["segments"] = collision_to_segments(dktheta)

        dkphi = np.transpose(np.array(ast.literal_eval(d["KEEPOUT_P"])))
        dpolyphi = dict()
        dpolyphi["circles"] = list()
        dpolyphi["segments"] = collision_to_segments(dkphi)

        kp = dict(excl["default"])
        kp["theta"] = dpolytheta
        kp["phi"] = dpolyphi
        pname = hash_exclusion(kp)[:16]

        if pname not in excl:
            excl[pname] = kp

        state["EXCLUSION"][r] = pname

        # Now set rest of the state table

        state["TIME"][r] = state_time_str
        state["LOCATION"][r] = d["PETAL_LOC"] * 1000 + d["DEVICE_LOC"]
        state["STATE"][r] = valid_states["OK"]
        if d["DEVICE_CLASSIFIED_NONFUNCTIONAL"]:
            state["STATE"][r] |= valid_states["STUCK"]
        if not d["FIBER_INTACT"]:
            state["STATE"][r] |= valid_states["BROKEN"]
        if d["CLASSIFIED_AS_RETRACTED"]:
            # This positioner is retracted.  Set the exclusion to the retracted
            # one and also limit the phi angle range.
            state["STATE"][r] |= valid_states["RESTRICT"]
            state["EXCLUSION"][r] = "retracted"
            state["MIN_P"] = restricted_positioner_phi(
                outer_clear_rotation + fp["LENGTH_R1"][r],
                fp["LENGTH_R1"][r],
                fp["LENGTH_R2"][r],
                fp["OFFSET_P"][r],
                fp["MIN_P"][r],
                fp["MAX_P"][r],
            )
        else:
            state["MIN_P"] = d["MIN_P"]
        # If the device is NOT good, track its current estimated location.
        if state["STATE"][r] == valid_states["OK"]:
            state["POS_P"][r] = 0.0
            state["POS_T"][r] = 0.0
        else:
            state["POS_P"][r] = d["POS_P"]
            state["POS_T"][r] = d["POS_T"]

    return (fp, state, excl, state_time_str)


def create_from_calibs(
    calib_file, testdir=None, fibermaps=None, force=False, sim=False
):
    """Construct a DESI focalplane from a calibration dump.

    This uses a dump from the online system and compares it to the current
    latest focalplane model and state.

    Args:
        calib_file (str):  Path to the calibration dump.
        testdir (str):  Override the output directory for testing.
        fibermaps (list):  Override list of tuples (DocDB number,
            DocDB version, DocDB csv file) of where to find the petal mapping
            files.
        force (bool):  If True, ignore the current focalplane model and
            create a new model from this calibration dump.  Default compares
            the new focalplane to the old and looks for changes in device
            state.  These changes are appended to the existing log.
        sim (bool):  If True, clear all transient state issues and set hardware
            to be as "good as possible", for use in simulations.

    Returns:
        None

    """
    log = get_logger()

    outdir = testdir
    if outdir is None:
        outdir = os.path.join(datadir(), "focalplane")

    # Get the model from the calib file
    log.info("Loading calibration dump from %s ...", calib_file)
    fpcal = load_fp_calibs(calib_file)

    log.info("Converting calibration format ...")
    fp, state, excl, date_str = convert_fp_calibs(fpcal)

    log.info("Calibration data retrieval date = %s", date_str)

    if force:
        # Ignore any previous focalplane info and dump out what we have

        log.info("Writing new focalplane model ...")
        out_fp_file = os.path.join(outdir, "desi-focalplane_{}.ecsv".format(date_str))
        out_excl_file = os.path.join(outdir, "desi-exclusion_{}.yaml".format(date_str))
        out_state_file = os.path.join(outdir, "desi-state_{}.ecsv".format(date_str))

        fp.write(out_fp_file, format="ascii.ecsv", overwrite=True)

        state.write(out_state_file, format="ascii.ecsv", overwrite=True)

        with open(out_excl_file, "w") as pf:
            yaml.dump(excl, pf, default_flow_style=False)

    else:
        # Load the current focalplane and just update the state

        # Get the focalplane from one second before the current datestamp
        cur_date = datetime.datetime.strptime(date_str, "%Y-%m-%dT%H:%M:%S")
        dtime = datetime.timedelta(seconds=1)
        oldtime = cur_date - dtime

        oldfp, oldexcl, oldstate, oldtmstr = load_focalplane(oldtime)

        log.info("Comparing generated focalplane to one from %s", oldtmstr)

        # Compare the old and new.
        checkcols = set(fp.colnames)
        checkcols.remove("LOCATION")
        diff = device_compare(oldfp, fp, list(checkcols))

        device_printdiff(diff)

        if len(diff) > 0:
            msg = (
                "Existing focalplane device properties have changed."
                "  Use the 'force' option to start with a new focalplane model."
            )
            raise RuntimeError(msg)

        # We got this far, which means that the focalplane data agrees.  Now
        # Look for differences in the state.
        checkcols = [
            "STATE",
            "MIN_P",
            "POS_P",
            "POS_T",
            "EXCLUSION",
        ]
        state_diff = device_compare(oldstate, state, checkcols)

        if len(state_diff) > 0:
            # We have some changes, load the full table and append.
            state_file = os.path.join(outdir, "desi-state_{}.ecsv".format(oldtmstr))
            st = Table.read(state_file, format="ascii.ecsv")

            excl_file = os.path.join(outdir, "desi-exclusion_{}.ecsv".format(oldtmstr))

            tmp_state = "{}.tmp".format(state_file)
            prev_state = "{}.previous".format(state_file)
            tmp_excl = "{}.tmp".format(excl_file)
            prev_excl = "{}.previous".format(excl_file)

            need_excl_update = False

            for loc, df in state_diff.items():
                if df["old"] is None or df["new"] is None:
                    # This should never happen, since it means that the LOCATION
                    # value does not exist in either the previous or current
                    # state.  We already checked for that above.
                    msg = "LOCATION {} missing from old or new state.  Should never happen!".format(
                        loc
                    )
                    raise RuntimeError(msg)

                new_st = df["new"]["STATE"]
                new_excl = str(
                    df["new"]["EXCLUSION"].tobytes().rstrip(b"\x00"), encoding="utf-8"
                )
                new_minp = df["new"]["MIN_P"]
                new_posp = df["new"]["POS_P"]
                new_post = df["new"]["POS_T"]

                msg = "Updating state for location {}:".format(loc)
                msg += "\n  old: state = {}, POS_T = {}, POS_P = {}, MIN_P = {}, excl = {}".format(
                    df["old"]["STATE"],
                    df["old"]["POS_T"],
                    df["old"]["POS_P"],
                    df["old"]["MIN_P"],
                    str(
                        df["old"]["EXCLUSION"].tobytes().rstrip(b"\x00"),
                        encoding="utf-8",
                    ),
                )
                msg += "\n  new: state = {}, POS_T = {}, POS_P = {}, MIN_P = {}, excl = {}".format(
                    new_st, new_post, new_posp, new_minp, new_excl
                )
                log.info(msg)

                st.add_row(
                    [date_str, loc, new_st, new_post, new_posp, new_minp, new_excl]
                )

                if new_excl not in oldexcl:
                    oldexcl[new_excl] = excl[new_excl]
                    need_excl_update = True

            # Write to temp file then move into place
            st.write(tmp_state, format="ascii.ecsv", overwrite=True)
            shutil.copy2(state_file, prev_state)
            os.rename(tmp_state, state_file)

            # If we updated any exclusions, write a new file
            if need_excl_update:
                with open(temp_excl, "w") as pf:
                    yaml.dump(oldexcl, pf, default_flow_style=False)
                shutil.copy2(excl_file, prev_excl)
                os.rename(tmp_excl, excl_file)
