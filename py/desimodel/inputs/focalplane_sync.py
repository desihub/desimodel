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
import gzip

import subprocess as sp

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


def convert_fp_calibs(fpcal, sim=False):
    """Convert the online system information.

    This returns a tuple containing the focalplane, the current state, and the set of
    unique exclusion polygons found in the file.

    Args:
        fpcal (Table):  The table loaded from a dump from the online system.
        sim (bool):  If True, clear all transient state issues and set hardware
            to be as "good as possible", for use in simulations.

    Returns:
        (tuple):  The (focalplane, state, exclusions, time string) loaded from
            the cal file.

    """
    # Parse the calibration time
    cal_time_str = fpcal.meta["DATE_RETRIEVED"]
    cal_time_raw = datetime.datetime.strptime(cal_time_str, "%Y-%m-%dT%H:%M:%S%z")
    cal_time = cal_time_raw.replace(tzinfo=None)

    state_time_str = cal_time.isoformat(timespec="seconds")

    # Parse other metadata

    eo_phi = fpcal.meta["Eo_phi"]
    eo_radius = fpcal.meta["Eo_radius_with_margin"]

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

    # Also make a set of exclusions for "retracted" positioners.  In this case place
    # a circle at the the center of the theta axis.

    retrct = dict(excl["default"])
    retrct["theta"]["segments"] = list()
    retrct["theta"]["circles"] = [[[0.0, 0.0], eo_radius]]
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
            if d["DEVICE_LOC"] in fmap[d["PETAL_ID"]]:
                # This is a POS or ETC device
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
        # Even if we are simulating, we do want to mark broken fibers
        if not d["FIBER_INTACT"]:
            state["STATE"][r] |= valid_states["BROKEN"]
        if not sim:
            # We want the true state...
            if d["DEVICE_CLASSIFIED_NONFUNCTIONAL"]:
                state["STATE"][r] |= valid_states["STUCK"]
            if d["CLASSIFIED_AS_RETRACTED"]:
                # This positioner is retracted.  Set the exclusion to the retracted
                # one and also limit the phi angle range.
                state["STATE"][r] |= valid_states["RESTRICT"]
                state["EXCLUSION"][r] = "retracted"
                # The focalplane cal information defines the Eo_Phi angle to be the
                # minimum Phi angle relative to the coordinate axis, not the offset.
                # So to get MIN_P we must subtract the offset.
                state["MIN_P"][r] = eo_phi - d["OFFSET_P"]
            else:
                state["MIN_P"][r] = d["MIN_P"]
        # The other positioner angles in the state are just the same as nominal
        state["MAX_P"][r] = d["MAX_P"]
        state["MIN_T"][r] = d["MIN_T"]
        state["MAX_T"][r] = d["MAX_T"]
        # If the device is NOT good, track its current estimated location.
        if state["STATE"][r] == valid_states["OK"]:
            state["POS_P"][r] = 0.0
            state["POS_T"][r] = 0.0
        else:
            state["POS_P"][r] = d["POS_P"]
            state["POS_T"][r] = d["POS_T"]

    return (fp, state, excl, state_time_str)


def create_from_calibs(
    calib_file, out_dir=None, reset=False, sim_good=False, commit=False, fibermaps=None
):
    """Construct a DESI focalplane from a calibration dump.

    This uses a dump from the online system and compares it to the current
    latest focalplane model and state.

    Args:
        calib_file (str):  Path to the calibration dump.
        out_dir (str):  Override the output directory for testing.  Default writes to
            $DESIMODEL/data/focalplane/
        reset (bool):  If True, ignore the current focalplane model and
            create a new model from this calibration dump.  Default compares
            the new focalplane to the old and looks for changes in device
            state.  These changes are appended to the existing log.
        sim_good (bool):  If True, clear all transient state issues and set hardware
            to be as "good as possible", for use in simulations.
        commit (bool):  If True, attempt to commit the result.
        fibermaps (list):  Override list of tuples (DocDB number,
            DocDB version, DocDB csv file) of where to find the petal mapping
            files.

    Returns:
        None

    """
    log = get_logger()

    if out_dir is None:
        out_dir = os.path.join(datadir(), "focalplane")
        test_svn = os.path.join(datadir(), ".svn")
        if not os.path.isdir(test_svn):
            test_svn = os.path.join(os.path.dirname(datadir()), ".svn")
            if not os.path.isdir(test_svn):
                msg = "Output data directory:  {}".format(out_dir)
                msg += "\nis not inside an svn checkout.  You will not be able to"
                msg += "\ncommit these changes without copying them to a checkout."
                log.warning(msg)
    else:
        log.warning("Using debug output directory %s", out_dir)
        log.warning("Files cannot be used until placed in $DESIMODEL/data/focalplane")

    # Get the model from the calib file
    log.info("Loading calibration dump from %s ...", calib_file)
    fpcal = load_fp_calibs(calib_file)

    log.info("Converting calibration format ...")
    fp, state, excl, date_str = convert_fp_calibs(fpcal, sim=sim_good)

    log.info("Calibration data retrieval date = %s", date_str)

    if reset:
        # Ignore any previous focalplane info and dump out what we have

        log.info("Writing new focalplane model- ignoring previous ones...")
        out_fp_file = os.path.join(out_dir, "desi-focalplane_{}.ecsv".format(date_str))
        out_excl_file = os.path.join(
            out_dir, "desi-exclusion_{}.yaml.gz".format(date_str)
        )
        out_state_file = os.path.join(out_dir, "desi-state_{}.ecsv".format(date_str))

        fp.write(out_fp_file, format="ascii.ecsv", overwrite=True)

        state.write(out_state_file, format="ascii.ecsv", overwrite=True)

        with gzip.open(out_excl_file, "wb") as pf:
            yaml.dump(excl, stream=pf, encoding="utf-8", default_flow_style=False)

        if commit:
            cmesg = "Creating new focalplane model from DB sync {}".format(date_str)
            sp.check_call(["svn", "update"], cwd=out_dir)
            sp.check_call(
                [
                    "svn",
                    "add",
                    out_fp_file,
                    out_excl_file,
                    out_state_file,
                ],
                cwd=out_dir,
            )
            sp.check_call(["svn", "commit", "-m", cmesg], cwd=out_dir)
            sp.check_call(["svn", "update"], cwd=out_dir)
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
        diff = device_compare(oldfp, fp, list(checkcols))

        # device_printdiff(diff)

        if len(diff) > 0:
            msg = (
                "Existing focalplane device properties have changed."
                "  Use the 'reset' option to start with a new focalplane model."
            )
            log.error(msg)
            return

        # We got this far, which means that the focalplane data agrees.  Now
        # Look for differences in the state.
        checkcols = set(state.colnames)
        checkcols.remove("TIME")
        state_diff = device_compare(oldstate, state, list(checkcols))

        if len(state_diff) == 0:
            log.info("New focalplane state is identical, no action needed.")
            return

        # We must have some changes, load the full table and append.
        state_file = os.path.join(out_dir, "desi-state_{}.ecsv".format(oldtmstr))
        st = Table.read(state_file, format="ascii.ecsv")

        excl_file = os.path.join(out_dir, "desi-exclusion_{}.yaml.gz".format(oldtmstr))

        tmp_state = "{}.tmp".format(state_file)
        prev_state = "{}.previous".format(state_file)
        tmp_excl = "{}.tmp".format(excl_file)
        prev_excl = "{}.previous".format(excl_file)

        need_excl_update = False

        device_printdiff(state_diff)

        for loc, df in state_diff.items():
            if df["old"] is None or df["new"] is None:
                # This should never happen, since it means that the LOCATION
                # value does not exist in either the previous or current
                # state.  We already checked for that above.
                msg = "LOCATION {} missing from old or new state.  Should never happen!".format(
                    loc
                )
                raise RuntimeError(msg)

            # new_st = df["new"]["STATE"]

            # new_minp = df["new"]["MIN_P"]
            # new_posp = df["new"]["POS_P"]
            # new_post = df["new"]["POS_T"]
            #
            # msg = "Updating state for location {}:".format(loc)
            # msg += "\n  old: state = {}, POS_T = {}, POS_P = {}, MIN_P = {}, excl = {}".format(
            #     df["old"]["STATE"],
            #     df["old"]["POS_T"],
            #     df["old"]["POS_P"],
            #     df["old"]["MIN_P"],
            #     str(
            #         df["old"]["EXCLUSION"].tobytes().rstrip(b"\x00"),
            #         encoding="utf-8",
            #     ),
            # )
            # msg += "\n  new: state = {}, POS_T = {}, POS_P = {}, MIN_P = {}, excl = {}".format(
            #     new_st, new_post, new_posp, new_minp, new_excl
            # )
            # log.info(msg)

            row = [df["new"][col] for col in df["new"].dtype.names]
            st.add_row(row)

            new_excl = str(
                df["new"]["EXCLUSION"].tobytes().rstrip(b"\x00"), encoding="utf-8"
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
            with gzip.open(temp_excl, "wb") as pf:
                yaml.dump(oldexcl, pf, default_flow_style=False)
            shutil.copy2(excl_file, prev_excl)
            os.rename(tmp_excl, excl_file)

        if commit:
            cmesg = "Appending DB sync {} to focalplane model {}".format(
                date_str, oldtmstr
            )
            sp.check_call(["svn", "update"], cwd=out_dir)
            sp.check_call(
                [
                    "svn",
                    "commit",
                    "-m",
                    cmesg,
                ],
                cwd=out_dir,
            )
            sp.check_call(["svn", "update"], cwd=out_dir)
