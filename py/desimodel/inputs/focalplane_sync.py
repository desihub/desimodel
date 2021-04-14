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
    fpcal = Table.read(path, format="ascii.ecsv")
    return fpcal


# Rotation matrices.  These few lines of code are copied from desimeter,
# since it is the only thing that is needed from that package and there
# is nothing even desi-specific here.


def Rx(angle):  # all in radians
    Rx = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, np.cos(angle), -np.sin(angle)],
            [0.0, np.sin(angle), np.cos(angle)],
        ]
    )
    return Rx


def Ry(angle):  # all in radians
    Ry = np.array(
        [
            [np.cos(angle), 0.0, np.sin(angle)],
            [0.0, 1.0, 0.0],
            [-np.sin(angle), 0.0, np.cos(angle)],
        ]
    )
    return Ry


def Rz(angle):  # all in radians
    Rz = np.array(
        [
            [np.cos(angle), -np.sin(angle), 0.0],
            [np.sin(angle), np.cos(angle), 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    return Rz


def Rxyz(alpha, beta, gamma):  # yaw-pitch-roll system, all in radians
    return Rz(gamma) @ Ry(beta) @ Rx(alpha)  # @ is matrix multiplication


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
    cal_time = datetime.datetime.strptime(cal_time_str, "%Y-%m-%dT%H:%M:%S%z")

    state_time_str = cal_time.isoformat(timespec="seconds")

    # Parse other metadata

    eo_phi = fpcal.meta["Eo_phi"]
    eo_radius = fpcal.meta["Eo_radius_with_margin"]
    alignments = fpcal.meta["PETAL_ALIGNMENTS"]

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

    # Build the transformation information for each petal alignment
    petal_rot = dict()
    petal_trans = dict()
    for petal_id, align in alignments.items():
        petal_rot[petal_id] = Rxyz(align["alpha"], align["beta"], align["gamma"])
        petal_trans[petal_id] = np.array([align["Tx"], align["Ty"], align["Tz"]])

    n_rows = len(fpcal)

    fp, state = create_tables(n_rows)

    # We only want to track the POS_P and POS_T values for positioners which are
    # non-functional (stuck) or have a broken fiber, since these are the positioners
    # we cannot move.
    stuck_or_broken = valid_states["BROKEN"] | valid_states["STUCK"]

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

        # The OFFSET_T angle is in petal coordinates.  We express this in global
        # coordinates by transforming 2 points and computing the angle from that.

        prot = petal_rot[d["PETAL_ID"]]
        ptrans = petal_trans[d["PETAL_ID"]]
        offset_pt_rad = np.radians(d["OFFSET_T"])
        offset_pt_pnts = np.array(
            [
                [d["OFFSET_X_CS5"], d["OFFSET_Y_CS5"], 0.0],
                [
                    3.0 * np.cos(offset_pt_rad) + d["OFFSET_X_CS5"],
                    3.0 * np.sin(offset_pt_rad) + d["OFFSET_Y_CS5"],
                    0.0,
                ],
            ]
        )
        offset_fp_pnts = prot.dot(offset_pt_pnts.T).T + np.tile(ptrans, 2).reshape(
            (2, -1)
        )
        offset_fp_rad = np.arctan2(
            offset_fp_pnts[1, 1] - offset_fp_pnts[0, 1],
            offset_fp_pnts[1, 0] - offset_fp_pnts[0, 0],
        )
        offset_fp_deg = np.degrees(offset_fp_rad)
        if offset_fp_deg < -180.0:
            offset_fp_deg += 360.0
        if offset_fp_deg > 180.0:
            offset_fp_deg -= 360.0

        # # Sanity check that the alignment-transformed offset is "close" to the
        # # expected value.
        # petalrot_check = (float(7 + fp["PETAL"][r]) * 36.0) % 360.0
        # petalrot_check += d["OFFSET_T"]
        # if petalrot_check < -180.0:
        #     petalrot_check += 360.0
        # if petalrot_check > 180.0:
        #     petalrot_check -= 360.0
        # print(
        #     "device {}, petal {}, offset_t = {}, check = {}".format(
        #         fp["DEVICE"][r], fp["PETAL"][r], offset_fp_deg, petalrot_check
        #     ),
        #     flush=True,
        # )

        fp["OFFSET_T"][r] = offset_fp_deg

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

        # Starting point is "good" with nominal MIN_P
        state["STATE"][r] = valid_states["OK"]
        state["MIN_P"][r] = d["MIN_P"]

        # Even if we are simulating, we want to mark both broken fibers and
        # non-movable positioners.
        if not d["FIBER_INTACT"]:
            state["STATE"][r] |= valid_states["BROKEN"]
        if d["DEVICE_CLASSIFIED_NONFUNCTIONAL"]:
            state["STATE"][r] |= valid_states["STUCK"]

        if not sim:
            # We want to also check for retracted positioners
            if d["CLASSIFIED_AS_RETRACTED"]:
                # This positioner is retracted.  Set the exclusion to the retracted
                # one and also limit the phi angle range.
                state["STATE"][r] |= valid_states["RESTRICT"]
                state["EXCLUSION"][r] = "retracted"
                # The focalplane cal information defines the Eo_Phi angle to be the
                # minimum Phi angle relative to the coordinate axis, not the offset.
                # So to get MIN_P we must subtract the offset.
                state["MIN_P"][r] = eo_phi - d["OFFSET_P"]

        # The other positioner angles in the state are just the same as nominal
        state["MAX_P"][r] = d["MAX_P"]
        state["MIN_T"][r] = d["MIN_T"]
        state["MAX_T"][r] = d["MAX_T"]
        # If the device is not movable, track its current estimated location.
        if state["STATE"][r] & stuck_or_broken:
            state["POS_P"][r] = d["POS_P"]
            state["POS_T"][r] = d["POS_T"]
        else:
            state["POS_P"][r] = 0.0
            state["POS_T"][r] = 0.0

    return (fp, state, excl, state_time_str)


def create_from_calibs(
    calib_file,
    out_dir=None,
    reset=False,
    sim_good=False,
    commit=False,
    test=False,
    fibermaps=None,
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
        test (bool):  If True, perform all operations but do not update any files.
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

        log.info("  Output focalplane:  %s", out_fp_file)
        log.info("  Output state     :  %s", out_state_file)
        log.info("  Output exclusions:  %s", out_excl_file)

        if test:
            log.info("Running with test==True, skipping file writes.")
            return

        fp.write(out_fp_file, format="ascii.ecsv", overwrite=True)

        state.write(out_state_file, format="ascii.ecsv", overwrite=True)

        with gzip.open(out_excl_file, "wb") as pf:
            yaml.dump(excl, stream=pf, encoding="utf-8", default_flow_style=False)

        if commit:
            cmesg = "Creating new focalplane model from DB sync {}".format(date_str)
            try:
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
            except sp.CalledProcessError:
                log.error("svn update / commit returned an error")
    else:
        # Load the current focalplane and just update the state

        # Get the focalplane from one second before the current datestamp
        cur_date = datetime.datetime.strptime(date_str, "%Y-%m-%dT%H:%M:%S%z")
        dtime = datetime.timedelta(seconds=1)
        oldtime = cur_date - dtime

        oldfp, oldexcl, oldstate, oldtmstr = load_focalplane(oldtime)

        log.info("Comparing generated focalplane to one from %s", oldtmstr)

        # Compare the old and new.
        checkcols = set(fp.colnames)
        diff = device_compare(oldfp, fp, list(checkcols))

        if len(diff) > 0:
            device_printdiff(diff)
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

        # If needed, promote the "TIME" state column to the new column type
        if len(st["TIME"][0]) < 21:
            newt = Column(
                name="TIME",
                length=len(st["TIME"]),
                dtype=np.dtype("a30"),
                description="The timestamp of the event (UTC, ISO format)",
            )
            for row, oldt in enumerate(st["TIME"]):
                nt = datetime.datetime.strptime(oldt, "%Y-%m-%dT%H:%M:%S")
                nt = nt.replace(tzinfo=datetime.timezone.utc)
                newt[row] = nt.isoformat()
            st.replace_column("TIME", newt)

        excl_file = os.path.join(out_dir, "desi-exclusion_{}.yaml.gz".format(oldtmstr))

        tmp_state = "{}.tmp".format(state_file)
        prev_state = "{}.previous".format(state_file)
        tmp_excl = "{}.tmp".format(excl_file)
        prev_excl = "{}.previous".format(excl_file)

        device_printdiff(state_diff)

        n_new_states = 0
        n_new_excl = 0
        for loc, df in state_diff.items():
            if df["old"] is None or df["new"] is None:
                # This should never happen, since it means that the LOCATION
                # value does not exist in either the previous or current
                # state.  We already checked for that above.
                msg = "LOCATION {} missing from old or new state.  Should never happen!".format(
                    loc
                )
                raise RuntimeError(msg)

            row = [df["new"][col] for col in df["new"].dtype.names]
            st.add_row(row)
            n_new_states += 1

            new_excl = str(
                df["new"]["EXCLUSION"].tobytes().rstrip(b"\x00"), encoding="utf-8"
            )
            if new_excl not in oldexcl:
                oldexcl[new_excl] = excl[new_excl]
                n_new_excl += 1

        log.info("Updating focalplane:  %s", oldtmstr)
        log.info("  State log appending %d rows", n_new_states)
        if n_new_excl > 0:
            log.info("  Adding %d new exclusion shapes", n_new_excl)
        else:
            log.info("  No new exclusion shapes, not updating file")

        if test:
            log.info("Running with test==True, skipping file writes.")
            return

        # Write to temp file then move into place
        st.write(tmp_state, format="ascii.ecsv", overwrite=True)
        shutil.copy2(state_file, prev_state)
        os.rename(tmp_state, state_file)

        # If we updated any exclusions, write a new file
        if n_new_excl > 0:
            with gzip.open(temp_excl, "wb") as pf:
                yaml.dump(oldexcl, pf, default_flow_style=False)
            shutil.copy2(excl_file, prev_excl)
            os.rename(tmp_excl, excl_file)

        if commit:
            cmesg = "Appending DB sync {} to focalplane model {}".format(
                date_str, oldtmstr
            )
            try:
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
            except sp.CalledProcessError:
                log.error("svn update / commit returned an error")
