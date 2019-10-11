.. code-block:: console

    usage: desi_generate_focalplane [-h] [--pos_settings POS_SETTINGS]
                                    [--collision COLLISION]
                                    [--exclusion EXCLUSION]
                                    [--petal_id2loc PETAL_ID2LOC]
                                    [--startvalid STARTVALID] [--fillfake]
                                    [--test TEST] [--fakeoffset] [--fakefiberpos]

    optional arguments:
      -h, --help            show this help message and exit
      --pos_settings POS_SETTINGS
                            The directory containing the positioner settings
                            files. For example, DESI svn fp_settings/pos_settings.
                            The file names are assumed to have the form
                            unit_<name>.conf. Only devices with an assigned
                            PETAL_ID and DEVICE_LOC are considered.
      --collision COLLISION
                            The text config file containing the exclusion polygons
                            to use. For example,
                            '_collision_settings_DEFAULT.conf'
      --exclusion EXCLUSION
                            The name of the default exclusion polygon to use (e.g.
                            'legacy', 'default', etc)
      --petal_id2loc PETAL_ID2LOC
                            Mapping of petal ID to focalplane location, specified
                            as '<id0>:<loc0>,<id1>:<loc1>,...'. If not specified,
                            petals will be placed in ID order starting at location
                            0.
      --startvalid STARTVALID
                            Optional start date (default is current date/time)
                            when the focalplane becomes valid. Format is YYYY-MM-
                            DDTHH:mm:ss in UTC time.
      --fillfake            For simulations, fill empty device locations with fake
                            positioners
      --test TEST           Override the output directory for testing.
      --fakeoffset          Set theta / phi offsets to zero for fake positioners.
                            This is *only* for consistency tests with legacy
                            fiberassignment and should not be used for anything
                            else.
      --fakefiberpos        Use the old fiberpos file for the device mapping. Only
                            for testing.
