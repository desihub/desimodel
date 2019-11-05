.. code-block:: console

    usage: desi_update_focalplane_exclusion [-h]
                                            [--exclusion [EXCLUSION [EXCLUSION ...]]]
                                            [--time TIME]

    optional arguments:
      -h, --help            show this help message and exit
      --exclusion [EXCLUSION [EXCLUSION ...]]
                            One or more text config files containing some
                            exclusion polygons to use. For example,
                            '$DESIMODEL/exclusions.conf'. This file should contain
                            a parameter 'NAME' with the name to use for this set
                            of exclusions.
      --time TIME           Optional date/time (default is current date/time) to
                            use when selecting the focalplane.Format is YYYY-MM-
                            DDTHH:mm:ss in UTC time.
