.. code-block:: console

    usage: desi_update_focalplane_log [-h] [--location LOCATION] [--petal PETAL]
                                      [--device DEVICE] [--state STATE]
                                      [--exclusion EXCLUSION] [--time TIME]

    optional arguments:
      -h, --help            show this help message and exit
      --location LOCATION   The device location (petal * 1000 + device loc)
                            modified by this event.
      --petal PETAL         The petal (focalplane location, not petal ID) modified
                            by this event (--device must also be used)
      --device DEVICE       The device location (--petal must also be given)
                            modified by this event.
      --state STATE         The new state to assign to the device.
      --exclusion EXCLUSION
                            The new exclusion polygon to assign to the device
                            (e.g. 'legacy', 'default', etc)
      --time TIME           Optional date/time (default is current date/time) when
                            this event happens. Format is YYYY-MM-DDTHH:mm:ss in
                            UTC time.
