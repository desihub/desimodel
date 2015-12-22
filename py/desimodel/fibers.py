# License information goes here
# -*- coding: utf-8 -*-
"""
====================
desimodel.fibers
====================

Provides fibers class
"""

from astropy.io import fits
import numpy as np
import os


class Fibers(object):
    """
    Keeps the relevant information to position fibers on the focal plane
    and its relationship to possible targets.

    Attributes:
        The properties initialized in the __init__ procedure:
        x_focal (float) : array for the x_positions on the focal plane, in mm
        y_focal (float) : array for the y_positions on the focal plane, in mm        
        z_focal (float) : array for the y_positions on the focal plane, in mm        
        fiber_id (int) :
        positioner_id (int) : 
        spectrograph_id (int) : 
        neighbors (int) : 2D array of shape (n_fibers, 6) holding the fiber of the 6 nearest fibers.
        n_fiber (int) : total number of fibers
    """

    def __init__(self):
        filename = os.getenv('DESIMODEL')+'/data/focalplane/fiberpos.fits'

        hdulist = fits.open(filename)        
        self.filename = filename
        self.x_focal = hdulist[1].data['x']
        self.y_focal = hdulist[1].data['y']
        self.z_focal = hdulist[1].data['z']
        self.fiber_id = hdulist[1].data['fiber']
        self.positioner_id = hdulist[1].data['positioner']
        self.spectrograph_id = hdulist[1].data['spectrograph']
        self.neighbors = np.zeros((np.size(self.x_focal), 6), dtype='i4') 
        self.n_fiber = np.size(self.x_focal)

        for i in range(self.n_fiber):
            x = self.x_focal[i]
            y = self.y_focal[i]
            radius = np.sqrt((self.x_focal -x )** 2 + (self.y_focal - y)**2)
            ids = radius.argsort()
            self.neighbors[i,:] = ids[1:7]
        

        # This section is related to targets
        self.available_targets = [None] * self.n_fiber
        self.n_targets = np.zeros(self.n_fiber)
        self.target = -1 * np.ones(self.n_fiber, dtype=np.int64)

    def set_available(self, position, ID_list):
        """
        Set-up the list and number of available targets to this fiber[position]
         
        Args:
             position (int): position in the list to be updated
             ID_list (int): array of available galaxies
             
        """
        self.available_targets[position] = ID_list.copy()
        self.n_targets[position] = np.size(self.available_targets)

    def reset_available(self, position):
        """
        Reset the list and number of available targets to this positioner
         
        Args:
             position (int): position in the list to be updated
        """
        self.available_targets[position] = None
        self.n_targets[position] = 0

    def reset_all_available(self):
        """
        Resets the list and number of available targets to this positioner
        """
        self.available_targets = [None] * self.n_fiber
        self.n_targets = np.zeros(self.n_fiber)

    def set_target(self, position, target_id):
        """
        Sets the id of the target assigned to this positioner
        Args:
             position (int): position in the list to be updated
            target_id (int): id of the target assigned to this positioner
        """
        self.target[position]  = target_id


    def reset_target(self, position):
        """
        resets the id of the target assigned to this positioner
        Args:
             position (int): position in the list to be updated
            target_id (int): id of the target assigned to this positioner
        """
        self.target[position]  = -1

    def reset_all_targets(self):
        """
        resets the id of the target assigned to this positioner
        """
        self.target = -1 * np.ones(self.n_fiber,  dtype=np.int64)

