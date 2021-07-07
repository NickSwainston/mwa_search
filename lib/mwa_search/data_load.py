"""
Loads all the data required by mwa_search from the data directory.
"""

import os

datadir = os.path.join(os.path.dirname(__file__), 'data')

# Hard code the path of the MWA tile receiver temperature file
SMART_POWER_FILE = os.path.join(datadir, 'SMART_obs_data.npy')