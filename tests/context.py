import os
import sys

#
# Provide local path to the ore_algebra package if not installed globally
# 

# sys.path.insert(0, "..")

#
# Provide local path to the volumes package if not installed globally
#

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../src/')))

import exact_convex_volumes