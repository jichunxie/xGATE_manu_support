#!/usr/bin/env python3
"""verify_xgate.py
Try importing the xGATE package under common module names.
"""
import sys

try:
    import xgate as x
    print("xGATE imported as 'xgate' module")
    sys.exit(0)
except Exception:
    pass

try:
    import xGATE as x
    print("xGATE imported as 'xGATE' module")
    sys.exit(0)
except Exception:
    pass

print("xGATE import failed — ensure 'pip install git+https://github.com/jichunxie/xGATE.git' succeeded.")
sys.exit(2)
