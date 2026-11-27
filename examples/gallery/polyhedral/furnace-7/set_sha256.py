"""Copy a duct partition and set its source_sha256 to the hash of an STL.

    python3 set_sha256.py SURFACE.stl PARTITION_IN.json PARTITION_OUT.json
"""

import hashlib
import json
import sys

stl, src, out = sys.argv[1:4]
partition = json.load(open(src))
partition["source_sha256"] = hashlib.sha256(open(stl, "rb").read()).hexdigest()
with open(out, "w") as f:
    json.dump(partition, f, indent=2)
    f.write("\n")
