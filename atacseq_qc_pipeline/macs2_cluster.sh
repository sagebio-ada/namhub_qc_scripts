#!/bin/bash
# MACS2 2.2.9.1 via the cluster module, with its environment scoped to this one command.
# The bioconda macs2 build fails on this system (undefined symbol __log_finite), and
# `module load macs2` in the calling shell would set PYTHONPATH to the cluster's python
# site-packages, which breaks other conda tools (idr fails on a numpy/scipy mismatch).
# Use as: atacseq_idr_concordance.py --macs2 macs2_cluster.sh
type module >/dev/null 2>&1 || source /etc/profile.d/modules.sh
module load macs2/2.2.9
exec macs2 "$@"
