#!/bin/bash
# Wrapper script for AF3_summary_results.py that can be run from anywhere

# Resolve the folder where this script itself is located
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Call the Python script inside the same folder, forwarding all arguments
python3 "$SCRIPT_DIR/AF3_summary_results.py" "$@"