import os

# Path to the mispr/Gaussian results database config file (db.json).
# Set this via an environment variable rather than hardcoding a path:
# export EBATT_DB_FILE=/path/to/your/db.json
DB_FILE = os.environ.get("EBATT_DB_FILE")
