#!/bin/bash
set -e
./khoca.py
exec "$@"
