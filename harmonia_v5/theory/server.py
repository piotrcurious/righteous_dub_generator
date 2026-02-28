#!/usr/bin/env python3
import sys
import os

# Add current directory to path so we can import as a package
sys.path.append(os.path.dirname(os.path.abspath(__file__)) + "/..")
from theory.main import main

if __name__ == "__main__":
    main()
