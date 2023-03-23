"""Utilities and helper functions
"""

import sys

def progress_bar(message, step, max, divisions=10):
    progress = int(divisions*(step + 1)/max)
    bar = "#"*progress + " "*(divisions - progress)
    sys.stdout.write(f"\r{message}: [{bar}] {step + 1}/{max}")
    sys.stdout.flush()
