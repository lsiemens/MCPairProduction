"""Utilities and helper functions
"""

import sys
import numpy

def progress_bar(message, step, max, divisions=10):
    progress = int(divisions*(step + 1)/max)
    bar = "#"*progress + " "*(divisions - progress)
    sys.stdout.write(f"\r{message}: [{bar}] {step + 1}/{max}")
    sys.stdout.flush()

def nice_ticks(i, num, symbol="\\pi"):
    num, i = num//numpy.gcd(num, i), i//numpy.gcd(num, i)
    if (i == 0):
        return "$0$"
    if (num == 1) and (i == 1):
        return f"${symbol}$"
    if (num == 1):
        return f"${i} {symbol}$"
    if (i == 1):
        return f"${symbol}/{num}$"
    return f"${i} {symbol}/{num}$"
