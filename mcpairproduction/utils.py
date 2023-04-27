"""Utilities and helper functions
"""

import sys
import numpy

def progress_bar(message, step, max, divisions=10):
    """Update command line progress bar

    Parameters
    ----------
    message : string
        Loading bar message.
    step : integer
        Current progress to display.
    max : integer
        Maximum number of steps.
    divisions : integer
        Number of subdivisions of the loading bar.
    """
    progress = int(divisions*(step + 1)/max)
    bar = "#"*progress + " "*(divisions - progress)
    sys.stdout.write(f"\r{message}: [{bar}] {step + 1}/{max}")
    sys.stdout.flush()

def nice_ticks(numerator, denominator, symbol="\\pi"):
    """Get LaTeX labels for ticks

    Parameters
    ----------
    numerator : integer
        Numerator in the fractional graph tick.
    denominator : integer
        Denominator in the fractional graph tick.
    symbol : string
        Latex math symbol of common factor.

    Returns
    -------
    string
        Tick label as simplified fraction.
    """
    denominator, numerator = (denominator//numpy.gcd(denominator, numerator),
                              numerator//numpy.gcd(denominator, numerator))
    if (numerator == 0):
        return "$0$"
    if (denominator == 1) and (numerator == 1):
        return f"${symbol}$"
    if (denominator == 1):
        return f"${numerator} {symbol}$"
    if (numerator == 1):
        return f"${symbol}/{denominator}$"
    return f"${numerator} {symbol}/{denominator}$"
