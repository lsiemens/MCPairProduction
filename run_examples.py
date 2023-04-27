#!/usr/bin/env python3
"""Run examples and tests

This script runs examples and tests for validating mcpairproduction
against the code for the first part of the project.
"""

import mcpairproduction.examples as examples

if __name__ == "__main__":
    examples.use_totcross.muon_Z_pole()

    examples.reproduce_totcross.compare_energy_distributions()
    examples.reproduce_totcross.compare_angular_distributions()
