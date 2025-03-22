####################################################################################
# TODO:
#    Utilize the ability to run multiple jobs in parallel.
#    Perhaps continue the golden-section search for both possible sections (and for
#    each subsection continue the search in the two new possible subsubsections and
#    so forth) and collaps child sections if the parent section resolved to not be
#    a minimum, or once a minimum is found.
####################################################################################

# This tunes the value of the range seperating parameter, omega (in units of 1/1000 bohr^-1)
# for long range corrected hybrid functionals for a given system, by enforcing DFT's
# version of Koopman's theorem, IP = -HOMO. This is done via a golden-section search
# algorithm to minimize the value of (IP + E_homo)^2.
# Note that this requires the premise that (IP + E_homo)^2 is unimodal with respect
# to omega for the range of omegas under investigation.

# An input file for the system we are optimizing named N.in must be in the
# same directory as this script. Additionally, the input for a +1 cation of
# the system with the same functional and basis (if the original system's charge
# is +1 then this cation will be +2 and so forth) should be named P.in and must
# be in the aforementioned directory. For both input files, replace the value for
# OMEGA in the $rem section with the string 'xxx'.
# The qchem.sh file should be in the same directory as this script as well.

# Example for N.in for H20 and the LRC-wPBEh functional:
"""
$molecule
    0 1
    O        -0.0589468554    0.2678084446    0.0000000000
    H        -0.9158854857   -0.1614307441    0.0000000000
    H         0.5742936224   -0.4515745405    0.0000000000
$end

$rem
    jobtype     sp
    exchange    gen
    lrc_dft     true
    omega       xxx
    basis       def2-TZVPP
$end

$xc_functional
    C   PBE     1.0
    X   wPBE    0.8
    K           0.2
$end
"""

# Example for P.in for H2O and the LRC-wPBEh functional:
"""
$molecule
    1 2
    O        -0.0589468554    0.2678084446    0.0000000000
    H        -0.9158854857   -0.1614307441    0.0000000000
    H         0.5742936224   -0.4515745405    0.0000000000
$end

$rem
    jobtype     SP
    exchange    gen
    lrc_dft     TRUE
    omega       xxx
    basis       def2-TZVPP
$end

$xc_functional
    C   PBE     1.0
    X   wPBE    0.8
    K           0.2
$end
"""

import datetime
import math
import os
import subprocess
import sys
from argparse import ArgumentParser
from typing import Callable

OUTPUT_FILE = "tuning.out"


def write(text: str):
    """
    Append text to the output file.
    """
    with open(OUTPUT_FILE, "a+") as output:
        output.write(text)


def gss(f: Callable[[int], float], a: float, b: float, tolerance: float) -> int:
    """
    Golden-Section search - find the minimum of f on [a,b]
    * f must be strictly unimodal on [a,b]
    This implementation is for a function that takes in only intigers (f:Z->R).
    """

    INVPHI = (math.sqrt(5) - 1) / 2  # 1/phi where phi is the golden ratio

    itteration = 1
    while b - a > tolerance:
        c = b - (b - a) * INVPHI
        d = a + (b - a) * INVPHI

        write(
            f"\n----------------------------------------------------\n"
            f"ITTERATION {itteration}\n"
            f"GSS values: a:{a}, b:{b}, c:{c}, d:{d}\n"
        )

        if f(round(c)) < f(round(d)):
            b = d
        else:
            a = c
        itteration += 1

    return round((b + a) / 2)


def IP_plus_HOMO_squared(omega: int) -> float:
    """
    Return the sqaure of the sum of the IP and the homo energies for a given omega.
    """
    write(f"\nStarting Calculations for omega={omega}:\n")

    # create new input files with omega
    if os.path.isfile(f"scratch/input/N_{omega}.in") and os.path.isfile(
        f"scratch/input/P_{omega}.in"
    ):
        write("Input files already exist.\n")
    else:
        write("Creating input files...\n")
        create_input_files(omega)

    # run qchem jobs for both N and P
    if os.path.isfile(f"scratch/output/N_{omega}.out") and os.path.isfile(
        f"scratch/output/P_{omega}.out"
    ):
        write("Output files already exist.\n")
    else:
        write("Running calculations...\n")
        run_qchem_job(f"scratch/input/N_{omega}.in", f"scratch/output/N_{omega}.out")
        run_qchem_job(f"scratch/input/P_{omega}.in", f"scratch/output/P_{omega}.out")
        write("Calculations done.\n")

    # extract data from the output
    write("\nExtracting data from output...\n")
    E_N, E_homo = extract_energies(f"scratch/output/N_{omega}.out")
    E_P, _ = extract_energies(f"scratch/output/P_{omega}.out")
    write(f"E(N): {E_N}, E(N-1): {E_P}, HOMO: {E_homo}\n")

    IP = E_P - E_N
    squared_sum = (IP + E_homo) ** 2
    write(f"(IP + HOMO)^2 = {squared_sum}\n")

    return squared_sum


def create_input_files(omega: int) -> None:
    """
    Create input files for the N and P calculations with the given omega.
    """
    with open(f"N.in", "r") as N:
        N_input = N.read()

    N_input = N_input.replace("xxx", f"{omega}")

    with open(f"scratch/input/N_{omega}.in", "w") as N_omega:
        N_omega.write(N_input)

    with open(f"P.in", "r") as P:
        P_input = P.read()

    P_input = P_input.replace("xxx", f"{omega}")

    with open(f"scratch/input/P_{omega}.in", "w") as P_omega:
        P_omega.write(P_input)


def run_qchem_job(input: str, output: str) -> None:
    """
    Run a qchem calculation, this blocks the script untill the job has finished.
    """
    subprocess.call(["bash", "qchem.sh", f"{input}", f"{output}"])


def extract_energies(file: str) -> tuple[float, float]:
    """
    Extract Total energy and homo energy from a qchem single-point calculation output file.
    """
    # initializtion and flags
    E_total = 0
    E_homo = 0
    restricted = False
    passed_alpha = False
    passed_beta = False
    symmetry_on_MO = False

    # read the output file
    with open(os.path.abspath(file), "r") as N_output:
        lines = N_output.readlines()
        for i, line in enumerate(lines):
            if "Orbital Energies (a.u.) and Symmetries" in line:
                symmetry_on_MO = True
            elif "Total energy in the final basis set" in line:
                E_total = float(line.split()[-1])
            elif "A restricted SCF calculation" in line:
                restricted = True
            elif "Alpha MOs" in line:
                passed_alpha = True
            elif "Beta MOs" in line:
                passed_beta = True
            elif "-- Virtual --" in line:
                if passed_alpha and not passed_beta:
                    if symmetry_on_MO:
                        E_homo = float(lines[i - 2].split()[-1])
                    else:
                        E_homo = float(lines[i - 1].split()[-1])
                    if restricted:
                        break
                elif passed_alpha and passed_beta:
                    if symmetry_on_MO:
                        E_beta = float(lines[i - 2].split()[-1])
                    else:
                        E_beta = float(lines[i - 1].split()[-1])
                    if E_beta > E_homo:
                        E_homo = E_beta
                    break
                else:
                    write("ERROR: could not find alpha orbital energies")

    # If E_total or E_homo remain zero (as they were initialized) we must have not found them in the ouput.
    if E_total == 0:
        write("ERROR: E_total could not be found\n")
        raise Exception("ERROR: E_total could not be found")
    if E_homo == 0:
        write("ERROR: E_homo could not be found\n")
        raise Exception("ERROR: E_homo could not be found")

    return (E_total, E_homo)


def main(argv: list[str]) -> int:
    # parse arguments
    parser = ArgumentParser(
        prog="tuning.py",
        description="Tune the range seperating parmater, omega, for long range corrected hybrid functionals.",
    )
    parser.add_argument(
        "-b",
        "--bounds",
        nargs=2,
        action="store",
        default=[200, 800],
        type=float,
        dest="bounds",
        help="Custom bounds for omega, must be two intigers between 0 and 1000.",
    )
    parser.add_argument(
        "-t",
        "--tolerance",
        action="store",
        default=0.2,
        type=float,
        dest="tolerance",
        help="When the size of range containing omega reaches the tolerance, rounding the center of the range is given as omega.",
    )
    args = parser.parse_args(argv[1:])

    # extract bounds for omega from the arguments
    if args.bounds[0] < args.bounds[1]:
        lower_bound, upper_bound = args.bounds
    else:
        upper_bound, lower_bound = args.bounds
    if lower_bound < 0 or upper_bound > 1000:
        parser.error("Bounds for omega must be between 0 and 1000.")

    # extract convergance tolerance from the arguments
    tolerance = args.tolerance
    print(tolerance, upper_bound - lower_bound)
    if tolerance > upper_bound - lower_bound:
        parser.error("Tolerance must smaller than the starting bounds of omega.")

    # create directories for files created by this script under the "scratch" directory
    if not os.path.exists("scratch/input"):
        os.makedirs("scratch/input")
    if not os.path.exists("scratch/output"):
        os.makedirs("scratch/output")

    # create an output file, it will be appended with each calculation.
    with open(OUTPUT_FILE, "w") as output:
        output.write(
            f"Output file for omega tuning.\n"
            f"\nJob Specifications:\n"
            f"-------------------\n"
            f"Bounds for omega: [{lower_bound}, {upper_bound}]\n"
            f"Tolerance: {tolerance}\n"
        )
        # copy N.in and P.in to the output file
        with open("N.in", "r") as N:
            output.write(
                f"\n------------------------ N.in ----------------------\n"
                f"{N.read()}"
                f"----------------------------------------------------\n"
            )
        with open("P.in", "r") as P:
            output.write(
                f"\n------------------------ P.in ----------------------\n"
                f"{P.read()}"
                f"----------------------------------------------------\n"
            )

    # run the optimization
    optimal_omega = gss(IP_plus_HOMO_squared, lower_bound, upper_bound, tolerance)

    # write the final value for omega
    write(
        f"\n----------------------------------------------------\n"
        f"|                   CONVERGED                      |\n"
        f"----------------------------------------------------\n"
        f"\nTuned omega: {optimal_omega}\n"
    )

    # move scratch directory into a archive directory so that it doesn't interfere with the next job
    write("\nArchiving files...\n")
    if not os.path.exists("archive"):
        os.makedirs("archive")
    timestamp = str(datetime.datetime.now()).replace(" ", "_").replace(":", "-")
    subprocess.call(["mv", "scratch", f"archive/{timestamp}"])

    write("\nFinished.")
    return 0


if __name__ == "__main__":
    main(sys.argv)
