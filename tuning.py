####################################################################################
# TODO:
# 1) Seperate the use of slurm from this script, such that this script does not rely
#    on slurm and only on qchem, then this script can be run on the cluster using
#    slurm if so wished.
# 2) Utilize the ability to run multiple jobs in parallel.
#    Perhaps continue the golden-section search for both possible sections (and for
#    each subsection continue the search in the two new possible subsubsections and
#    so forth) and collaps child sections if the parent section resolved to not be
#    a minimum, or once a minimum is found.
####################################################################################

# This script finds the optimal omega value (in units of 1/1000 bohr^-1) for long
# range corrected hybrid functionals for a given system, by enforcing DFT's version
# of Koopman's theorem, IP = -E_homo. This is done via a golden-section search
# algorithm to minimize the value of (IP + E_homo)^2.
# Note that this requires the premise that (IP + E_homo)^2 is unimodal for
# the range of omegas under investigation.

# An input file for the system we are optimizing named N.in must be in the
# same directory as this script. Additionally, the input for a +1 cation of
# the system (if the original system's charge is +1 then this cation will be
# +2 and so forth) should be named P.in and must be in the aforementioned
# directory. For both input files replace the value for OMEGA in the $rem
# section with the string 'xxx'.
# The single_point.sh file should be in the same directory as this script as well.

# example for N.in using H20 and the LRC-wPBEh functional:
"""
$molecule
    0 1
    O        -0.0589468554    0.2678084446    0.0000000000
    H        -0.9158854857   -0.1614307441    0.0000000000
    H         0.5742936224   -0.4515745405    0.0000000000
$end

$rem
    JOBTYPE     SP
    EXCHANGE    gen
    LRC_DFT     TRUE
    OMEGA       xxx
    BASIS       def2-TZVPP
$end

$xc_functional
    C   PBE     1.0
    X   wPBE    0.8
    K           0.2
$end
"""

# example for P.in using H2O and the LRC-wPBEh functional:
"""
$molecule
    1 2
    O        -0.0589468554    0.2678084446    0.0000000000
    H        -0.9158854857   -0.1614307441    0.0000000000
    H         0.5742936224   -0.4515745405    0.0000000000
$end

$rem
    JOBTYPE     SP
    EXCHANGE    gen
    LRC_DFT     TRUE
    OMEGA       xxx
    BASIS       def2-TZVPP
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


def write(text: str):
    """
    Append text to the output file.
    """
    with open("tuning.out", "a+") as output:
        output.write(text)


# golden-section search algorithm
def gss(f: Callable[[int], float], a: int, b: int, tolerance: int) -> int:
    """
    Golden-Section search - find the minimum of f on [a,b]
    * f must be strictly unimodal on [a,b]
    This implementation uses intigers for x-values since omega must be an intiger between 0 and 1000.
    """

    INVPHI = (math.sqrt(5) - 1) / 2  # 1/phi where phi is the golden ratio

    f_a = f_b = 0
    while b - a > tolerance:
        c = round(b - (b - a) * INVPHI)
        d = round(a + (b - a) * INVPHI)

        write(
            f"\n----------------------------------------------------\n"
            f"GSS values: a:{a}, b:{b}, c:{c}, d:{d}\n"
        )

        f_c = f(c)
        f_d = f(d)
        if f_c < f_d:
            b = d
            f_b = f_d
        if f(c) < f(d):
            b = d
            f_b = f_d
        else:
            a = c
            f_a = f_c

    if f_a < f_b:
        return math.floor((b + a) / 2)
    else:
        return math.ceil((b + a) / 2)


def square_difference(omega: int) -> float:
    """
    Return the sqaure of the sum of th IP and the homo energies for a given omega.
    """
    write(f"\nStarting Calculations for omega={omega}:\n")

    E_N, E_P, E_homo = single_point(omega)
    IP = E_P - E_N
    squared_sum = (IP + E_homo) ** 2

    write(f"(IP + E_homo)^2 = {squared_sum}\n")
    return squared_sum


# run single point calculations and extract energies from output
def single_point(omega: int) -> tuple[float, float, float]:
    """
    Run single point calculations for both neutral and cation species
    and return the tuple (neutral energy, cation energy, homo of neutral species).
    """
    # create new input files with the correct omega
    if os.path.isfile(f"scratch/input/N_{omega}.in") and os.path.isfile(
        f"scratch/input/P_{omega}.in"
    ):
        write("Input files already exist.\n")
    else:
        write("Creating input files...\n")
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

    # run qchem jobs for both N and P
    if os.path.isfile(f"scratch/output/N_{omega}.out") and os.path.isfile(
        f"scratch/output/P_{omega}.out"
    ):
        write("Output files already exist.\n")
    else:
        write("Running calculations...\n")
        subprocess.call(["sbatch", "--wait", "single_point.sh", f"{omega}"])
        write("Calculations done.\n")

    # extract data from the output
    write("\nExtracting data from output...\n")
    E_N, E_homo = extract_energies(f"scratch/output/N_{omega}.out")
    E_P, _ = extract_energies(f"scratch/output/P_{omega}.out")
    write(f"E_N: {E_N}, E_P: {E_P}, E_homo: {E_homo}\n")

    return (E_N, E_P, E_homo)


def extract_energies(file: str) -> tuple[float, float]:
    """
    Extract Total energies and homo energies from single-point calculation output.
    """

    # initializtion and flags
    E_total = 0
    E_homo = 0
    restricted = False
    passed_alpha = False
    passed_beta = False
    symmetry_on_MO = False

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

    if E_total == 0:
        write("ERROR: E_total could not be found")
    if E_homo == 0:
        write("ERROR: E_homo could not be found")

    return (E_total, E_homo)


def main(argv: list[str]) -> int:
    # parse arguments
    parser = ArgumentParser(
        prog="tuning.py",
        description="Find the optimal omega value for long range corrected hybrid functionals.",
    )
    parser.add_argument(
        "-b",
        "--bounds",
        nargs=2,
        action="store",
        default=[200, 800],
        type=int,
        dest="bounds",
        help="Custom bounds for omega, must be two intigers between 0 and 1000.",
    )
    parser.add_argument(
        "-t",
        "--tolerance",
        action="store",
        default=1,
        type=int,
        dest="tolerance",
        help="When the size of range containing omega reaches the tolerance, the center of the range is given as omega. Must be an intiger between 1 and 999.",
    )
    args = parser.parse_args(argv[1:])

    # extract bounds for omega
    if args.bounds[0] < args.bounds[1]:
        lower_bound, upper_bound = args.bounds
    else:
        upper_bound, lower_bound = args.bounds
    if lower_bound < 0 or upper_bound > 1000:
        parser.error("Bounds for omega must be between 0 and 1000.")

    tolerance = args.tolerance
    if tolerance < 1 or tolerance > 999:
        parser.error("Tolerance must be an intiger between 1 and 999.")

    # create directories for files created by this script under the "scratch" directory
    if not os.path.exists("scratch/input"):
        os.makedirs("scratch/input")
    if not os.path.exists("scratch/output"):
        os.makedirs("scratch/output")
    if not os.path.exists("scratch/slurm"):
        os.makedirs("scratch/slurm")

    # create an output file, it will be appended with each calculation.
    with open("tuning.out", "w") as output:
        with open("N.in", "r") as N:
            with open("P.in", "r") as P:
                output.write(
                    f"Output file for omega tuning.\n"
                    f"\n------------------------ N.in ----------------------\n"
                    f"{N.read()}"
                    f"----------------------------------------------------\n"
                    f"\n------------------------ P.in ----------------------\n"
                    f"{P.read()}"
                    f"----------------------------------------------------\n"
                )

    # run the optimization
    optimal_omega = gss(square_difference, lower_bound, upper_bound, tolerance)
    write(
        f"\n----------------------------------------------------\n"
        f"|                   CONVERGED                      |\n"
        f"----------------------------------------------------\n"
        f"\noptimal omega: {optimal_omega}"
    )

    # move scratch directory into a history directory so that it doesn't interfere with the next job
    print("\nArchiving files...")
    if not os.path.exists("history"):
        os.makedirs("history")
    timestamp = str(datetime.datetime.now()).replace(" ", "_").replace(":", "-")
    subprocess.call(["mv", "scratch", f"history/{timestamp}"])

    print("\nFinished.")
    return 0


if __name__ == "__main__":
    main(sys.argv)
