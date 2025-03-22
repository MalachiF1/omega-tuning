# Omega Tuning

This script tunes the range seperating parameter, $\omega$, (in units of 10$^{-3}$ bohr$^{-1}$) for a given system and RSH functional. The tuning is done by enforcing DFT's version of koopman's theorem, finding the value of $\omega$ such that: $\text{IP} = -\epsilon_{\text{homo}}$.\\
In practice, this is done via a golden-section search algorithm to minimize the value of $(\text{IP} + \epsilon_{\text{homo}})^2$ with respect to $\omega$. Note that this requires the premise that $(\text{IP} + \epsilon_{\text{homo}})^2$ is unimodal with respect to $\omega$ for the range of $\omega$ under investigation

## How to Use:

This script uses Qchem, and was tested on Qchem version 5.4. \\
Create an input file for your system and RSH functional, named `N.in`, and add `omega xxx` to the `$rem` section.
Create another input file for the cation of your system, with the same functional and basis (if your original system is a +1 cation this input file will be for a +2 cation and so forth) and name it `P.in`. This file should have `omega xxx` in the `$rem` section as well. The job types should be `job sp` (single-point). With the input files ready, run the python script.

Example for `N.in` for H$\_2$0 and the LRC-$\omega$PBEh functional:

```
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
```

Example for `P.in` for H$_2$O and the LRC-$\omega$PBEh functional:

```
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
```

A sample run for these two input files is included in the `sample_run` directory.
