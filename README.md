# HAVMOL
Early PC's program for ab initio Hartree -  Fock calculations
Title: HAVMOL.

Description: HAVMOL is a simple computer program designed to perform calculations of molecular struc-ture and properties by a quantum mechanical approach, according to the ab initio SCF MO the-ory of Roothaan1 and Hall2.

The pedigree of this program is the following: HAVMOL is a version of MICROMOL, an adaption made by Dr. Susan Colwell of the Cambridge Analytic Derivatives Package (CADPAC). That is a suite of programs developed in Cambridge, United Kingdom, between 1980 and 1983 by Dr. Roger D. Amos. Those, in turn, originated from former versions of Dupuis and King's HONDO program.

HAVMOL includes an intensive use of memory and resources of personal computers to in-crease portability and calculations in very modest CPU’s, input simplicity and easiness of use by both initiated and non experienced scientists, anywhere. It allows to perform top level calcu-lations of molecular structure and other properties in any place where personal computers exist, including home and experimental laboratories without mainframe or work station terminals. Substantial  enhancements have been made to code and integral management in certain critical routines and improvements have been introduced to the input and output of information with respect to MICROMOL. A very easing feature is added to this version which allow a simple in-put from either Cartesian or internal coordinate files obtained by TC-HABANA program out-puts, like semiempirical SCF-MO HAVPAC, molecular mechanics. Internal coordinates in the widespread MOPAC-AMPAC, or Gaussian's Z matrix formats are also allowed as input. A con-tracted two electron integral external file is built to save disk space, respect to former versions.

QCPE Origin: This program was among the Quantum Chemistry Program Exchange repository at the University of Indiana, Bloomington, IN underthe code QCPE 182, 1997. As QCPE is no longer active ist is being relocated in GitHub. Most of the compiling facilities could be no longer useful.

Features:

This program is able to calculate at the Roothaan's SCF-MO level both close and open shell systems, analytic first derivatives of the energy, and hence it can perform geometry optimiza-tions and calculate numerical force constants. Thermodynamic quantities on the grounds of the statistical approach for ideal gases at standard and other desired temperatures, rotational con-stants, dipole and quadrupole moments and isotopic effects are also possible to be calculated for either input and optimized molecular structures. It can take into account isotopic substitu-tion effects.

Briefly, program abilities are the following:

	Evaluation of 1 and 2 electron integrals over contracted cartesian gaussian basis func-tions of type s, p and d.
	SCF iterative procedure for closed and high spin open shell wave functions. Ghost terms can be used to perform counterpoise basis set superposition effect corrections and float-ing function calculations.
	Calculation of the energy gradients, and derivatives of both one and two electron inte-grals.
	Use of the gradients for automatic geometry optimization by Murtagh-Sargent minimi-zation routine.
	Use of the gradients for the calculation of force constants by numerical differentiation.
	Force constants matrix is used to calculate molecular vibrational energies and partition function  derived quantities. It can include isotopic effects on vibrational terms.
	Calculation of rotational constants and partition function derived quantities of the opti-mized and other given geometries of the  system. It can include isotopic effects.
	Calculation of some electrical properties of the desired system, including dipole and quadrupole moments, Mulliken gross charges and electric field at nuclei.
	Calculation of translational, rotational, and vibrational thermodynamic constants at 298.16 K and other desired temperatures as ideal gases.
	When geometry optimization is performed, an atomic cartesian coordinate output file is produced, to be compatible with other TC HABANA personal computer program inputs.

Installation Instructions:

The programs can be used as Windows CMD executable and no further installation is required.
The programs were developed for IBM PC/XT/AT and PS/n, or compatibles, and they run under patched DOS v. 3.2 or higher versions. Eventually, they has been compiled and linked with overlays, using a MICROSOFT optimizing FORTRAN compiler version 5.1. In the case of HAVMOLE, the overlay structure forces to set the program in the default directory or in DOS path and not to change the executable file name. See in your operating system manuals the DOS environment related instructions, such as PATH and SET for further  information. Typical HAVMOL.PIF and HAVMOLE.PIF command programs are provided in the package for run-ning the corresponding DOS programs in MS-WINDOWS environment, taking advantage of multitasking abilities of this i80386 protected mode operating system.

The DOS configuration must allow at least 10 files to be opened simultaneously (see in your DOS manual reference to the CONFIG.SYS file). An approx. 512 kb virtual disk (see in your DOS manual reference to the VDISK.SYS or RAMDRIVE.SYS device drivers) is recommended to be used for allocating an scratch file (see Section 10) and saving much computing time in long calculations. Larger virtual disk capacities could be required in certain cases.

Usage Instructions:

HAVMOL.EXE which has been optimized to be run in any  personal  computer  (min-imum 640 KB RAM) with or without the i 8087/287/387 arithmetic coprocessors. It is dimensioned for:
	(a) a maximum of 63 basis functions.
	(b) a maximum of 30 shells.
	(c) a maximum of 12 atoms.
	(d) a maximum of 10 primitive  gaussians  in  any  contracted function.
	(e) a maximum of 110 unique primitives in total.

This is the fastest and most portable version.


  	HAVMOLE.EXE which has been optimized to be run in any personal computer (with a minimum of 640 KB RAM) with an installed i 8087/287/387 arithmetic coprocessor, 486 DX processors or compatible higher in the series. Dimensioned for:
	(a) a maximum of 100 basis functions.
	(b) a maximum of 60 shells.
	(c) a maximum of 24 atoms.
	(d) a maximum of 10 primitive gaussians in any contracted function.
	(e) a maximum of 200 unique primitives in total.

This version is dimensioned for full use of the available DOS memory which can calculate larger molecular systems with larger basis sets. It uses an optimized overlay structure to save code room in memory.

  	HAVMOLW.EXE which has been optimized to be run in any 386 or higher personal computer under MS WINDOWS environment. Dimensioned for:
	(a) a maximum of 127 basis functions.
	(b) a maximum of 100 shells.
	(c) a maximum of 30 atoms.
	(d) a maximum of 10 primitive gaussians in any contracted function.
	(e) a maximum of 508 unique primitives in total.

This version is dimensioned for full use of the available RAM, and eventually, virtual memory from hard disk. It can calculate the largest molecular systems with largest basis sets in this se-ries of programs.

Dependencies:

The DOS configuration must allow at least 10 files to be opened simultaneously (see in your DOS manual reference to the CONFIG.SYS file). An approx. 512 kb virtual disk (see in your DOS manual reference to the VDISK.SYS or RAMDRIVE.SYS device drivers) is recommended to be used for allocating an scratch file (see Section 10) and saving much computing time in long calculations. Larger virtual disk capacities could be required in certain cases.

Contact Information: lmc@fq.uh.cu

Citation Information:

Montero, L.A., HAVMOL: Calculations of Molecular Structure and Properties by a Quantum Mechanical Approach, According to the ab initio SCF MO Theory of Roothaan. Quantun Chemistry Program Exchange program: QCMP 182, University of Indiana, Bloomington, 1997.
