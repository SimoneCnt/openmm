#
# Atomic radius derived from solvent electrostatic charge distribution
# Tested with free energy perturbation with explicit solvent
# Authors:  Mafalda Nina, Dmitrii Belogv, and Benoit Roux
# University of Montreal, June 1996.
# M. Nina and B. Roux. Atomic Radii for Continuum Electrostatics Calculations based on 
# Molecular Dynamics Free Energy Simulations. J. Phys. Chem. B 101: 5239-5248 (1997).
#
# UPDATES:
# --------
#
# GLU and ASP modified December 1st, 1998 by Mafalda Nina
#
# Protonated histidine HSP has been included, January 1999 by Mafalda Nina
# dG_elec = -68.15 kcal/mol (MD/FES) or -68.10 kcal/mol (PBEQ)
#
# Tetramethylamonium (TEA) and ions added by Benoit Roux, January 1999.
#
# Sodium added by Benoit Roux, November 2000
#
# 1. PMF guided optimization of input radii for charged residues (2005)
# J. Chen, W. Im and C. L. Brooks III. Balancing Solvation and Intramolecular 
# Interactions:  Toward a Consistent Generalized Born Force Field. JACS (2006).
#

def get_gbsw_radius(atype, rname):
    """ 
        Return GBWS optimized radius from atom name and residue.
        Returned radius are in angstrom.
    """

    rad = float('nan')

    if atype[0]=='H': rad=0.0
    if atype[0]=='C': rad=2.3
    if atype[0]=='O': rad=1.8
    if atype[0]=='N': rad=2.3
    if atype[0]=='S': rad=2.3

    if atype=='ZN'  : rad=1.30
    if atype=='CAY' : rad=2.67
    if atype=='CAT' : rad=2.06
    if atype=='CY'  : rad=2.04
    if atype=='OY'  : rad=1.52
    if atype=='OT1' : rad=1.40
    if atype=='OT2' : rad=1.40
    if atype=='C'   : rad=2.04
    if atype=='O'   : rad=1.52
    if atype=='N'   : rad=2.03
    if atype=='NT'  : rad=2.03
    if atype=='NL'  : rad=2.03
    if atype=='NR'  : rad=2.03
    if atype=='CA'  : rad=2.86
    if atype=='CB'  : rad=2.67
    if atype=='CG'  : rad=2.46
    if atype=='CG1' : rad=2.46
    if atype=='CG2' : rad=2.46
    if atype=='CD'  : rad=2.44
    if atype=='CD1' : rad=2.44
    if atype=='CD2' : rad=2.44
    if atype=='CE'  : rad=2.10
    if atype=='OG'  : rad=1.64
    if atype=='OG1' : rad=1.64
    if atype=='CGY' : rad=2.77
    if atype=='CDY' : rad=1.98
    if atype=='OEY1': rad=1.40
    if atype=='OEY2': rad=1.40

    if rname=='ARG':
        if atype=='CZ'  : rad=2.20
        if atype=='NH1' : rad=1.70
        if atype=='NH2' : rad=1.70
        if atype=='NE'  : rad=1.70

    if rname=='ASN':
        if atype=='CG'  : rad=1.98
        if atype=='OE1' : rad=1.42
        if atype=='OE2' : rad=1.42
        if atype=='OD1' : rad=1.60
        if atype=='OD2' : rad=1.42
        if atype=='ND2' : rad=2.00

    if rname=='ASP':
        if atype=='CG'  : rad=1.98
        if atype=='OE1' : rad=1.40
        if atype=='OE2' : rad=1.40
        if atype=='OD1' : rad=1.40
        if atype=='OD2' : rad=1.40

    if rname=='CYS':
        if atype=='SG'  : rad=2.00

    if rname=='GLN':
        if atype=='CD'  : rad=1.98
        if atype=='OE1' : rad=1.60
        if atype=='OE2' : rad=1.42
        if atype=='OD1' : rad=1.42
        if atype=='OD2' : rad=1.42
        if atype=='NE2' : rad=2.00

    if rname=='GLU':
        if atype=='CG'  : rad=2.77
        if atype=='CD'  : rad=1.98
        if atype=='OE1' : rad=1.40
        if atype=='OE2' : rad=1.40
        if atype=='OD1' : rad=1.40
        if atype=='OD2' : rad=1.40

    if rname=='GLY':
        if atype=='CA'  : rad=2.38

    if rname=='HIS':
        if atype=='CE1' : rad=1.98
        if atype=='CD2' : rad=1.98
        if atype=='NE2' : rad=1.80
        if atype=='ND1' : rad=1.80

    if rname=='HSD':
        if atype=='CE1' : rad=1.98
        if atype=='CD2' : rad=1.98
        if atype=='NE2' : rad=1.80
        if atype=='ND1' : rad=1.90

    if rname=='HSP':
        if atype=='CE1' : rad=1.98
        if atype=='CD2' : rad=1.98
        if atype=='NE2' : rad=1.90
        if atype=='ND1' : rad=1.90

    if rname=='LYS':
        if atype=='CE'  : rad=2.80
        if atype=='NZ'  : rad=1.80

    if rname=='MET':
        if atype=='SD'  : rad=2.00

    if rname=='PHE':
        if atype=='CE1' : rad=2.00
        if atype=='CE2' : rad=2.00
        if atype=='CD1' : rad=2.00
        if atype=='CD2' : rad=2.00
        if atype=='CZ'  : rad=2.00
        
    if rname=='PRO':
        if atype=='CB'  : rad=1.98
        if atype=='CG'  : rad=1.98
        if atype=='CD'  : rad=1.98

    if rname=='TRP':
        if atype=='CE2' : rad=2.00
        if atype=='CE3' : rad=2.00
        if atype=='CD1' : rad=2.00
        if atype=='CD2' : rad=2.00
        if atype=='CZ2' : rad=2.00
        if atype=='CZ3' : rad=2.00
        if atype=='CH2' : rad=2.00
        if atype=='NE1' : rad=1.85

    if rname=='TYR':
        if atype=='CE1' : rad=2.00
        if atype=='CE2' : rad=2.00
        if atype=='CD1' : rad=2.00
        if atype=='CD2' : rad=2.00
        if atype=='CZ'  : rad=2.00
        if atype=='OH'  : rad=1.85
    
    if rname=='TIP3':
        if atype=='OH2' : rad=2.20

    if rname=='TEA':
        if atype=='N'   : rad=2.15
        if atype=='C1'  : rad=2.30
        if atype=='C2'  : rad=2.30
        if atype=='C3'  : rad=2.30
        if atype=='C4'  : rad=2.30
        if atype=='C5'  : rad=2.30
        if atype=='C6'  : rad=2.30
        if atype=='C7'  : rad=2.30
        if atype=='C8'  : rad=2.30

    if rname=='POT' : rad=2.035
    if rname=='CLA' : rad=2.035
    if rname=='SOD' : rad=1.66

    return rad


