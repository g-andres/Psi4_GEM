import numpy as np
import psi4
import sys
import math
from psi4.driver import constants

molA = psi4.geometry("""
 O  0.0 0.0 -0.125477814671
 H  1.45357733502   -0.000188972612457  0.996263612871
 H  -1.45376630763  0.000188972612457   0.995885667646
 symmetry c1
 no_reorient
 units bohr
 no_com
""")

molB = psi4.geometry("""
 O -0.108659252162 0.0 6.43678512549
 H 1.58963761598 -0.000188972612457 5.73890926769
 H 0.135682335744 0.0 8.25659138345
 symmetry c1
 no_reorient
 units bohr
 no_com
""")

methodname = 'WB97X-D'
basisname = 'aug-cc-pVTZ'
gembasisname = 'gemherm'
K_exchange = 7.075662


# The singular value threshold when inverting the metric
epsilon = 1e-12
# The denominator shift to be applied when inverting the metric
regularizer = 0.000001
# Whether the charge should be constrained
constrain = True

# Force Cartesian, so I can compare with the reference outputs
psi4.set_options({'puream' : False,
                  'print' : 1,
                  'scf_type' : 'df'})

#### Get the unperturbed wavefunctions and energies
eA, wfnA = psi4.energy(methodname+'/'+basisname, molecule=molA, return_wfn=True)
eB, wfnB = psi4.energy(methodname+'/'+basisname, molecule=molB, return_wfn=True)
psi4.core.print_out("------ FINISHED ENERGY OF MONOMER ------\n\n")

def invert(mat):
    if epsilon == 0:
        # Simple inversion via Cholesky factorization - no filtering
        c = np.linalg.inv(np.linalg.cholesky(mat))
        return np.dot(c.T,c)
    else:
        # Eigendecompose, filter eigenvalues, invert, then reconstruct
        evals, evecs = np.linalg.eigh(mat)
        evals = np.where(evals < epsilon, 0.0, 1.0/evals)
        return np.einsum('ik,k,jk->ij', evecs, evals, evecs)

def build_gem_field(wfn, add_coulomb, add_exchange):
    """ Builds external potential objects for Coulomb and exchange (via overlap) terms """
    psi4.core.print_out("Computing GEM representation")
    orbital_basis = wfn.basisset()
    molecule = orbital_basis.molecule()
    aux_basis = psi4.core.BasisSet.build(molecule, "ORBITAL", gembasisname)
    aux_basis.apply_hermite_normalization()
    zero_basis = psi4.core.BasisSet.zero_ao_basis_set()
    mints = psi4.core.MintsHelper(orbital_basis)

    field = psi4.QMMM().extern

    # J electronic terms
    J_PQ = np.squeeze(mints.ao_eri(aux_basis, zero_basis, aux_basis, zero_basis))
    J_PQinv = invert(J_PQ + regularizer*np.eye(aux_basis.nbf()))
    J_Pmn = np.squeeze(mints.ao_eri(aux_basis, zero_basis, orbital_basis, orbital_basis))

    d =  2 * np.einsum('Qmn,mn->Q', J_Pmn, wfn.Da())
    fit_coefficients = np.einsum('PQ,Q->P', J_PQinv, d)

    if constrain:
        # Make sure the aux basis holds the correct charge using the method described in
        # https://doi.org/10.1080/00268976.2010.518982 particularly equations 16 and 17.
        def dfact(n):
            """ A simple way to compute double factorials using built-in functions. """
            if n <= 1:
                return 1.0
            if n%2 == 1:
                # odd case
                k = (n+1)//2
                return math.factorial(2*k) / (2**k * math.factorial(k))
            else:
                # odd case
                k = n//2
                return 2**k * math.factorial(k)

        # Find q; the integral of each basis function 
        nshell = aux_basis.nshell()
        q = []
        for sh in range(nshell):
            shell = aux_basis.shell(sh)
            nprim = shell.nprimitive
            if nprim > 1:
                raise("This code currently assumes uncontracted GEM basis sets")
            L = shell.am
            ex = shell.exp(0)
            coef = shell.coef(0)
            prefac = coef * np.sqrt((np.pi**3)/(2**L * ex**(L+3)))
            for lx in range(L,-1,-1):
                for lz in range(L-lx+1):
                    ly = L - lx - lz
                    # Only even parity functions have a nonzero integral
                    if (lx%2==0) and (ly%2==0) and (lz%2==0):
                        q.append(prefac*dfact(lx-1)*dfact(ly-1)*dfact(lz-1))
                    else:
                        q.append(0.0)
        # Fine the expected number of electrons
        natom = molecule.natom()
        Q = 0
        for i in range(natom):
            Q += molecule.Z(i)
        Q -= molecule.molecular_charge()

        old_nelec = np.dot(fit_coefficients, q)
        psi4.core.print_out(f"Number of electrons before constraining: {old_nelec:.6f}")

        # Figure out the Lagrange multiplier and constrain the fit coefficients
        qdot = np.dot(J_PQinv, q)
        lam = (Q - np.dot(qdot,  d)) / np.dot(qdot, q)
        fit_coefficients += lam * qdot
        new_nelec = np.dot(fit_coefficients, q)
        psi4.core.print_out(f"Number of electrons after constraining: {new_nelec:.6f}")

    coefs_vec = psi4.core.Vector(aux_basis.nbf())
    for n, coef in enumerate(fit_coefficients): coefs_vec.set(n, coef)
    field.addBasis(aux_basis, coefs_vec)

    # K exchange terms
    if add_exchange:
        K_coefs_vec = psi4.core.Vector(aux_basis.nbf())
        for n, coef in enumerate(fit_coefficients): K_coefs_vec.set(n, K_exchange*coef)
        field.addExchangeBasis(aux_basis, K_coefs_vec)

    # Nuclear terms
    for n in range(molecule.natom()):
        field.addCharge(molecule.Z(n), molecule.x(n), molecule.y(n), molecule.z(n))
    return field

# Coulomb field
field = build_gem_field(wfnA, add_coulomb=True, add_exchange=False)
psi4.core.set_global_option_python('EXTERN', field)
eB_Jperturbed = psi4.energy(methodname+'/'+basisname, molecule=molB)

# Coulomb+exchange field
field = build_gem_field(wfnA, add_coulomb=True, add_exchange=True)
psi4.core.set_global_option_python('EXTERN', field)
eB_JKperturbed = psi4.energy(methodname+'/'+basisname, molecule=molB)

psi4.core.print_out("\nCoulomb energy: {:6} kcal/mol\n".format(constants.hartree2kcalmol*(eB_Jperturbed-eB)))
psi4.core.print_out("\nExchange energy: {:6} kcal/mol\n".format(constants.hartree2kcalmol*(eB_JKperturbed-eB_Jperturbed)))
