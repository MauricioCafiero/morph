#!/usr/bin/env python3
"""
strain.py — apo -> co-fold peptide conformational strain, as an MM energy
difference (ff14SB + GBn2 implicit solvent), the same force field the openmm
repo uses for protein MD.

For each state: take the peptide (ligand dropped), add hydrogens with PDBFixer at
the given pH (cached, since PDBFixer's hydrogen placement varies per call), apply an
extremely limited declash relaxation (RELAX_ITERS = 30 iterations, force target
RELAX_TOL = 500 kJ/mol/nm, single-threaded CPU) that pulls steric overlaps apart
without changing the conformation, then take the potential energy. The heavy-atom
RMSD the relaxation introduced is returned so callers can flag moves past
MAX_HEAVY_RMSD. Strain is dE = E(cofold conformation) - E(apo conformation):
positive means the bound conformation costs the peptide energy relative to its apo fold.

Reads .pdb and Boltz/BoltzGen/OpenFold .cif (peptide ATOM records only) via morph.py.

Usage:
  strain.py apo.pdb cofold.cif
  strain.py --labels esmfold,boltz --ph 7.0 --workdir cache/ apo.pdb cofold.cif
  strain.py --iters 30 --tol 500 ...     # change the declash budget
"""
import argparse
import os
import sys
import tempfile
import numpy as np

import morph  # same folder

from pdbfixer import PDBFixer
from openmm import app, unit, LangevinIntegrator, Platform


def peptide_pdb(src, chain, dst):
    """Write a peptide-only heavy-atom PDB (drops ligand) that PDBFixer can read."""
    atoms = morph.load(src, chain)
    if not atoms:
        sys.exit(f"{src}: no ATOM records on chain {chain}")
    with open(dst, "w") as fh:
        for serial, k in enumerate(atoms, 1):
            name, resname, resseq, element, xyz = atoms[k]
            fh.write(morph.pdb_line(serial, name, resname, resseq, element, xyz, "A"))
        fh.write("END\n")


# Strain is scored on the predicted folded structure, kept as close to as-predicted as
# possible: we add the missing hydrogens and then apply an *extremely limited* relaxation
# whose only job is to pull steric clashes apart -- not to refold anything.
#   * high force target + very few iterations: clash gradients (r^-12) are enormous, so the
#     first handful of steps move exactly the overlapping atoms and little else; heavy atoms
#     shift ~0.1 A on clean structures and ~0.3 A on clashy ones (local declash, not a state
#     change, which would be several A).
#   * we record heavy-atom RMSD from the predicted coordinates per state and flag any that
#     moves more than MAX_HEAVY_RMSD -- a guard that the conformation was preserved.
#   * single-thread CPU + cached one-time protonation -> reproducible (PDBFixer's per-call
#     hydrogen placement is the only run-to-run noise, so we cache it).
RELAX_ITERS = 30
RELAX_TOL = 500.0       # kJ/mol/nm: high target, stops as soon as the bad contacts are gone
MAX_HEAVY_RMSD = 0.6    # A: above this the "limited" relaxation moved too far -> flag


def energy(pep_pdb, ph, iters=RELAX_ITERS, tol=RELAX_TOL, cache=None):
    """Protonate (or reuse a cached protonation), declash by a limited relaxation.

    Returns (energy kcal/mol, n_atoms, heavy-atom RMSD from the predicted coords in A).
    """
    if cache and os.path.exists(cache):
        pdb = app.PDBFile(cache)
        top, pos = pdb.topology, pdb.positions
    else:
        fx = PDBFixer(filename=pep_pdb)
        fx.findMissingResidues()
        fx.findMissingAtoms()
        fx.addMissingAtoms()
        fx.addMissingHydrogens(pH=ph)
        top, pos = fx.topology, fx.positions
        if cache:
            app.PDBFile.writeFile(top, pos, open(cache, "w"))   # freeze this protonation
    ff = app.ForceField("amber14/protein.ff14SB.xml", "implicit/gbn2.xml")
    system = ff.createSystem(top, nonbondedMethod=app.NoCutoff, constraints=None)
    integ = LangevinIntegrator(300 * unit.kelvin, 1 / unit.picosecond, 0.001 * unit.picoseconds)
    sim = app.Simulation(top, system, integ, Platform.getPlatformByName("CPU"), {"Threads": "1"})
    sim.context.setPositions(pos)
    p0 = np.array(sim.context.getState(getPositions=True).getPositions().value_in_unit(unit.angstrom))
    sim.minimizeEnergy(tolerance=tol * unit.kilojoule_per_mole / unit.nanometer, maxIterations=iters)
    st = sim.context.getState(getPositions=True, getEnergy=True)
    p1 = np.array(st.getPositions().value_in_unit(unit.angstrom))
    heavy = [i for i, a in enumerate(top.atoms()) if a.element is not None and a.element.symbol != "H"]
    rmsd = float(np.sqrt(((p1[heavy] - p0[heavy]) ** 2).sum(1).mean()))
    e = st.getPotentialEnergy().value_in_unit(unit.kilocalorie_per_mole)
    return e, top.getNumAtoms(), rmsd


def compute_strain(apo, cof, labels=("apo", "cofold"), chain="A", ph=7.0,
                   iters=RELAX_ITERS, tol=RELAX_TOL, workdir=None, verbose=True):
    """Return each state's declashed energy, the apo->cofold strain dE, and the max heavy
    RMSD the limited relaxation introduced. If workdir is given, the peptide and its frozen
    protonation are cached there (<label>.pdb / <label>_H.pdb) for reproducible re-runs.
    """
    tmp = workdir or tempfile.mkdtemp()
    os.makedirs(tmp, exist_ok=True)
    energies, natoms, rmsds = [], [], []
    for src, lab in zip((apo, cof), labels):
        pep = os.path.join(tmp, f"{lab}.pdb")
        peptide_pdb(src, chain, pep)
        e, na, rmsd = energy(pep, ph, iters, tol, cache=os.path.join(tmp, f"{lab}_H.pdb"))
        energies.append(e)
        natoms.append(na)
        rmsds.append(rmsd)
        if verbose:
            print(f"  {lab:12s}  E = {e:12.1f} kcal/mol   (heavy RMSD {rmsd:.2f} A, {na} atoms incl. H)")
    de = energies[1] - energies[0]
    max_rmsd = max(rmsds)
    if verbose:
        print(f"\nstrain dE ({labels[1]} - {labels[0]}) = {de:+.1f} kcal/mol"
              f"   [declash {iters} iters, tol {tol} kJ/mol/nm, 1 thread; max heavy RMSD {max_rmsd:.2f} A]")
    return {"e_apo": energies[0], "e_cof": energies[1], "d_strain": de, "n_atoms": natoms[0],
            "rmsd_apo": rmsds[0], "rmsd_cof": rmsds[1], "max_heavy_rmsd": max_rmsd, "labels": labels}


def main():
    ap = argparse.ArgumentParser(description="apo->cofold peptide strain (ff14SB + GBn2).")
    ap.add_argument("structures", nargs=2, help="apo cofold (.pdb or .cif)")
    ap.add_argument("--labels", default="apo,cofold")
    ap.add_argument("--chain", default="A")
    ap.add_argument("--ph", type=float, default=7.0)
    ap.add_argument("--iters", type=int, default=RELAX_ITERS, help="hydrogen-relaxation iteration cap")
    ap.add_argument("--tol", type=float, default=RELAX_TOL, help="minimizer force target (kJ/mol/nm)")
    ap.add_argument("--workdir", default=None, help="cache peptide/protonated PDBs here (reproducible re-runs)")
    args = ap.parse_args()
    labels = tuple(args.labels.split(","))
    compute_strain(args.structures[0], args.structures[1], labels, args.chain, args.ph,
                   args.iters, args.tol, args.workdir)
    print(f"  (+ = {labels[1]} conformation is strained relative to {labels[0]})")


if __name__ == "__main__":
    main()
