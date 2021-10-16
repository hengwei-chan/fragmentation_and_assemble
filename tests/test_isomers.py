import pytest

from rdkit import Chem

import datamol as dm


def test_enumerate_tautomers():
    mol = dm.to_mol("OC1=CC2CCCCC2[N:1]=C1")

    mols = dm.enumerate_tautomers(mol, n_variants=10)

    assert {dm.to_smiles(m) for m in mols} == {"O=C1C=[N:1]C2CCCCC2C1", "OC1=CC2CCCCC2[N:1]=C1"}


def test_enumerate_stereo():
    mol = dm.to_mol("OC1=CC2CCCCC2[N:1]=C1")

    mols = dm.enumerate_stereoisomers(mol, n_variants=10)

    assert {dm.to_smiles(m) for m in mols} == {
        "OC1=C[C@@H]2CCCC[C@@H]2[N:1]=C1",
        "OC1=C[C@@H]2CCCC[C@H]2[N:1]=C1",
        "OC1=C[C@H]2CCCC[C@@H]2[N:1]=C1",
        "OC1=C[C@H]2CCCC[C@H]2[N:1]=C1",
    }


def test_enumerate_structural():
    mol = dm.to_mol("CCCCC")  # pentane has only three structural isomers
    mols_iso = dm.enumerate_structisomers(
        mol,
        n_variants=5,
        allow_cycle=False,
        depth=2,
        allow_double_bond=False,
        allow_triple_bond=False,
    )
    mols_cyclo_iso = dm.enumerate_structisomers(mol, n_variants=5, depth=2, allow_cycle=True)

    assert {dm.to_smiles(m) for m in mols_iso} == {"CCC(C)C", "CC(C)(C)C"}
    # expect 3 molecules with cycles
    assert sum([Chem.rdMolDescriptors.CalcNumRings(x) == 1 for x in mols_cyclo_iso]) == 3

    # mols_cyclo_iso_double = dm.enumerate_structisomers(
    #     mol, n_variants=10, allow_cycle=True, allow_double_bond=True
    # )
    # should have mol with double link
    # assert sum(["=" in dm.to_smiles(x) for x in mols_cyclo_iso_double]) > 0
