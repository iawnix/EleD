"""Regression cases from the descriptor audit, with nonorthogonal densities."""
from dataclasses import replace
from copy import deepcopy

import numpy as np
import pytest
from rdkit import Chem

from eled import (XTBResult, CalcPiDescriptors, CalcDescriptorReport,
                  CalcMayerBondDescriptors, CalcFrontierBondDescriptors,
                  DescriptorSensitivity)
from eled.basis import NativePDirections
from eled.input import _input_metadata, _xyz_from_mol
from eled.options import CalculationSettings


def pi_system():
    mol = Chem.AddHs(Chem.MolFromSmiles('c1ccccc1'))
    conf = Chem.Conformer(12); conf.Set3D(True)
    for i in range(6):
        for j, radius in [(i, 1.4), (i+6, 2.4)]:
            conf.SetAtomPosition(j, (radius*np.cos(i*np.pi/3), radius*np.sin(i*np.pi/3), 0.))
    mol.AddConformer(conf)
    mapping = tuple(np.repeat(np.arange(6), 3)) + tuple(range(6, 12))
    angular = (1,)*6 + (0,)*6
    n = len(mapping)
    s = np.eye(n); s[2,5] = s[5,2] = .2
    eig, u = np.linalg.eigh(s)
    c = (u/np.sqrt(eig)) @ u.T
    meta = _input_metadata(mol, -1, 0, 0, _xyz_from_mol(mol, -1))
    meta.update(CalculationSettings=CalculationSettings(), engine={'sha256':'fixture'}, calculation_status='success')
    r = XTBResult('GFN2-xTB', None, -1., 1., n, 30, 0,
        orbital_energies_ev=tuple(range(n)), orbital_occupations=tuple([2.,2.,1.]*6+[0.]*6),
        partial_charges=(0.,)*12, wiberg_bond_orders=tuple(map(tuple, np.zeros((12,12)))),
        orbital_coefficients=tuple(map(tuple,c)), overlap_matrix=tuple(map(tuple,s)),
        ao_atom_indices=mapping, ao_shell_indices=mapping, shell_atom_indices=tuple(range(12)),
        shell_angular_momenta=angular, ao_p_directions=NativePDirections(mapping, angular), metadata=meta)
    return mol, r


@pytest.mark.parametrize('scheme', ['mulliken','lowdin'])
def test_pi_metric_charge_reference_and_missing_rows(scheme):
    mol, r = pi_system()
    report = CalcPiDescriptors(mol, r, scheme=scheme)
    assert report['atoms'][0]['pi_n_e'] == pytest.approx(1.)
    assert report['atoms'][0]['pi_q_e'] == pytest.approx(0.)
    assert len(report['atoms']) == 12 and len(report['bonds']) == mol.GetNumBonds()
    assert report['atoms'][6]['pi_n_e'] is None
    assert report['atoms'][6]['missing']['pi_n_e'] == 'atom_is_not_sp2_or_aromatic'
    assert report['bonds'][6]['pi_bo'] is None
    explicit = CalcPiDescriptors(mol, r, pi_ref={0: 2.})
    assert explicit['atoms'][0]['pi_q_e'] == pytest.approx(1.)
    radical = CalcPiDescriptors(mol, replace(r, number_of_unpaired_electrons=1))
    assert radical['bonds'][0]['pi_bo'] is None
    assert radical['bonds'][0]['missing']['pi_bo'] == 'spin_density_unavailable_for_open_shell'
    assert radical['atoms'][0]['pi_n_e'] is not None


@pytest.mark.parametrize('scheme', ['mulliken','lowdin'])
def test_pi_invariant_to_ao_permutation_sign_phase_and_rigid_motion(scheme):
    mol, r = pi_system()
    baseline = CalcPiDescriptors(mol, r, scheme=scheme)
    rng = np.random.default_rng(71)
    perm = rng.permutation(len(r.ao_atom_indices)); signs = rng.choice([-1,1],len(perm))
    c = np.asarray(r.orbital_coefficients); s = np.asarray(r.overlap_matrix)
    changed = replace(r, orbital_coefficients=tuple(map(tuple,c[perm]*signs[:,None]*signs[None,:])),
        overlap_matrix=tuple(map(tuple,s[np.ix_(perm,perm)]*signs[:,None]*signs[None,:])),
        ao_atom_indices=tuple(np.array(r.ao_atom_indices)[perm]),
        ao_shell_indices=tuple(np.array(r.ao_shell_indices)[perm]),
        ao_p_directions=tuple(map(tuple,np.asarray(r.ao_p_directions)[perm]*signs[:,None])))
    permuted = CalcPiDescriptors(mol, changed, scheme=scheme)
    rotation, _ = np.linalg.qr(rng.normal(size=(3,3)))
    transform = np.eye(len(perm))
    for i in range(6): transform[i*3:i*3+3,i*3:i*3+3] = rotation
    rotated = Chem.Mol(mol)
    for i, pos in enumerate(mol.GetConformer().GetPositions()):
        rotated.GetConformer().SetAtomPosition(i, rotation @ pos + [3.,-2.,7.])
    moved = replace(r, orbital_coefficients=tuple(map(tuple,transform@c)),
        overlap_matrix=tuple(map(tuple,transform@s@transform.T)),
        metadata={**r.metadata, **_input_metadata(rotated,-1,0,0,_xyz_from_mol(rotated,-1))})
    motion = CalcPiDescriptors(rotated, moved, scheme=scheme)
    for other in [permuted,motion]:
        for entity, fields in [('atoms',['pi_n_e','pi_q_e']),('bonds',['pi_bo'])]:
            for a,b in zip(baseline[entity],other[entity]):
                for field in fields:
                    if a[field] is None: assert b[field] is None
                    else: assert b[field] == pytest.approx(a[field],abs=1e-10)


def test_pi_bond_is_squared_density_projection_not_overlap_magnitude():
    mol,r = pi_system()
    c = np.eye(24); c[np.ix_([2,5],[2,5])] = np.array([[1,1],[1,-1]])/np.sqrt(2)
    occ=list(r.orbital_occupations);occ[2]=2.;occ[5]=0.
    r=replace(r,orbital_coefficients=tuple(map(tuple,c)),overlap_matrix=tuple(map(tuple,np.eye(24))),orbital_occupations=tuple(occ))
    for scheme in ['lowdin','mulliken']:
        report=CalcPiDescriptors(mol,r,scheme=scheme)
        assert report['bonds'][0]['pi_bo'] == pytest.approx(1.)  # S_AB=0, D_AB=1


@pytest.mark.parametrize('geometry,reason', [('collinear','local_plane_degenerate'),('warped','local_plane_nonplanar')])
def test_pi_rejects_undefined_plane(geometry,reason):
    mol,r=pi_system()
    if geometry=='collinear':
        for i in range(mol.GetNumAtoms()):mol.GetConformer().SetAtomPosition(i,(i,0.,0.))
    else: mol.GetConformer().SetAtomPosition(0,(1.4,0.,2.))
    r=replace(r,metadata={**r.metadata,**_input_metadata(mol,-1,0,0,_xyz_from_mol(mol,-1))})
    report=CalcPiDescriptors(mol,r)
    assert all(row['pi_n_e'] is None and row['missing']['pi_n_e']==reason for row in report['atoms'][:6])


def test_missing_component_metadata_never_guesses_cartesian_order():
    mol,r=pi_system()
    report=CalcPiDescriptors(mol,replace(r,ao_p_directions=()))
    assert report['atoms'][0]['pi_n_e'] is None
    assert report['atoms'][0]['missing']['pi_n_e']=='p_component_directions_unavailable'


@pytest.mark.parametrize('func',[CalcMayerBondDescriptors,CalcFrontierBondDescriptors,CalcPiDescriptors])
def test_mapping_rejects_changed_order_geometry_and_absent_provenance(func):
    mol,r=pi_system()
    renumbered=Chem.RenumberAtoms(mol,[6,1,2,3,4,5,0,7,8,9,10,11])
    with pytest.raises(ValueError,match='structure or atom order'):func(renumbered,r)
    with pytest.raises(ValueError,match='provenance'):func(mol,replace(r,metadata={}))
    mol.GetConformer().SetAtomPosition(0,(9,9,9))
    with pytest.raises(ValueError,match='coordinates'):func(mol,r)


def test_report_forwards_all_frontier_group_parameters_and_propagates_invalid_maps():
    mol,r=pi_system()
    # Two filled orbitals separated by .05 eV must group with a .1 eV threshold.
    energies=tuple([-20.+i for i in range(16)]+[-1.,-.95,1.,2.,3.,4.,5.,6.])
    r=replace(r,orbital_energies_ev=energies,orbital_occupations=(2.,)*18+(0.,)*6,number_of_electrons=36)
    report=CalcDescriptorReport(mol,r,degeneracy_ev=.1)
    assert report['frontier']['homo']['orb_idx']==[16,17]
    assert report['frontier_bonds']['bonds'][0]['orb_idx']==[16,17]
    bad=replace(r,ao_shell_indices=tuple([999]*24))
    with pytest.raises(ValueError,match='shell map'): CalcDescriptorReport(mol,bad)


@pytest.mark.parametrize('change',['geometry','charge','engine','method','settings','analysis','temperature'])
def test_sensitivity_rejects_unmatched_comparisons(change):
    mol,r=pi_system(); a=CalcDescriptorReport(mol,r); b=deepcopy(a)
    if change=='geometry':b['provenance']['geometry_sha256']='other'
    elif change=='charge':b['provenance']['charge']=1
    elif change=='engine':b['provenance']['engine']['sha256']='other'
    elif change=='method':b['values']['xtb_method']='GFN1-xTB'
    elif change=='settings':b['provenance']['CalculationSettings']['accuracy']=.1
    elif change=='analysis':b['analysis_settings']['degeneracy_ev']=.5
    # Unchanged 300 K report supplied as 1000 K must also fail.
    with pytest.raises(ValueError):DescriptorSensitivity(a,temperature_1000=b)


def test_sensitivity_includes_atom_and_bond_fields():
    mol,r=pi_system(); a=CalcDescriptorReport(mol,r)
    b=CalcDescriptorReport(mol,replace(r,method='GFN1-xTB'))
    out=DescriptorSensitivity(a,gfn1=b,gfn2=a)
    assert out['atoms']['0']['shell_populations.mul_p_e']['gfn12_delta']==pytest.approx(0.)
    assert out['bonds']['0']['frontier_bonds.homo.mul_ovlp_e']['gfn12_delta']==pytest.approx(0.)


def test_degenerate_rotations_preserve_density_and_frontier_bonds():
    mol, r = pi_system()
    occ = (2.,)*18 + (0.,)*6
    r = replace(r, orbital_occupations=occ, number_of_electrons=36,
                orbital_energies_ev=(-1.,)*18+(1.,)*6)
    rotation = np.eye(24)
    rng = np.random.default_rng(72)
    for block in (slice(0,18), slice(18,24)):
        size = len(range(*block.indices(24)))
        rotation[block,block], _ = np.linalg.qr(rng.normal(size=(size,size)))
    changed = replace(r, orbital_coefficients=tuple(map(tuple,np.asarray(r.orbital_coefficients)@rotation)))
    for func, entity, fields in [(CalcPiDescriptors,'atoms',['pi_n_e','pi_q_e']),
                                (CalcPiDescriptors,'bonds',['pi_bo']),
                                (CalcFrontierBondDescriptors,'bonds',['mul_ovlp_e']),
                                (CalcMayerBondDescriptors,'bonds',['mbo'])]:
        for a,b in zip(func(mol,r)[entity],func(mol,changed)[entity]):
            for key in fields:
                if a[key] is None: assert b[key] is None
                else: assert a[key] == pytest.approx(b[key],abs=1e-12)


@pytest.mark.integration
@pytest.mark.parametrize('method',['gfn0','gfn1','gfn2'])
@pytest.mark.parametrize('smiles',['c1ccccc1','c1ccncc1','C=C','c1ccoc1'])
def test_real_density_interface_parity_and_pi_rotation(method,smiles):
    import os
    if os.environ.get('ELED_EXPECT_BASIS') != '1': pytest.skip('set ELED_EXPECT_BASIS=1 for native density tests')
    from rdkit.Chem import AllChem
    from helpers import calculate_xtb
    from eled import CalcDensityDescriptors
    from eled.runtime import _find_xtb_library
    mod = pytest.importorskip('eled.rdXTBTools')
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    assert AllChem.EmbedMolecule(mol,randomSeed=73)==0
    AllChem.MMFFOptimizeMolecule(mol)
    lib = _find_xtb_library()
    result = calculate_xtb(mol,method=method,keep_overlap=True,xtb_library=lib)
    ok, compiled = mod.RunMol(mol,method=method,keepOverlap=True,xtbLibrary=lib)
    assert ok
    a,b = CalcDensityDescriptors(mol,result),CalcDensityDescriptors(mol,compiled)
    def assert_same(left,right):
        if isinstance(left,dict):
            assert left.keys()==right.keys()
            for key in left: assert_same(left[key],right[key])
        elif isinstance(left,(list,tuple)):
            assert len(left)==len(right)
            for x,y in zip(left,right): assert_same(x,y)
        elif isinstance(left,float): assert left==pytest.approx(right,abs=1e-8)
        else: assert left==right
    assert_same(a,b)
    assert abs(a['density']['n_err_e'])<1e-8
    for scheme in ['mulliken','lowdin']:
        assert a['shell_populations']['shell_err_e'][scheme]<1e-8
    for i,row in enumerate(a['shell_populations']['atoms']):
        if method!='gfn0': assert row['mul_qrec_e']==pytest.approx(result.partial_charges[i],abs=1e-7)
    # libxtb's closed-shell Wiberg implementation uses the same PS pair contraction.
    for row in a['mayer_bonds']['bonds']: assert row['mbo']==pytest.approx(row['wbo'],abs=1e-7)
    rotation,_=np.linalg.qr(np.random.default_rng(74).normal(size=(3,3)))
    moved=Chem.Mol(mol)
    for i,pos in enumerate(mol.GetConformer().GetPositions()):
        moved.GetConformer().SetAtomPosition(i,rotation@pos+[2.,-3.,7.])
    rotated=calculate_xtb(moved,method=method,keep_overlap=True,xtb_library=lib)
    for scheme in ['mulliken','lowdin']:
        baseline=CalcPiDescriptors(mol,result,scheme=scheme)
        transformed=CalcPiDescriptors(moved,rotated,scheme=scheme)
        assert any(row['pi_n_e'] is not None for row in baseline['atoms'])
        for entity,fields in [('atoms',['pi_n_e','pi_q_e']),('bonds',['pi_bo'])]:
            for x,y in zip(baseline[entity],transformed[entity]):
                for key in fields:
                    if x[key] is None: assert y[key] is None
                    else: assert x[key]==pytest.approx(y[key],abs=2e-6)


def test_matched_temperature_sensitivity_has_units_and_invariant_fields_only():
    mol,r=pi_system(); a=CalcDescriptorReport(mol,r)
    b=CalcDescriptorReport(mol,replace(r,total_energy_hartree=-.9,
        metadata={**r.metadata,'CalculationSettings':CalculationSettings(electronic_temperature=1000.)}))
    out=DescriptorSensitivity(a,temperature_300=a,temperature_1000=b)
    assert out['values']['xtb_energy_ha']['et300_1000_delta']==pytest.approx(.1)
    assert out['units']['xtb_energy_ha']=='Hartree'
    assert out['atom_units']['pi.pi_n_e']=='e'
    assert out['bond_units']['pi.pi_bo']=='dimensionless'
    assert not any(k in out['values'] for k in ['xtb_dip_x_ebohr','homo_cx_a','xtb_homo_idx'])


def test_unavailable_frontier_group_preserves_bond_rows():
    mol,r=pi_system()
    out=CalcFrontierBondDescriptors(mol,replace(r,orbital_occupations=(2.,)*24,number_of_electrons=48))
    virtual=[row for row in out['bonds'] if row['frontier_group']=='lumo']
    assert len(virtual)==mol.GetNumBonds()
    assert all(row['mul_ovlp_e'] is None and row['orb_idx']==[] for row in virtual)
    assert all(row['missing']['mul_ovlp_e']=='frontier_group_unavailable' for row in virtual)


def test_heteroatom_pi_population_has_no_implicit_charge_reference():
    from rdkit.Chem import AllChem
    mol=Chem.AddHs(Chem.MolFromSmiles('c1ccncc1'))
    assert AllChem.EmbedMolecule(mol,randomSeed=75)==0
    # Keep the analytical planar geometry, remove the hydrogen absent on pyridine N.
    _,base=pi_system()
    source_positions=pi_system()[0].GetConformer().GetPositions()
    positions=[*source_positions[:6],*[source_positions[6+i] for i in [0,1,2,4,5]]]
    for i,pos in enumerate(positions):mol.GetConformer().SetAtomPosition(i,pos)
    s=np.asarray(base.overlap_matrix)[:23,:23]
    eig,u=np.linalg.eigh(s);c=(u/np.sqrt(eig))@u.T
    shell_map=base.ao_shell_indices[:23]
    r=replace(base,orbital_energies_ev=base.orbital_energies_ev[:23],
        orbital_occupations=base.orbital_occupations[:23],number_of_molecular_orbitals=23,
        orbital_coefficients=tuple(map(tuple,c)),overlap_matrix=tuple(map(tuple,s)),
        ao_atom_indices=base.ao_atom_indices[:23],ao_shell_indices=shell_map,
        shell_atom_indices=tuple(range(11)),shell_angular_momenta=(1,)*6+(0,)*5,
        ao_p_directions=base.ao_p_directions[:23],partial_charges=(0.,)*11,
        metadata={**base.metadata,**_input_metadata(mol,-1,0,0,_xyz_from_mol(mol,-1))})
    row=CalcPiDescriptors(mol,r)['atoms'][3]
    assert row['pi_n_e']==pytest.approx(1.) and row['pi_q_e'] is None
    assert row['missing']['pi_q_e']=='pi_reference_unspecified'
    assert CalcPiDescriptors(mol,r,pi_ref={3:1.})['atoms'][3]['pi_q_e']==pytest.approx(0.)
