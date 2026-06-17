import pytest

from rdkit import Chem

LOOP_COUNT = 1000000


def test_set_prop(loop_count=LOOP_COUNT):
    mol = Chem.MolFromSmiles('n1ccccc1')
    atom = mol.GetAtomWithIdx(0)
    for i in range(loop_count):
        atom.SetProp('test_prop', 'test_val')
    assert atom.GetProp('test_prop') == 'test_val'


if __name__ == '__main__':
    import sys
    sys.exit(pytest.main([__file__]))
