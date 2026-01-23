"""Code to generate orthogonal set partitions with a uniform block size."""
# Copyright (C) 2023 Anton Pirogov, licensed under the MIT License
from itertools import product, repeat
import numpy as np

def prime_mols(n: int):
    """Generate complete set of MOLS for prime number n."""
    return [[[(k*i + j) % n for j in range(n)] for i in range(n)] for k in range(1, n)]

def mols_to_mops(mols):
    """Convert MOLS into orthogonal set partitions via orthogonal arrays.

    Returns n+2 set partitions of n^2 elements from a set of n MOLS.
    The returned list will always start with the trivial partition
    (i.e., counting up in sequential order, separated into groups).

    See https://en.wikipedia.org/wiki/Mutually_orthogonal_Latin_squares#Nets
    """
    oa = [  # convert to OA
        [r, c] + [mols[i][r][c] for i in range(len(mols))]
        for r in range(len(mols[0]))
        for c in range(len(mols[0][0]))
    ]
    num_mols, ls_sz, oa_rows = len(mols), len(mols[0]), len(oa)
    rc_j = [ # partitions based on rows/cols of OA
        [
            [ ls_sz * oa[k][0] + oa[k][1]
              for k in range(oa_rows) if oa[k][row_or_col] == j
            ] for j in range(ls_sz)
        ] for row_or_col in range(2)
    ]
    l_ij = [ # partitions induced by latin square columns
        [
            [
                ls_sz * oa[k][0] + oa[k][1]
                for k in range(oa_rows) if oa[k][2 + i] == j
            ] for j in range(ls_sz)
        ] for i in range(num_mols)
    ]
    return rc_j + l_ij

def base_rgdd(n: int):
    """Return a base RGDD for prime number n.

    The RGDD consists of n mutually orthogonal, non-trivial partitions
    of n^2 elements, each with uniform block size n.
    """
    return np.array(mols_to_mops(prime_mols(n))[1:])

# ----

def base_rgdd2(n: int):
    """Like base_rgdd, but direct and avoiding MOLS construction."""
    std = np.arange(n**2).reshape(n, n)
    return np.array([
        np.transpose(np.array([np.roll(std[i], k*i) for i in range(n)]))
        for k in range(1, n+1)
    ])

def get_tree(b: int, i: int):
    """Get standard index tree for b^i elements."""
    return np.arange(b**i).reshape(*repeat(b,i))

def to_partition(tree):
    """Convert an index tree into a set partition."""
    return tree.reshape(-1, tree.shape[-1])

def permutate_subtrees(tree, lvl: int, perm_mat):
    """Apply a permutation (given as a matrix) to all tree nodes at some level."""
    dim, deg = tree.ndim, len(tree)
    old_shape = list(repeat(deg, dim))# save old shape
    # build numpy indexing expression for subtree permutation
    perm_expr = list(repeat(slice(None),lvl)) + [perm_mat]
    # target shape = "flatten" the currently focused tree level
    flat_lvl = list(repeat(deg, lvl)) + [deg**2] + list(repeat(deg, dim - lvl - 2))
    # apply permutation in terms of numpy operations
    return tree.reshape(*flat_lvl)[*perm_expr].reshape(old_shape)

def permutate_tree(tree, perm_word):
    """Apply a sequence of permutations to an index tree."""
    assert tree.ndim == len(perm_word)+1
    for lv, lv_perm in enumerate(perm_word):
        tree = permutate_subtrees(tree, lv, lv_perm)
    return tree

def exp_partitions(perm_mats, i: int):
    """Lift a set of non-trivial orthogonal partitions on b^2 points to b^i points."""
    assert i > 1
    assert perm_mats[0].ndim == 2 and perm_mats[0].shape[0] == perm_mats[0].shape[1]
    block_size = perm_mats[0].shape[1]
    std_tree = get_tree(block_size, i)

    perms = [p.reshape(-1).tolist() for p in perm_mats]
    perm_words = product(*repeat(perms, i-1))
    return [to_partition(permutate_tree(std_tree, word)) for word in perm_words]

def gen_partitions(b: int, i: int):
    """For a prime b, generate b^(i-1) orthogonal partitions of b^i points."""
    return exp_partitions(base_rgdd2(b), i)

# ----

def partitions_orthogonal(p1, p2) -> bool:
    """Check whether two partitions are orthogonal."""
    return all(map(lambda b: len(set(b[0]).intersection(set(b[1]))) <= 1, product(p1, p2)))

def is_uniform_partition(p, b_size=None) -> bool:
    """Check that p is a partition of numbers 0 upto a value with uniform block size."""
    b_size = b_size or len(p[0])
    ns = p.reshape(-1).tolist()
    is_partition = set(ns) == set(range(len(ns)))  # is partition on 0..max
    is_uniform = all(map(lambda b: len(b) == b_size, p))  # uniform block size
    return is_partition and is_uniform

def sym_pairs(*vs):
    """Return all unordered pairs (i.e. outputs 1,2 but not 2,1)."""
    s = list(vs)
    return ((s[i], s[j]) for i in range(len(s)) for j in range(i+1, len(s)))

def check_partitions(ps):
    """Check that the computed partitions are uniform and mutually orthogonal."""
    b_size = len(ps[0][0])
    assert all(map(lambda p: is_uniform_partition(p, b_size), ps))
    for i, j in sym_pairs(*range(len(ps))):
        assert partitions_orthogonal(ps[i], ps[j])
