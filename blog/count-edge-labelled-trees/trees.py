# Naive implementation of a natural edge-labeled tree encoding into numbers.
# Copyright (C) 2024 Anton Pirogov, Licensed under MIT License

def primes(n):
    """Generate all prime numbers up to some value n."""
    sieve = [True] * (n+1)
    for p in range(2, n+1):
        if (sieve[p]):
            yield p
            for i in range(p, n+1, p):
                sieve[i] = False

# Primes up to some arbitrary constant, in-order + inverse lookup table
Pi = dict(map(lambda x: (x[0]+1, x[1]), enumerate(primes(300000))))
Pi[0] = 1
Np = {v: k for k, v in Pi.items()}

def is_prime(n: int) -> bool:
    """Check whether a number is prime."""
    return n in Np

def factorize(n: int):
    """Compute prime factorization of n."""
    ret, cur, i = {}, 0, 1
    while n > 1:
        p = Pi[i]
        if n % p == 0:
            n /= p
            if i != cur:
                cur = i
                ret[p] = 0
            ret[p] += 1
        else:
            i += 1
    return ret

def fac_tree(n):
    """Compute factorization tree of n."""
    return { p: fac_tree(exp) for p, exp in factorize(n).items() }

def fac_to_el_tree(t):
    """Decode an edge-labeled tree from a factorization tree."""
    if not len(t):
        return []

    ps = list(t.keys())
    ns = [Np[ps[0]]]
    for i in range(1, len(ps)):
        ns.append(Np[ps[i]]-Np[ps[i-1]])

    return [(ns[i], fac_to_el_tree(t[ps[i]])) for i in range(len(ps))]

def to_tree(n: int):
    """Get the edge-labeled tree associated with number n.

    A tree node is represented as an ordered list of pairs (l, s)
    where l is an edge label value and s is another tree node.
    """
    return fac_to_el_tree(fac_tree(n))

# TODO: Exercise for the reader:
# Implement reverse mapping (tree to number).

def edge_sum(t) -> int:
    """Return sum of all edges in the edge-labelled tree."""
    es = map(lambda x: x[0], t)
    ts = map(lambda x: x[1], t)
    return sum(es) + sum(map(edge_sum, ts))

def height(t) -> int:
    if not t:
        return 0
    return 1+max(map(lambda x: height(x[1]), t))

print("N Height Sum Tree")
for n in range(1, 100):
    t = to_tree(n)
    print(n, height(t), edge_sum(t), t)
