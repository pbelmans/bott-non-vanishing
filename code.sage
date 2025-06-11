def levi_type(G, k):
    r"""
    Compute the Dynkin type of the Levi factor
    """
    G = CartanType(G)
    n = G.rank()

    assert 1 <= k and k <= n, "k is invalid"

    if G.type() == "A":
        return CartanType(["A", k - 1], ["A", n - k])

    elif G.type() == "B":
        if n == k:  # because Sage doesn't do B_0
            return CartanType(["A", n - 1])
        return CartanType(["A", k - 1], ["B", n - k])

    elif G.type() == "C":
        if n == k:  # because Sage doesn't do C_0
            return CartanType(["A", n - 1])
        return CartanType(["A", k - 1], ["C", n - k])

    elif G.type() == "D":
        if n == 4:
            if k in [1, 3, 4]:
                return CartanType(["A", 3])
            else:
                return CartanType(["A", 1], ["A", 1], ["A", 1])
        else:
            if k == n - 1 or k == n:
                return CartanType(["A", n - 1])
            return CartanType(["A", k - 1], ["D", n - k])

    elif G.type() == "E":
        if k == 1:
            return CartanType(["D", n - 1])
        elif k == 2:
            return CartanType(["A", n - 1])
        elif k == 3:
            return CartanType(["A", 1], ["A", n - 2])
        elif k == 4:
            return CartanType(["A", 2], ["A", 1], ["A", n - 4])
        elif k == 5:
            return CartanType(["A", 4], ["A", n - 5])
        elif k == 6:
            return CartanType(["D", 5], ["A", n - 6])
        elif k == 7:
            return CartanType(["E", 6], ["A", n - 7])
        elif k == 8:
            return CartanType(["E", 7])

    elif G.type() == "F":
        if k == 1:
            return CartanType("C3")
        elif k == 2:
            return CartanType("A1", "A2")
        elif k == 3:
            return CartanType("A2", "A1")
        else:
            return CartanType("B3")

    elif G.type() == "G":
        return CartanType("A1")


def relabeling(G, k):
    G = CartanType(G)
    n = G.rank()

    # these are exceptional relabelings
    if G.type() == "D" and n == 4:
        if k == 1:
            return {1: 3, 2: 2, 3: 4}
        elif k == 3:
            return {1: 1, 2: 2, 3: 4}
        elif k == 4:
            return {1: 1, 2: 2, 3: 3}

    if G.type() in ["A", "B", "C", "D"]:
        return {i: i for i in range(1, k)} | {i: i + 1 for i in range(k, n)}

    elif G.type() == "E":
        if k == 1:
            return {i: n + 1 - i for i in range(1, n)}
        elif k == 2:
            return {i: i for i in range(1, k)} | {i: i + 1 for i in range(k, n)}
        elif k == 3:
            return {1: 1, 2: 2} | {i: i + 1 for i in range(3, n)}
        elif k == 4:
            return {1: 1, 2: 3, 3: 2} | {i: i + 1 for i in range(4, n)}
        elif k == 5:
            return {1: 1, 2: 3, 3: 4, 4: 2} | {i: i + 1 for i in range(5, n)}
        elif k == 6:
            return {1: 1, 2: 3, 3: 4, 4: 5, 5: 2} | {i: i + 1 for i in range(6, n)}
        elif k == 7:
            return {i: i for i in range(1, k)} | {i: i + 1 for i in range(k, n)}
        elif k == 8:
            return {i: i for i in range(1, n)}

    elif G.type() == "F":
        if k == 1:
            return {1: 4, 2: 3, 3: 2}
        elif k == 2:
            return {1: 1, 2: 3, 3: 4}
        elif k == 3:
            return {1: 1, 2: 2, 3: 4}
        else:
            return {1: 1, 2: 2, 3: 3}

    elif G.type() == "G":
        return {1: 3 - k}


def restriction(GP):
    LG = RootSystem(CartanType(GP[0])).weight_lattice()
    k = GP[1]

    # constructing the restriction matrix: using simple roots
    M = matrix(
        [alpha.dense_coefficient_list() for alpha in LG.simple_roots()]
    ).transpose()
    # using the ith fundamental weight for the basis of this subspace
    M[:, k - 1] = matrix(LG.fundamental_weight(k).dense_coefficient_list()).transpose()

    return M.inverse().row(k - 1)


# combine every summand with the line bundle
def combine(GP, summands, line_bundle):
    LG = RootSystem(CartanType(GP[0])).weight_lattice()

    # the Levi type and its weight lattice
    L = levi_type(GP[0], GP[1])
    LL = RootSystem(L).weight_lattice()

    result = []

    for summand in summands:
        # replacing baromega with omega
        A = RootSystem(L).ambient_space()
        w = A(list(summand._monomial_coefficients.keys())[0])

        summand = sum(
            [
                LL.fundamental_weight(i) * int(A.simple_coroot(i).dot_product(w))
                for i in LL.index_set()
            ]
        )
        lift = sum(
            [
                LG.fundamental_weight(relabeling(GP[0], GP[1])[i])
                * summand.coefficient(i)
                for i in L.index_set()
            ]
        )

        # untwisting the replacement
        twist = vector(lift.dense_coefficient_list()).dot_product(restriction(GP))

        # contribution from the tensor power
        twist = twist - line_bundle
        assert twist.is_integer(), "twist is not an integer"
        twist = int(twist)

        # end result
        result.append(lift - twist * LG.fundamental_weight(GP[1]))

    return result


def twist_by_line_bundle(GP, w, i):
    r"""
    Twist the tensor product of the vector bundle described by the weight `w` by `O(i)`
    """
    LG = RootSystem(CartanType(GP[0])).weight_lattice()
    LG.print_options(prefix="w")

    return w + i * LG.fundamental_weights()[GP[1]]


def exterior_power(GP, w, p):
    r"""
    Compute the `p`th exterior power of vector bundle described by the weight `w`
    """
    LG = RootSystem(CartanType(GP[0])).weight_lattice()
    LG.print_options(prefix="w")

    if w not in LG:
        w = sum([a * b for (a, b) in zip(w, LG.fundamental_weights())])

    # projection on semisimple part: element of weight lattice of Levi factor
    L = levi_type(CartanType(GP[0]), GP[1])
    LL = RootSystem(L).weight_lattice()
    wA = sum(
        [
            LL.fundamental_weight(i) * w.coefficient(relabeling(GP[0], GP[1])[i])
            for i in L.index_set()
        ]
    )

    # projection on torus: scalar
    wB = sum(
        [
            w.coefficient(i) * restriction(GP)[i - 1]
            for i in CartanType(GP[0]).index_set()
        ]
    )

    # summands for the exterior power restricted to A
    RL = WeylCharacterRing(L)
    # compute exterior power on semisimple part
    exterior_power = RL(wA).exterior_power(p)

    summands = [
        summand
        for summand in exterior_power.monomials()
        for i in range(exterior_power.multiplicity(summand))
    ]

    # putting in the tensor power of the line bundle
    line_bundle = p * wB

    return combine(GP, summands, line_bundle)


def rank(GP, w):
    r"""
    Compute the rank of the vector bundle described by the weight `w`
    """
    LG = RootSystem(CartanType(GP[0])).weight_lattice()
    LG.print_options(prefix="w")

    if w not in LG:
        w = sum([a * b for (a, b) in zip(w, LG.fundamental_weights())])
    # the Levi type and its weight lattice
    L = levi_type(GP[0], GP[1])
    LL = RootSystem(L).weight_lattice()

    # projecting onto the Levi
    vA = sum(
        [
            LL.fundamental_weight(i) * w.coefficient(relabeling(GP[0], GP[1])[i])
            for i in L.index_set()
        ]
    )
    return WeylCharacterRing(L)(vA).degree()


def is_singular(w):
    r"""
    Determine whether the weight `w` is a singular weight
    """
    w = w.to_dominant_chamber()
    return w in w.simple_reflections()


def cohomology(GP, w):
    r"""
    Compute the cohomology of the vector bundle described by the weight `w`

    If the weight is singular, return None.
    If the weight is regular, return the degree of non-zero cohomology,
    the weight of the representation, its dimension,
    and the weight of the bundle.
    """
    LG = RootSystem(CartanType(GP[0])).weight_lattice()
    LG.print_options(prefix="w")

    if w not in LG:
        w = sum([a * b for (a, b) in zip(w, LG.fundamental_weights())])

    if is_singular(w + LG.rho()):
        return None

    dominant, word = (w + LG.rho()).to_dominant_chamber(reduced_word=True)
    dominant = dominant - LG.rho()
    degree = len(word)
    dimension = WeylDim(
        RootSystem(CartanType(GP[0])), dominant.dense_coefficient_list()
    )

    return (degree, dominant, dimension, w)
