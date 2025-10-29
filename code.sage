def levi_type(G, k):
    r"""
    Compute the Dynkin type of the Levi factor `L`

    - `G` is the Dynkin type of the ambient simple algebraic group
    - `k` is the index of the vertex defining the maximal parabolic subgroup
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
    r"""
    Compute how fundamental weights of `G` restrict to those of `L`

    - `G` is the Dynkin type of the ambient simple algebraic group
    - `k` is the index of the vertex defining the maximal parabolic subgroup
    """
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
    r"""
    Determine the restriction of fundamental weights

    - `GP` is a pair `(D, k)` of a Dynkin type and an index `k` encoding `G/P`
    """
    (D, k) = GP
    LG = RootSystem(CartanType(D)).weight_lattice()

    # constructing the restriction matrix: using simple roots
    M = matrix(
        [alpha.dense_coefficient_list() for alpha in LG.simple_roots()]
    ).transpose()
    # using the ith fundamental weight for the basis of this subspace
    M[:, k - 1] = matrix(LG.fundamental_weight(k).dense_coefficient_list()).transpose()

    return M.inverse().row(k - 1)


# combine every summand with the line bundle
def combine(GP, summands, line_bundle):
    r"""
    Combine the computation in the center and derived subgroup

    - `GP` is a pair `(D, k)` of a Dynkin type and an index `k` encoding `G/P`
    - `summands` is a list of fundamental weights for the derived subgroup
    - `line_bundle` is the twist defined by the center
    """
    (D, k) = GP
    LG = RootSystem(CartanType(D)).weight_lattice()

    # the Levi type and its weight lattice
    L = levi_type(D, k)
    LL = RootSystem(L).weight_lattice()

    result = []

    for summand in summands:
        # lifting to weight lattice of the Levi
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
                LG.fundamental_weight(relabeling(D, k)[i])
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
        result.append(lift - twist * LG.fundamental_weight(k))

    return result


def twist_by_line_bundle(GP, weight, i):
    r"""
    Twist the tensor product of the vector bundle described by the weight `weight` by `O(i)`

    - `GP` is a pair `(D, k)` of a Dynkin type and an index `k` encoding `G/P`
    - `weight` is a fundamental `L`-dominant weight
    - `i` the twist
    """
    (D, k) = GP
    LG = RootSystem(CartanType(D)).weight_lattice()
    LG.print_options(prefix="w")

    return weight + i * LG.fundamental_weights()[k]


def exterior_power(GP, weight, p):
    r"""
    Compute the `p`th exterior power of vector bundle described by the weight `weight`

    - `GP` is a pair `(D, k)` of a Dynkin type and an index `k` encoding `G/P`
    - `weight` is a fundamental `L`-dominant weight
    - `p` is the exponent of the exterior power
    """
    (D, k) = GP
    LG = RootSystem(CartanType(D)).weight_lattice()
    LG.print_options(prefix="w")

    # coerce to weight lattice if necessary
    if weight not in LG:
        weight = sum([a * b for (a, b) in zip(weight, LG.fundamental_weights())])

    # projection on semisimple part: element of weight lattice of Levi factor
    L = levi_type(CartanType(D), k)
    LL = RootSystem(L).weight_lattice()
    wA = sum(
        [
            LL.fundamental_weight(i) * weight.coefficient(relabeling(D, k)[i])
            for i in L.index_set()
        ]
    )

    # projection on torus: scalar
    wB = sum(
        [
            weight.coefficient(i) * restriction(GP)[i - 1]
            for i in CartanType(D).index_set()
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


def rank(GP, weight):
    r"""
    Compute the rank of the vector bundle described by the weight `weight`

    - `GP` is a pair `(D, k)` of a Dynkin type and an index `k` encoding `G/P`
    - `weight` is a fundamental `L`-dominant weight
    """
    (D, k) = GP
    LG = RootSystem(CartanType(D)).weight_lattice()
    LG.print_options(prefix="w")

    # coerce to weight lattice if necessary
    if weight not in LG:
        weight = sum([a * b for (a, b) in zip(weight, LG.fundamental_weights())])

    # the Levi type and its weight lattice
    L = levi_type(D, k)
    LL = RootSystem(L).weight_lattice()

    # projecting onto the Levi
    vA = sum(
        [
            LL.fundamental_weight(i) * weight.coefficient(relabeling(D, k)[i])
            for i in L.index_set()
        ]
    )
    return WeylCharacterRing(L)(vA).degree()


def is_singular(weight):
    r"""
    Determine whether the weight `weight` is a singular weight
    """
    weight = weight.to_dominant_chamber()
    return weight in weight.simple_reflections()


def cohomology(GP, weight):
    r"""
    Compute the cohomology of the vector bundle described by the weight `weight`

    This encodes the Borel-Weil-Bott theorem.

    If the weight is singular, return None.
    If the weight is regular, return the degree of non-zero cohomology,
    the weight of the representation, its dimension,
    and the weight of the bundle.

    - `GP` is a pair `(D, k)` of a Dynkin type and an index `k` encoding `G/P`
    - `weight` is a fundamental `L`-dominant weight
    """
    (D, _) = GP
    LG = RootSystem(CartanType(D)).weight_lattice()
    LG.print_options(prefix="w")

    # coerce to weight lattice if necessary
    if weight not in LG:
        weight = sum([a * b for (a, b) in zip(weight, LG.fundamental_weights())])

    # the singular case
    if is_singular(weight + LG.rho()):
        return None

    # the regular case
    dominant, word = (weight + LG.rho()).to_dominant_chamber(reduced_word=True)
    dominant = dominant - LG.rho()
    degree = len(word)
    dimension = WeylDim(
        RootSystem(CartanType(D)), dominant.dense_coefficient_list()
    )

    return (degree, dominant, dimension, weight)
