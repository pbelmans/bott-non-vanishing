load("code.sage")


def verify(D, k, A, B):
    r"""
    Verify that there is non-zero cohomology in degree 1 for the given setting

    `D` is the Dynkin type
    `k` is the marked vertex
    `A` is the subbundle of T_X(-1)
    `B` is the quotient bundle of T_X(-1)
    """
    GP = (D, k)
    LG = RootSystem(CartanType(GP[0])).weight_lattice()
    LG.print_options(prefix="w")

    if not cohomology(GP, A)[0] == 1 or cohomology(GP, B) != None:
        raise Exception

    print(f"Verifying the claim for adjoint of type {D}")
    print("using the tangent bundle twisted by -1")

    T = [("weight of bundle", "rank", "H^i", "weight of cohomology", "dimension")]
    cohomologies = [cohomology(GP, w) for w in [A, B]]
    for summand, H in zip([A, B], cohomologies):
        w = sum([a * b for (a, b) in zip(summand, LG.fundamental_weights())])

        if H is None:
            T.append((w, rank(GP, summand), "", "", ""))
        else:
            T.append((w, rank(GP, summand), H[0], H[1], H[2]))

    print()
    print(table(T, header_row=True))
    print()

    H = cohomology(GP, A)

    print("Cohomology in degree 1 is")
    print(f"- vector bundle of weight {H[3]} giving {H[1]}")
    print()
    print()

    return True


for n in range(3, 10):
    verify(f"C{n}", 2, (1, -2, 1), (2, -1))

verify("F4", 4, (0, 0, 1, -2), (1, 0, 0, -1))
