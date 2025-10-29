load("code.sage")


def E(D):
    """
    Return the highest weight of the vector bundle `E` of an adjoint variety

    - `D` is the Dynkin type
    """
    G = D[0]
    n = int(D[1:])

    if G == "B":
        if n == 3:
            return (1, -1, 2)
        if n >= 4:
            return (1, -1, 1)
    if G == "D":
        if n == 4:
            return (1, -1, 1, 1)
        elif n >= 5:
            return (1, -1, 1)
    if G == "E" and n == 6:
        return (0, -1, 0, 1)
    if G == "E" and n == 7:
        return (-1, 0, 1)
    if G == "E" and n == 8:
        return (0, 0, 0, 0, 0, 0, 1, -1)
    if G == "F" and n == 4:
        return (-1, 1)
    if G == "G" and n == 2:
        return (3, -1)

    raise Exception


def verify(D, k, q):
    r"""
    Verify that there is non-zero cohomology in degree 1 for the given setting

    `D` is the Dynkin type
    `k` is the marked vertex
    `q` is the exterior power of the tangent bundle
    """
    GP = (D, k)

    print(f"Verifying the claim for adjoint of type {D}")
    print(f"using the {q}th exterior power of the tangent bundle twisted by -1")

    summands = [
        twist_by_line_bundle(GP, summand, -1) for summand in exterior_power(GP, E(D), q)
    ] + exterior_power(GP, E(D), q - 1)

    cohomologies = [cohomology(GP, summand) for summand in summands]
    T = [("weight of bundle", "rank", "H^i", "weight of cohomology", "dimension")]
    for summand, H in zip(summands, cohomologies):
        if H is None:
            T.append((summand, rank(GP, summand), "", "", ""))
        else:
            T.append((summand, rank(GP, summand), H[0], H[1], H[2]))

    print()
    print(table(T, header_row=True))
    print()

    cohomologies = list(filter(None, cohomologies))

    # type D_4 gives 3 representations in degree 1
    if D == "D4":
        weights = [H[1] for H in cohomologies if H[0] == 1]
        # too many components in degree 1 to be cancelled by degree 0
        # this will not happen for our choice
        if len(weights) < len([H for H in cohomologies if H[0] == 0]):
            raise Exception

        print(f"Cohomologies in degree 1 are")
        for H in cohomologies:
            if H[0] == 1:
                print(f"- vector bundle of weight {H[3]} giving {H[1]}")
        print()

        return True

    # look for the (unique) representation in degree 1 (if not of type D_4)
    weight = None
    for H in cohomologies:
        if H[0] == 0:
            continue

        # check there is a unique one
        if weight:
            raise Exception

        weight = (H[1], H[3])

    # verify that the representation is not cancelled by degree 0
    for H in cohomologies:
        if H[0] == 1:
            continue

        if H[1] == weight[0]:
            raise Exception

    print(f"Cohomology in degree 1 is")
    print(f"- vector bundle of weight {weight[1]} giving {weight[0]}")
    print()
    print()

    return True


# types B and D
for n in range(3, 10):
    verify(f"B{n}", 2, 3)
for n in range(4, 10):
    verify(f"D{n}", 2, 3)

# exceptional cases
verify("E6", 2, 5)
verify("E7", 1, 7)
verify("E8", 8, 11)
verify("F4", 1, 4)
verify("G2", 2, 2)
