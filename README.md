# Code accompanying "Failure of Bott vanishing for (co)adjoint partial flag varieties"

This repository contains the code supporting the computations in the paper
"Failure of Bott vanishing for (co)adjoint partial flag varieties",
by Pieter Belmans,
see [`arXiv:2506.09811`](https://arxiv.org/abs/2506.09811).

It implements a computation of exterior powers of completely reducible homogeneous vector bundles
and the Borel-Weil-Bott theorem, to compute the cohomology of twists of exterior powers of the tangent bundle
on certain generalized Grassmannians.

In order to cite this code, please use

```bibtex
@software{bott-non-vanishing,
  author = {Belmans, Pieter},
  title  = {Code accompanying ``Failure of Bott vanishing for (co)adjoint partial flag varieties''},
  url    = {https://github.com/pbelmans/bott-non-vanishing},
  year   = {2025},
}
```

## `adjoint.sage`

By running

```bash
sage adjoint.sage
```

one computes the results reported in the paper.

The code computes the indicated exterior power of
the semisimplification of the tangent bundle,
twisted by `O(-1)`, and verifies that

* there is cohomology in degree 1 which cannot be cancelled

The code was written and tested using **SageMath 10.5**.
It ran in 5 hours, on a MacBook Pro M4 with 32 GB of memory
(memory requirements are the main bottleneck).

The only time-consuming computation is that for E8,
all other types combined finishes in less than a minute.
The output can be found in [`output-adjoint.txt`](output-adjoint.txt).

## `coadjoint.sage`

The code uses the semisimplification of the tangent bundle,
twisted by `O(-1)`, and verifies that

* there is cohomology in degree 1 which cannot be cancelled

The code was written and tested using **SageMath 10.5**.
It runs in less than 3 seconds on a MacBook Pro M4 with 32 GB of memory.
The output can be found in [`output-coadjoint.txt`](output-coadjoint.txt).
