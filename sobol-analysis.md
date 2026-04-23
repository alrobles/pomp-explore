# Sobol Sequence Implementation in `pomp`

This document provides a thorough, cited analysis of how Sobol's algorithm for
generating low-discrepancy sequences in hypercubes is implemented in this
repository.

---

## Overview

The implementation spans three layers:

| Layer | File | Role |
|-------|------|------|
| R public API | `R/sobol_design.R` | Validation, scaling, user-facing wrapper |
| C entrypoint | `src/sobolseq.c` | Core sequence generation |
| Data tables | `src/soboldata.h` | Primitive polynomials and initial direction numbers |

The algorithm follows:

- P. Bratley and B. L. Fox, *Algorithm 659*, ACM Trans. Math. Soft. **14**(1),
  88–100, 1988.
- S. Joe and F. Y. Kuo, *Remark on Algorithm 659*, ACM Trans. Math. Soft.
  **29**(1), 49–57, 2003.

---

## Call Chain

```
sobol_design()          R/sobol_design.R:22
  └── sobol()           R/sobol_design.R:39
        └── .Call(P_sobol_sequence, d, n)    R/sobol_design.R:43
              └── sobol_sequence()           src/sobolseq.c:234
                    ├── nlopt_sobol_create() src/sobolseq.c:188
                    │     └── sobol_init()   src/sobolseq.c:122
                    ├── nlopt_sobol_skip()   src/sobolseq.c:225
                    └── sobol_gen() × (n-1) src/sobolseq.c:96
```

The `P_` prefix on `P_sobol_sequence` comes from `useDynLib(pomp, .registration
= TRUE, .fixes="P_")` (`NAMESPACE:206`, `R/package.R:92`). The native symbol is
registered in `src/init.c:34`.

---

## 1. Primitive Polynomials over GF(2)

**File:** `src/soboldata.h:37–146`

Each Sobol dimension (beyond the first) is associated with a distinct primitive
polynomial over **Z₂** (the binary field). The polynomial
`p(z) = a₀ + a₁z + a₂z² + … + a_{31}z^{31}` is encoded as a 32-bit integer:
bit *i* of `sobol_a[j]` is coefficient `a_i` of the *j*-th polynomial.

```c
// src/soboldata.h:37-40
/* successive primitive binary-coefficient polynomials p(z)
   = a_0 + a_1 z + a_2 z^2 + ... a_31 z^31, where a_i is the
     i-th bit of sobol_a[j] for the j-th polynomial. */
static const uint32_t sobol_a[MAXDIM-1] = { 3, 7, 11, 13, 19, 25, ... };
```

Constants:

```c
// src/soboldata.h:34-35
#define MAXDIM 1111   // maximum supported dimension
#define MAXDEG 12     // maximum polynomial degree
```

---

## 2. Initial Direction Numbers

**File:** `src/soboldata.h:148–end`

For each primitive polynomial of degree `d`, the first `d` direction numbers
`m₀, m₁, …, m_{d-1}` are tabulated in `sobol_minit`:

```c
// src/soboldata.h:148-150
/* starting direction #'s m[i] = sobol_minit[i][j] for i=0..d of the
 * degree-d primitive polynomial sobol_a[j]. */
static const uint32_t sobol_minit[MAXDEG+1][MAXDIM-1] = { ... };
```

These seeds are precomputed values from the Joe & Kuo (2003) tables.

---

## 3. Initialization: `sobol_init`

**File:** `src/sobolseq.c:122–175`

For each dimension `i > 0`:

1. **Compute polynomial degree** `d` by counting the bits of `sobol_a[i-1]`
   (`src/sobolseq.c:136–143`).
2. **Seed** the first `d` direction numbers from the table
   (`src/sobolseq.c:145–147`):
   ```c
   for (j = 0; j < d; ++j)
       sd->m[j][i] = sobol_minit[j][i-1];
   ```
3. **Extend** direction numbers beyond degree `d` using the recurrence
   (`src/sobolseq.c:149–157`):
   ```c
   for (j = d; j < 32; ++j) {
       a = sobol_a[i-1];
       sd->m[j][i] = sd->m[j - d][i];           // base term
       for (k = 0; k < d; ++k) {
           sd->m[j][i] ^= ((a & 1) * sd->m[j-d+k][i]) << (d-k);  // XOR accumulation
           a >>= 1;
       }
   }
   ```
   This is the standard GF(2) direction-number recurrence driven by the
   polynomial coefficients.

Dimension 0 is special-cased: all 32 direction numbers are set to 1
(`src/sobolseq.c:133`).

State arrays `x` (previous point) and `b` (fixed-point position) are zeroed to
start the sequence at index 0 (`src/sobolseq.c:160–172`).

---

## 4. Incremental Point Generation: `sobol_gen`

**File:** `src/sobolseq.c:93–118`

```c
static int sobol_gen (soboldata *sd, double *x)
{
    unsigned c = rightzero32(sd->n++);   // rightmost zero bit = Gray-code step
    for (i = 0; i < sdim; ++i) {
        b = sd->b[i];
        if (b >= c) {
            sd->x[i] ^= sd->m[c][i] << (b - c);
            x[i] = ((double) sd->x[i]) / (1U << (b+1));
        } else {
            sd->x[i] = (sd->x[i] << (c - b)) ^ sd->m[c][i];
            sd->b[i] = c;
            x[i] = ((double) sd->x[i]) / (1U << (c+1));
        }
    }
    return 1;
}
```

**Key ideas:**

- `rightzero32(n)` (`src/sobolseq.c:84–91`) finds the position of the
  least-significant **zero** bit of `n` using a 32-bit magic-number trick from
  Knuth *TAOCP* vol. 4A §7.1.3. This step index is the Gray-code-based
  "which direction number to XOR" selector.
- Each dimension's integer state `x[i]` is XOR-updated with a single
  direction-number entry `m[c][i]`, appropriately shifted by the fixed-point
  position `b[i]`. This corresponds to flipping exactly one bit of the Van der
  Corput/Sobol lattice.
- The floating-point output is obtained by scaling: `x[i] / 2^(b+1)`.
- The generation runs in **O(dim)** per point.

---

## 5. Sequence Start-Point Skipping

**File:** `src/sobolseq.c:221–231`

Joe & Kuo (2003) recommend starting the sequence at the largest power of 2
strictly less than `n`, so the sequence is better stratified:

```c
static void nlopt_sobol_skip(nlopt_sobol s, unsigned n, double *x)
{
    unsigned k = 1;
    while (k*2 < n) k *= 2;   // k = largest power of 2 < n
    while (k-- > 0) sobol_gen(s, x);
}
```

This was added in `pomp` v1.7.6 (`inst/NEWS:1523–1525`).

---

## 6. R-Level Entry Point: `sobol_sequence`

**File:** `src/sobolseq.c:234–249`

```c
SEXP sobol_sequence (SEXP dim, SEXP length)
{
    nlopt_sobol s = nlopt_sobol_create((unsigned int) d);
    if (s == 0) err("dimension is too high");
    PROTECT(data = allocMatrix(REALSXP, d, n));
    nlopt_sobol_skip(s, n, dp);               // skip warm-up points
    for (k = 1; k < n; k++) sobol_gen(s, dp + k*d);
    nlopt_sobol_destroy(s);
    UNPROTECT(1);
    return data;
}
```

The matrix returned is `d × n` (column-major in R, so each column is a point
in `ℝ^d`).

---

## 7. R Wrapper: `sobol_design` / `sobol`

**File:** `R/sobol_design.R:22–53`

`sobol_design` validates that `lower` and `upper` are named vectors of equal
length with matching names, then calls the internal `sobol()` helper.  
`sobol()` calls `.Call(P_sobol_sequence, d, n)` and affinely maps each column
from `(0,1)` to the user-supplied hyper-rectangle:

```r
vars[[k]][1L] + (vars[[k]][2L] - vars[[k]][1L]) * x[k,]
```

Limits checked in R (`R/sobol_design.R:41`):

- `n > 2^30` → `"too many points requested"` (tested in `tests/sobol.R:11`).
- `d > MAXDIM` → `"dimension is too high"` (tested in `tests/sobol.R:14`).

---

## 8. Reusable Practices

| Practice | Where |
|----------|-------|
| Separate immutable tables from algorithm code | `src/soboldata.h` vs `src/sobolseq.c` |
| Bit-level GF(2) arithmetic (XOR/shift) avoids heavy library dependencies | `src/sobolseq.c:153–155`, `108`, `112` |
| Stateful O(dim) incremental generator using `(x, b)` pair | `src/sobolseq.c:65–72`, `105–116` |
| Documented, standards-cited warm-start skip policy | `src/sobolseq.c:221–231` |
| Layered API: R validation+scaling → C generation → explicit symbol registration | `R/sobol_design.R`, `src/init.c:34`, `NAMESPACE:206` |

---

## References

- P. Bratley and B. L. Fox. Algorithm 659 Implementing Sobol's quasirandom
  sequence generator. *ACM Transactions on Mathematical Software* **14**,
  88–100, 1988. <https://doi.org/10.1145/42288.214372>
- S. Joe and F. Y. Kuo. Remark on algorithm 659: Implementing Sobol'
  quasirandom sequence generator. *ACM Transactions on Mathematical Software*
  **29**, 49–57, 2003. <https://doi.org/10.1145/641876.641879>
- D. E. Knuth. *The Art of Computer Programming*, vol. 4A, §7.1.3.
