# champCoeffsTable

Optimal interface coefficients for the CHAMP conjugate heat transfer scheme
(Meng et al. 2017), looked up from a precomputed table. Currently: the
second-order (BDF2) scheme.

## Usage (MATLAB)

```matlab
[rL,rR] = get_r_opt( theta,beta,lambda );          % r in [0,1]
[pL,pR] = get_p_opt( theta,beta,lambda );          % Robin weights p = r/(1-r)

[rL,rR,kappa,A,info] = get_r_opt( theta,beta,lambda,method );
```

* `theta = K_L/K_R`, `beta = D_L/D_R`, `lambda = dt*D_L/h^2` (scalars or arrays of equal size)
* interface condition on each side: `(1-r)*N[T] + r*D[T]`; `r = 0` Neumann, `r = 1`
  Dirichlet; `p = r/(1-r)` (`p = Inf` is the Dirichlet limit)
* `method`: `'linear'` (default) or `'nearest'`
* optional outputs:
  * `kappa`, `A`: wavenumbers and complex amplification factors of the (up to 4) tracked
    peaks of |A(kappa)| of the chosen solution, largest |A| first (estimates from the
    nearest table node)
  * `info.Amax`: verified max|A| over all roots and wavenumbers; `info.Astar`: the minimax
    value; `info.family`, `info.families`: which family of optima was chosen;
    `info.active`: peaks in the minimax active set; `info.method`: per query

The table covers `theta, beta` in `[1e-6, 1e6]`, `lambda` in `[1, 1e8]`. Outside that
range the values of the nearest table node (on the boundary of the table) are returned,
with `info.method = 'nearestExtrap'` and a warning. `NaN` is returned where no stable
coefficients exist at the table node (max|A| > 1 for every family; `info.method = 'none'`)
and for non-positive or non-finite input.

## How the lookup works

The table holds several families of locally optimal (minimax) coefficients. Each family is
interpolated on its own (trilinear in log10 of theta, beta, lambda; r in the measure
log(r+0.01)); a family is used only where it is stable at the surrounding table nodes. The
best family is then chosen: smallest minimax value, with near-ties decided by the smallest
verified max|A| and continuity with the neighbouring nodes. Interpolated coefficients are
not themselves verified; with `'nearest'` every returned value has been verified (all
roots, all wavenumbers) at the table node. Next to the stability limit, where no family can
be interpolated, the nearest-node value is returned.

## Table

`champ_order2.h5`: 0.25-decade grid in theta, beta, lambda. Per family (`/A`, `/Ap`, `/C`,
`/Cp`): `rL`, `rR`, `Atracked` (minimax value), `Amax` (verified), `kappaMax`, `status`
(1 optimal, 3 repaired, -1 unstable, 0 none), `peakKappa`, `peakAreal`, `peakAimag`,
`peakActive`. Global: `best` (chosen family per node), `lambdaStar` / `lambdaStarNext`
(largest stable lambda per theta, beta).
