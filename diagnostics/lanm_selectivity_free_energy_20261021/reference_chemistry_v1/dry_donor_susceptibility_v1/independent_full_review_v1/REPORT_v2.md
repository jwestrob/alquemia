# Independent full10 dry-donor susceptibility review

Full10 is numerically UNQUALIFIED:14/16 checks pass; both Dy odd-energy/zero-moment checks fail. The La scout receipt remains unchanged. Zero native calls by this review.

| Quantity | h, au | Model χE | Quantum χE | Quantum χμ | Model−quantum χE | Model−quantum χμ |
|---|---:|---:|---:|---:|---:|---:|
| La | .001 | 54.125748541 | 47.882389600 | 47.865829688 | 6.243358941 | 6.259918852 |
| La | .0005 | 54.125716669 | 47.881020000 | 47.862562672 | 6.244696669 | 6.263153996 |
| Dy | .001 | 49.102971785 | 55.102519190 | 55.119005865 | -5.999547405 | -6.016034079 |
| Dy | .0005 | 49.102939913 | 54.983321680 | 55.111944599 | -5.880381767 | -6.009004686 |
| Dy−La | .001 | -5.022776755 | 7.220129590 | 7.253176176 | -12.242906345 | -12.275952932 |
| Dy−La | .0005 | -5.022776755 | 7.102301680 | 7.249381927 | -12.125078435 | -12.272158682 |

Both Dy odd-slope residuals, at h=.001/.0005 au, are 0.00017150721560817942 / 0.0002506522106081794 ea0, exceeding the unchanged1e-4 threshold. These raw susceptibility comparisons are retained descriptively; the full block did not qualify. No averaging, extrapolation, changed field step or rerun is authorized.

Susceptibilities are in a0³. Both estimators/steps, raw tokens and propagated bounds are retained.

| Central projection, ea0 | Model μE at half step | Quantum zero μ |
|---|---:|---:|
| La | 0.114374815 | 0.580438037 |
| Dy | 0.100683834 | 0.635663281 |
| Dy−La | -0.013690981 | 0.055225244 |

All nuclear/electronic/total dipole, laboratory-origin, state/basis, field and electronic/final/D4 closures pass. Fixed-geometry D4 is field invariant. Full scheduler cost including ten native jobs and two collectors: 6347 CPU-s; exact per-job breakdown retained. Numerical consistency does not establish physical reference accuracy; failed native printed SCF criteria and frozen-f/basis/global-projection limitations remain. No unique mechanism or parameter repair is inferred.
