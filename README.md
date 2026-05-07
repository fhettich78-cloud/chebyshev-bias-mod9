# Chebyshev Bias modulo 9

Numerical study of the prime number race between quadratic residues and
quadratic non-residues modulo 9, up to N = 10⁷.

## Background

Chebyshev observed in 1853 that primes seem to prefer the residue class
3 (mod 4) over 1 (mod 4). This *Chebyshev bias* generalises to other
moduli: for a modulus q, primes split into the φ(q) classes coprime to
q, and the counts in quadratic-non-residue classes (QNR) tend to lead
those in quadratic-residue classes (QR), although Littlewood (1914)
proved the lead changes infinitely often.

For modulus q = 9:
- Classes coprime to 9: {1, 2, 4, 5, 7, 8}
- Quadratic residues (QR): {1, 4, 7}
- Quadratic non-residues (QNR): {2, 5, 8}

This script computes

    Δ(N) = π(N; QNR) − π(N; QR)

over all primes 5 ≤ p ≤ N and analyses its behaviour.

## Results at N = 10⁷

| Quantity | Value |
|---|---|
| Total primes (excluding 2, 3) | 664 577 |
| π(QR) | 332 194 |
| π(QNR) | 332 383 |
| Δ(10⁷) | +189 |
| Maximum of Δ | +363 (at p = 7 291 259) |
| Lead changes | **0** |
| Fraction of time QNR leads | 99.999 % |
| Scaling exponent α in \|Δ\| ~ k^α | 0.4613 |

### Per-class counts

| Class | Type | Count | Deviation from mean |
|---|---|---|---|
| 1 | QR | 110 772 | +9.17 |
| 2 | QNR | 110 835 | +72.17 |
| 4 | QR | 110 743 | −19.83 |
| 5 | QNR | 110 760 | −2.83 |
| 7 | QR | 110 679 | −83.83 |
| 8 | QNR | 110 788 | +25.17 |

Mean per class: 110 762.83.

## Observations

1. **No lead change up to N = 10⁷.** QNR classes lead from prime #3
   onwards; Δ never returns to negative values. For comparison, the
   first lead change in the classical Chebyshev race mod 4 occurs at
   p = 26 861. The first lead change for mod 9 lies further than 10⁷.

2. **Scaling exponent α ≈ 0.46.** Consistent with the Hardy–Littlewood
   heuristic Δ(N) ~ √(N/log N), whose effective log-log slope at
   N = 10⁷ is ≈ 0.46. Pure random walk would give 0.5, pure linear
   drift would give 1.0. The observed value sits between, reflecting
   stochastic fluctuations around a sub-√N drift.

3. **Substructure within QR/QNR.** The bias is not uniform across
   the three QR classes (and the three QNR classes). Class 7 is
   the largest underperformer (−83.8), class 2 the strongest leader
   (+72.2). Compare Granville & Martin (2006) §5 on subrace structure.

## Usage

```bash
pip install numpy sympy matplotlib
python3 chebyshev_race_mod9.py
```

Generates console output and saves `chebyshev_race_mod9.png`.

Runtime: ~1 minute on a modern laptop.

## Reference

Granville, A. & Martin, G. (2006).
[*Prime Number Races.*](https://arxiv.org/abs/math/0408319)
American Mathematical Monthly, 113(1), 1–33.

## License

MIT
