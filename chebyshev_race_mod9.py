#!/usr/bin/env python3
"""
Chebyshev Bias modulo 9 — Numerical Study of the Prime Number Race
====================================================================

Computes the race between primes in quadratic-residue (QR) and
quadratic-non-residue (QNR) classes modulo 9, up to N = 10^7.

Quadratic residues mod 9     : {1, 4, 7}
Quadratic non-residues mod 9 : {2, 5, 8}
(All classes coprime to 9; classes {0, 3, 6} are excluded.)

Output:
  - Lead change count (sign changes of Delta(N) = pi(N; QNR) - pi(N; QR))
  - Drift exponent alpha in |Delta| ~ N^alpha
  - Per-class statistics
  - Plot: chebyshev_race_mod9.png

Theoretical reference:
  Granville, A. & Martin, G. (2006). "Prime Number Races."
  American Mathematical Monthly, 113(1), 1-33.

Author : Franz Hettich
License: MIT
"""

import numpy as np
from sympy import primerange
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# -----------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------
N_MAX = 10**7
QR    = {1, 4, 7}      # quadratic residues mod 9
QNR   = {2, 5, 8}      # quadratic non-residues mod 9
OUT_PNG = "chebyshev_race_mod9.png"


# -----------------------------------------------------------------------
# Generate primes and class labels
# -----------------------------------------------------------------------
print(f"Generating primes up to {N_MAX:,} ...")
primes = np.array(list(primerange(5, N_MAX + 1)), dtype=np.int64)
mod9   = primes % 9
N_p    = len(primes)
print(f"  -> {N_p:,} primes (excluding 2 and 3)")


# -----------------------------------------------------------------------
# Per-class statistics at final N
# -----------------------------------------------------------------------
counts = {r: int((mod9 == r).sum()) for r in [1, 2, 4, 5, 7, 8]}
mu     = np.mean(list(counts.values()))
sigma_each = np.sqrt(mu * (1 - 1/6))   # multinomial 1-bin stdev

print(f"\nPer-class counts at N = {N_MAX:,} (mean = {mu:.2f})")
print(f"{'class':>6} | {'type':>4} | {'count':>10} | {'dev':>10} | {'sigma':>8}")
print("-" * 55)
for r in [1, 2, 4, 5, 7, 8]:
    typ = "QR " if r in QR else "QNR"
    dev = counts[r] - mu
    print(f"{r:>6} | {typ:>4} | {counts[r]:>10,} | {dev:>+10.2f} | {dev/sigma_each:>+8.3f}")

n_qr  = sum(counts[r] for r in QR)
n_qnr = sum(counts[r] for r in QNR)
print(f"\nTotal QR  = {n_qr:,}")
print(f"Total QNR = {n_qnr:,}")
print(f"Delta     = QNR - QR = {n_qnr - n_qr:+d}")


# -----------------------------------------------------------------------
# Cumulative race curve
# -----------------------------------------------------------------------
delta_step = np.where(np.isin(mod9, list(QNR)), 1, -1)
delta_cum  = np.cumsum(delta_step)
k_arr      = np.arange(1, N_p + 1)


# -----------------------------------------------------------------------
# Lead changes (sign reversals of delta_cum)
# -----------------------------------------------------------------------
sign = np.sign(delta_cum)
lead_changes = []
prev = 0
for i, s in enumerate(sign):
    if s != 0 and prev != 0 and s != prev:
        lead_changes.append((i, primes[i], int(delta_cum[i])))
    if s != 0:
        prev = s

print(f"\nLead changes up to N = {N_MAX:,}: {len(lead_changes)}")
if lead_changes:
    print("First 10 lead changes:")
    for i, p, d in lead_changes[:10]:
        print(f"  prime #{i+1:,}  p = {p:,}  Delta = {d:+d}")


# -----------------------------------------------------------------------
# Curve statistics
# -----------------------------------------------------------------------
print(f"\nCurve statistics:")
print(f"  Endpoint Delta(N_max) = {delta_cum[-1]:+d}")
print(f"  Maximum               = {delta_cum.max():+d}  at prime #{delta_cum.argmax()+1:,} = {primes[delta_cum.argmax()]:,}")
print(f"  Minimum               = {delta_cum.min():+d}  at prime #{delta_cum.argmin()+1:,} = {primes[delta_cum.argmin()]:,}")
print(f"  Fraction QNR leading  = {(delta_cum > 0).mean()*100:.3f} %")
print(f"  Fraction QR  leading  = {(delta_cum < 0).mean()*100:.3f} %")
print(f"  Fraction tied         = {(delta_cum == 0).mean()*100:.3f} %")


# -----------------------------------------------------------------------
# Drift / scaling analysis
# -----------------------------------------------------------------------
# Linear fit Delta(k) = a + b*k
b, a = np.polyfit(k_arr, delta_cum, 1)
print(f"\nLinear fit:  Delta(k) ~ {a:+.2f} + {b:.4e} * k")

# Power-law fit for |Delta| ~ k^alpha
mask = delta_cum > 10
if mask.sum() > 100:
    alpha, _ = np.polyfit(np.log(k_arr[mask]), np.log(delta_cum[mask]), 1)
    print(f"Scaling law: |Delta| ~ k^alpha,  alpha = {alpha:.4f}")
    print(f"  (alpha = 0.5 -> pure random walk, alpha = 1.0 -> linear drift)")
    print(f"  Hardy-Littlewood prediction:  Delta ~ sqrt(N/log N), "
          f"effective alpha at N=1e7: ~0.46")


# -----------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------
fig, axes = plt.subplots(2, 2, figsize=(13, 9))

# Full curve
axes[0, 0].plot(k_arr, delta_cum, lw=0.5)
axes[0, 0].axhline(0, color="k", lw=0.5)
axes[0, 0].set_xlabel("Prime index k")
axes[0, 0].set_ylabel(r"$\Delta(k) = \pi(\mathrm{QNR}) - \pi(\mathrm{QR})$")
axes[0, 0].set_title(f"Chebyshev race mod 9 up to $N = 10^7$")
axes[0, 0].grid(alpha=0.3)

# First 1000 primes
axes[0, 1].plot(k_arr[:1000], delta_cum[:1000], lw=0.7)
axes[0, 1].axhline(0, color="k", lw=0.5)
axes[0, 1].set_xlabel("Prime index k")
axes[0, 1].set_ylabel(r"$\Delta(k)$")
axes[0, 1].set_title("First 1000 primes")
axes[0, 1].grid(alpha=0.3)

# Delta vs sqrt(k)
axes[1, 0].plot(np.sqrt(k_arr), delta_cum, lw=0.5)
axes[1, 0].set_xlabel(r"$\sqrt{k}$")
axes[1, 0].set_ylabel(r"$\Delta(k)$")
axes[1, 0].set_title(r"$\Delta$ vs $\sqrt{k}$ (random walk reference)")
axes[1, 0].grid(alpha=0.3)

# Log-log
mask2 = delta_cum > 0
axes[1, 1].loglog(k_arr[mask2], delta_cum[mask2], lw=0.5, label=r"$\Delta(k)$")
axes[1, 1].loglog(k_arr[mask2], np.sqrt(k_arr[mask2]), "r--", lw=1,
                  label=r"$\sqrt{k}$ (random walk)")
axes[1, 1].loglog(k_arr[mask2], b * k_arr[mask2], "g--", lw=1,
                  label=f"linear drift {b:.2e}·k")
axes[1, 1].set_xlabel("k")
axes[1, 1].set_ylabel(r"$\Delta$")
axes[1, 1].set_title("Log-log scaling")
axes[1, 1].legend()
axes[1, 1].grid(alpha=0.3, which="both")

plt.tight_layout()
plt.savefig(OUT_PNG, dpi=120)
print(f"\nPlot saved to {OUT_PNG}")
