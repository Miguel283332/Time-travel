"""
Compara el efecto del rebote para distintos estados iniciales de las
perturbaciones y su ajuste a Planck (docs/resultados-cmb.md, sección 4).

Necesita haber ejecutado antes:
  python3 simulaciones/lqc_perturbaciones.py
  python3 simulaciones/lqc_perturbaciones.py --estado no-oscilante --ta=-1e4
  python3 simulaciones/lqc_perturbaciones.py --estado no-oscilante --ta=0 \
          --salida resultados/espectro_lqc_no_oscilante_tras_rebote.csv
  python3 simulaciones/cmb_lqc.py
  python3 simulaciones/cmb_lqc.py --espectro resultados/espectro_lqc_no_oscilante.csv \
          --etiqueta no_oscilante

Uso:  python3 simulaciones/comparar_estados.py
"""

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

R = "resultados"
ESTADOS = [
    ("Vacío adiabático en la contracción (t = -10⁴)", "espectro_lqc.csv", "barrido_kLQC.csv", "#c03030"),
    ("Vacío no oscilante en [-10⁴, 10⁴]", "espectro_lqc_no_oscilante.csv",
     "barrido_kLQC_no_oscilante.csv", "#e08020"),
    ("Vacío no oscilante en [0, 10⁴] (tras el rebote)", "espectro_lqc_no_oscilante_tras_rebote.csv",
     None, "#3060a0"),
]


def cruce(b, nivel=4.0):
    """k_LQC donde Δχ² cruza 'nivel' (interpolación logarítmica)."""
    k, d = b[:, 0], b[:, 2]
    i = np.argmax(d > nivel)
    if i == 0:
        return None
    f = (nivel - d[i - 1]) / (d[i] - d[i - 1])
    return k[i - 1] * (k[i] / k[i - 1]) ** f, b[i - 1, 1] * (b[i, 1] / b[i - 1, 1]) ** f


def main():
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
    for nombre, esp, barr, col in ESTADOS:
        d = np.loadtxt(os.path.join(R, esp), delimiter=",")
        ax[0].loglog(d[:, 0], d[:, -1], color=col, lw=1.5, label=nombre)
        if barr:
            b = np.loadtxt(os.path.join(R, barr), delimiter=",")
            ax[1].semilogx(b[:, 1], b[:, 2], "o-", ms=3, color=col, label=nombre)
            c = cruce(b)
            if c:
                print(f"{nombre}: k_LQC < {c[0]:.2e} Mpc⁻¹ (ℓ_LQC < {c[1]:.1f}, "
                      f"{c[0] / 2.247e-4:.1f} a₀H₀)")
        else:
            print(f"{nombre}: F - 1 < {np.abs(d[:, -1] - 1).max():.3f} en todo el rango; sin cota")
    ax[0].axhline(1, color="0.5", lw=0.8, ls="--")
    ax[0].set_xlabel(r"$k/k_{\rm LQC}$")
    ax[0].set_ylabel(r"$\langle F\rangle = P_{\rm LQC}/P_{\rm BD}$")
    ax[0].set_title("Efecto del rebote según el estado inicial")
    ax[0].legend(fontsize=8)
    ax[1].axhline(4, color="0.3", lw=0.8, ls="--")
    ax[1].axhline(0, color="0.6", lw=0.8)
    ax[1].set_ylim(-2, 40)
    ax[1].set_xlabel(r"$\ell_{\rm LQC}\approx k_{\rm LQC}D_A$")
    ax[1].set_ylabel(r"$\Delta\chi^2$ respecto a ΛCDM (Planck low-ℓ TT)")
    ax[1].set_title("Ajuste a Planck")
    ax[1].legend(fontsize=8)
    fig.tight_layout()
    ruta = os.path.join(R, "comparar_estados.png")
    fig.savefig(ruta, dpi=130)
    print(f"Figura: {ruta}")


if __name__ == "__main__":
    main()
