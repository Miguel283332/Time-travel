"""
Del espectro primordial con rebote (LQC) al fondo cósmico de microondas y
comparación con los datos de Planck 2018.

Fase B, pasos 2-4 de docs/gravedad-cuantica-bucles.md.

1. Lee F(k/k_LQC) = P_LQC / P_BD de resultados/espectro_lqc.csv
   (lo produce simulaciones/lqc_perturbaciones.py).
2. Espectro primordial físico:
       P(k) = A_s (k/k_p)^(n_s-1) · F(k / k_LQC_fis)
   El único parámetro nuevo es k_LQC_fis: la escala del rebote vista hoy
   (depende del número de e-folds entre el rebote y hoy, que la teoría no fija).
   Los demás parámetros se fijan al mejor ajuste ΛCDM de Planck 2018.
3. CAMB calcula los espectros angulares TT, EE y TE.
4. Verosimilitud TT de bajo ℓ (2 ≤ ℓ ≤ 29) de Planck 2018 (método de
   Gibbs/Commander "gaussianizado"; reimplementación de la versión en Python de
   E. Gjerløw distribuida con cobaya). Se calcula Δχ² = χ²(LQC) - χ²(ΛCDM).
5. Predicción para la polarización EE a ℓ < 30: cuánto se separa de ΛCDM y
   con qué significancia la vería un experimento limitado por la varianza
   cósmica (como LiteBIRD).

Datos: se descargan la primera vez en datos/planck_2018_lowT_native/
(publicados por la colaboración Planck; empaquetados por CobayaSampler).

Uso:  python3 simulaciones/cmb_lqc.py [--espectro resultados/espectro_lqc.csv]
Requiere: numpy, scipy, camb, matplotlib
"""

import argparse
import hashlib
import io
import os
import urllib.request
import zipfile

import numpy as np
from scipy.interpolate import InterpolatedUnivariateSpline

import camb

URL_DATOS = ("https://github.com/CobayaSampler/planck_native_data/releases/download/"
             "v1/planck_2018_lowT.zip")
DIR_DATOS = os.path.join("datos", "planck_2018_lowT_native")
DIR_RES = "resultados"

# Mejor ajuste Planck 2018 (TT,TE,EE+lowE+lensing)
COSMO = dict(H0=67.36, ombh2=0.02237, omch2=0.1200, mnu=0.06, tau=0.0544)
AS, NS, KPIVOT = 2.100e-9, 0.9649, 0.05
LMIN, LMAX = 2, 29


# ---------------------------------------------------------------- datos ----

def descargar_datos():
    if os.path.exists(os.path.join(DIR_DATOS, "cov.txt")):
        return
    print(f"Descargando datos de Planck 2018 low-ℓ TT desde\n  {URL_DATOS}")
    with urllib.request.urlopen(URL_DATOS, timeout=120) as r:
        contenido = r.read()
    print(f"  {len(contenido) / 1e6:.1f} MB, sha256 {hashlib.sha256(contenido).hexdigest()[:16]}…")
    os.makedirs(DIR_DATOS, exist_ok=True)
    with zipfile.ZipFile(io.BytesIO(contenido)) as z:
        for nombre in ("cov.txt", "mu.txt", "mu_sigma.txt", "cl2x_1.txt", "cl2x_2.txt"):
            with z.open(nombre) as fuente, open(os.path.join(DIR_DATOS, nombre), "wb") as dest:
                dest.write(fuente.read())


class PlanckLowTT:
    """Verosimilitud TT de bajo ℓ de Planck 2018 (aproximación gaussianizada).

    Cada D_ℓ se transforma en una variable x_ℓ mediante un spline calibrado
    con las muestras de Gibbs; las x_ℓ siguen una gaussiana multivariante.
    """

    def __init__(self, lmin=LMIN, lmax=LMAX):
        descargar_datos()
        sl = slice(lmin - 2, lmax - 1)
        self.lmin, self.lmax = lmin, lmax
        cov = np.loadtxt(os.path.join(DIR_DATOS, "cov.txt"))[sl, sl]
        self.cov = cov
        self.covinv = np.linalg.inv(cov)
        self.mu = np.loadtxt(os.path.join(DIR_DATOS, "mu.txt"))[sl]
        mu_sigma = np.loadtxt(os.path.join(DIR_DATOS, "mu_sigma.txt"))[sl]
        cl_grid = np.loadtxt(os.path.join(DIR_DATOS, "cl2x_1.txt"))[:, sl]
        x_grid = np.loadtxt(os.path.join(DIR_DATOS, "cl2x_2.txt"))[:, sl]
        self.splines, self.dsplines, self.limites = [], [], []
        self.cl_grid, self.x_grid = cl_grid, x_grid
        for i in range(cl_grid.shape[1]):
            j = 0
            while abs(x_grid[j, i] + 5) < 1e-4:
                j += 1
            lo = cl_grid[j + 2, i]
            j = cl_grid.shape[0] - 1
            while abs(x_grid[j, i] - 5) < 1e-4:
                j -= 1
            hi = cl_grid[j - 2, i]
            self.limites.append((lo, hi))
            s = InterpolatedUnivariateSpline(cl_grid[:, i], x_grid[:, i])
            self.splines.append(s)
            self.dsplines.append(s.derivative())
        self.offset = 0.0
        self.offset = self.loglike(mu_sigma)

    def loglike(self, dl):
        """dl: D_ℓ = ℓ(ℓ+1)C_ℓ/2π en μK² para ℓ = lmin..lmax."""
        logl = 0.0
        x = np.empty_like(dl)
        for i, d in enumerate(dl):
            lo, hi = self.limites[i]
            if not lo <= d <= hi:
                return -np.inf
            der = self.dsplines[i](d)
            if der < 0:
                return -np.inf
            logl += np.log(der)
            x[i] = self.splines[i](d)
        delta = x - self.mu
        return logl - 0.5 * delta @ self.covinv @ delta - self.offset

    def puntos(self):
        """D_ℓ 'medidos' y barras de error al 68 % para dibujar.

        Para cada ℓ por separado (ignorando las correlaciones entre ℓ, solo
        para la figura): máximo de L(D) ∝ x'(D) exp(-(x(D)-μ)²/2σ²) y el
        intervalo donde ln L cae 0.5 respecto al máximo.
        """
        centro, abajo, arriba = [], [], []
        sig = np.sqrt(np.diag(self.cov))
        for i in range(len(self.mu)):
            lo, hi = self.limites[i]
            d = np.linspace(lo, hi, 20000)
            der = np.clip(self.dsplines[i](d), 1e-300, None)
            lnl = np.log(der) - 0.5 * ((self.splines[i](d) - self.mu[i]) / sig[i]) ** 2
            j = int(np.argmax(lnl))
            dentro = d[lnl >= lnl[j] - 0.5]
            centro.append(d[j])
            abajo.append(dentro.min())
            arriba.append(dentro.max())
        return np.array(centro), np.array(abajo), np.array(arriba)


# ---------------------------------------------------------------- teoría ---

def cargar_F(ruta):
    d = np.loadtxt(ruta, delimiter=",")
    return d[:, 0], d[:, -1]


def espectro_fisico(k, x_tab, F_tab, k_lqc):
    """A_s (k/k_p)^(n_s-1) F(k/k_lqc). Fuera de la tabla: F = 1 (UV) y el
    valor del extremo inferior (IR, ver advertencias)."""
    base = AS * (k / KPIVOT) ** (NS - 1)
    if k_lqc is None:
        return base
    x = k / k_lqc
    F = np.interp(np.log(x), np.log(x_tab), F_tab, left=F_tab[0], right=1.0)
    return base * F


def espectros_cmb(x_tab, F_tab, k_lqc, lmax=200):
    p = camb.set_params(**COSMO, As=AS, ns=NS, lmax=lmax)
    k = np.geomspace(1e-7, 5.0, 6000)
    p.set_initial_power_table(k, espectro_fisico(k, x_tab, F_tab, k_lqc),
                              effective_ns_for_nonlinear=NS)
    p.WantTensors = False
    res = camb.get_results(p)
    cl = res.get_cmb_power_spectra(p, CMB_unit="muK", lmax=lmax)["total"]
    return cl, res.get_derived_params()["DAstar"] * 1000.0  # Mpc


def significancia_ee(cl_lqc, cl_ref, lmin=2, lmax=29, fsky=0.7):
    """Δχ² de Fisher para distinguir EE(LQC) de EE(ΛCDM) con un experimento
    limitado por la varianza cósmica en una fracción de cielo fsky."""
    l = np.arange(lmin, lmax + 1)
    a, b = cl_lqc[l, 1], cl_ref[l, 1]
    return float(np.sum(fsky * (2 * l + 1) / 2 * ((a - b) / b) ** 2))


# ---------------------------------------------------------------- main -----

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--espectro", default=os.path.join(DIR_RES, "espectro_lqc.csv"))
    ap.add_argument("--etiqueta", default="", help="sufijo para los archivos de salida")
    ap.add_argument("--n", type=int, default=36, help="puntos del barrido en k_LQC")
    args = ap.parse_args()
    suf = f"_{args.etiqueta}" if args.etiqueta else ""

    x_tab, F_tab = cargar_F(args.espectro)
    like = PlanckLowTT()
    l = np.arange(LMIN, LMAX + 1)

    cl_ref, DA = espectros_cmb(x_tab, F_tab, None)
    lnL_ref = like.loglike(cl_ref[l, 0])
    print(f"Distancia comóvil a la última dispersión D_A = {DA:.0f} Mpc")
    print(f"ΛCDM (Planck 2018): -2 lnL_lowTT = {-2 * lnL_ref:.2f}")

    kls = np.geomspace(1e-5, 3e-3, args.n)
    filas = []
    for kl in kls:
        cl, _ = espectros_cmb(x_tab, F_tab, kl)
        lnL = like.loglike(cl[l, 0])
        dchi2 = -2 * (lnL - lnL_ref)
        sig_ee = significancia_ee(cl, cl_ref)
        filas.append((kl, kl * DA, dchi2, sig_ee))
        print(f"  k_LQC = {kl:.2e} Mpc⁻¹ (ℓ ≈ {kl * DA:6.1f})   Δχ²_TT = {dchi2:+8.2f}   "
              f"Δχ²_EE(futuro) = {sig_ee:7.2f}")
    filas = np.array(filas)
    os.makedirs(DIR_RES, exist_ok=True)
    np.savetxt(os.path.join(DIR_RES, f"barrido_kLQC{suf}.csv"), filas, delimiter=",",
               header="k_LQC_Mpc^-1, ell_LQC, Delta_chi2_lowTT_Planck, Delta_chi2_EE_varianza_cosmica",
               fmt="%.6e")

    # Resumen
    i_best = int(np.argmin(filas[:, 2]))
    excl = filas[filas[:, 2] > 4]
    permitido = filas[filas[:, 2] <= 4]
    print("\nResumen")
    print(f"  Mejor ajuste: k_LQC = {filas[i_best, 0]:.2e} Mpc⁻¹ (ℓ ≈ {filas[i_best, 1]:.1f}), "
          f"Δχ² = {filas[i_best, 2]:+.2f} respecto a ΛCDM")
    if len(excl):
        print(f"  Excluido (Δχ² > 4): {len(excl)} de {len(filas)} valores; "
              f"el menor excluido es k_LQC = {excl[:, 0].min():.2e} Mpc⁻¹ (ℓ ≈ {excl[:, 1].min():.1f})")
    if len(permitido):
        kmax_ok = permitido[:, 0].max()
        print(f"  Mayor k_LQC permitido (Δχ² ≤ 4): {kmax_ok:.2e} Mpc⁻¹")

    graficas(like, x_tab, F_tab, filas, cl_ref, suf)


def graficas(like, x_tab, F_tab, filas, cl_ref, suf):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    l = np.arange(LMIN, LMAX + 1)
    fig, ax = plt.subplots(2, 2, figsize=(12, 9))

    a = ax[0, 0]
    a.semilogx(x_tab, F_tab, color="#3060a0", lw=1.2)
    a.axhline(1, color="0.5", lw=0.8, ls="--")
    a.set_xlabel(r"$k / k_{\rm LQC}$")
    a.set_ylabel(r"$F = P_{\rm LQC} / P_{\rm BD}$")
    a.set_yscale("log")
    a.set_title("Efecto del rebote en el espectro primordial")

    a = ax[0, 1]
    c, lo, hi = like.puntos()
    a.errorbar(l, c, yerr=[c - lo, hi - c], fmt="o", color="k", ms=4, label="Planck 2018 (ℓ<30)")
    a.plot(np.arange(2, 60), cl_ref[2:60, 0], color="0.2", lw=1.5, ls="--", label="ΛCDM", zorder=5)
    i_best = int(np.argmin(filas[:, 2]))
    elegidos = sorted({i_best, int(np.argmin(abs(filas[:, 2] - 4))), len(filas) - 1})
    colores = ["#c03030", "#e08020", "#3060a0"]
    for col, i in zip(colores, elegidos):
        cl, _ = espectros_cmb(x_tab, F_tab, filas[i, 0])
        a.plot(np.arange(2, 60), cl[2:60, 0], color=col, lw=1.2,
               label=fr"LQC $\ell_{{\rm LQC}}\approx{filas[i, 1]:.0f}$ ($\Delta\chi^2$={filas[i, 2]:+.1f})")
    a.set_xscale("log")
    a.set_xlabel(r"$\ell$")
    a.set_ylabel(r"$D_\ell^{TT}\ [\mu K^2]$")
    a.set_ylim(0, 3500)
    a.legend(fontsize=8)
    a.set_title("Temperatura a gran escala")

    a = ax[1, 0]
    a.semilogx(filas[:, 1], filas[:, 2], "o-", color="#3060a0", ms=3)
    a.axhline(0, color="0.5", lw=0.8)
    a.axhline(4, color="#c03030", lw=0.8, ls="--", label=r"$\Delta\chi^2 = 4$")
    a.set_xlabel(r"$\ell_{\rm LQC} \approx k_{\rm LQC} D_A$")
    a.set_ylabel(r"$\Delta\chi^2$ (Planck low-ℓ TT)")
    a.legend()
    a.set_title("Ajuste a Planck en función de la escala del rebote")

    a = ax[1, 1]
    for col, i in zip(colores, elegidos):
        cl, _ = espectros_cmb(x_tab, F_tab, filas[i, 0])
        a.plot(np.arange(2, 60), cl[2:60, 1] / cl_ref[2:60, 1], color=col, lw=1.2,
               label=fr"$\ell_{{\rm LQC}}\approx{filas[i, 1]:.0f}$")
    ll = np.arange(2, 60)
    cv = np.sqrt(2 / (2 * ll + 1) / 0.7)
    a.fill_between(ll, 1 - cv, 1 + cv, color="0.85", label="varianza cósmica (fsky=0.7)")
    a.axhline(1, color="0.4", lw=0.8)
    a.set_xscale("log")
    a.set_xlabel(r"$\ell$")
    a.set_ylabel(r"$C_\ell^{EE}({\rm LQC}) / C_\ell^{EE}(\Lambda{\rm CDM})$")
    a.legend(fontsize=8)
    a.set_title("Predicción: polarización E a gran escala")

    from matplotlib.ticker import FixedLocator, NullFormatter, ScalarFormatter
    for a in (ax[0, 1], ax[1, 1]):
        a.xaxis.set_major_locator(FixedLocator([2, 3, 5, 10, 20, 30, 50]))
        a.xaxis.set_major_formatter(ScalarFormatter())
        a.xaxis.set_minor_formatter(NullFormatter())
    fig.tight_layout()
    ruta = os.path.join(DIR_RES, f"cmb_lqc{suf}.png")
    fig.savefig(ruta, dpi=130)
    print(f"Figura: {ruta}")


if __name__ == "__main__":
    main()
