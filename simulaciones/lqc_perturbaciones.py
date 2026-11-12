"""
Espectro primordial de perturbaciones escalares con rebote cuántico (LQC).

Fase B, paso 1 de docs/gravedad-cuantica-bucles.md.

Modelo
------
Unidades de Planck (G = ħ = c = 1, m_Pl = 1). Inflatón con potencial de
Starobinsky
    V(φ) = 3M²/(32π) (1 - exp(-sqrt(16π/3) φ))²,  M = 2.56e-6
(M fijado para que A_s ≈ 2.1e-9 con N ≈ 55).

Fondo: dinámica efectiva de la cosmología cuántica de bucles
    H²  = (8π/3) ρ (1 - ρ/ρ_c),     ρ_c = sqrt(3)/(32π² γ³) ≈ 0.41
    Ḣ   = -4π φ̇² (1 - 2ρ/ρ_c)
    φ̈  + 3Hφ̇ + V'(φ) = 0
con el rebote en t = 0 (a_B = 1, ρ = ρ_c).

Perturbaciones: enfoque de "métrica vestida" (Agullo, Ashtekar y Nelson).
La perturbación del inflatón en el gauge espacialmente plano obedece
    Q̈ + 3H Q̇ + (k²/a² + 𝔘) Q = 0,
    𝔘 = V_φφ + f V + 2 s sqrt(f) V_φ,   f = 24π φ̇²/ρ,   s = signo(φ̇ H).
En relatividad general 𝔘 coincide exactamente con el potencial de
Mukhanov-Sasaki; escrito con ρ en lugar de H² queda finito en el rebote.
La perturbación de curvatura es R = H Q / φ̇ y P_R = k³|R|²/(2π²).

Estado inicial: vacío adiabático (orden 0) en la fase de contracción, en
t_i < 0, cuando la energía cinética domina. Es una elección: el espectro a
k < k_LQC depende del estado inicial (ver la sección de advertencias del
documento). El script permite cambiar t_i para medir esa dependencia.

Referencia: el mismo modo empezando en vacío de Bunch-Davies al inicio de la
inflación (sin rebote). F(k) = P_LQC / P_BD aísla el efecto del rebote.

Para k < k_LQC, F oscila en k con un periodo Δk/k ~ 1e-4, muchísimo más
fino de lo que el fondo cósmico puede resolver. Lo observable es su promedio,
que se obtiene de forma exacta con el coeficiente de Bogoliubov β_k medido
en la fase cinética tras el rebote:  <F> = 1 + 2|β_k|²  (modo por defecto).
La evolución completa de cada modo está en --exacto (lento; sirve de
validación).

Uso:  python3 simulaciones/lqc_perturbaciones.py [--ti=-1e4] [--exacto]
Salida: resultados/espectro_lqc.csv (promedio) o resultados/espectro_lqc_exacto.csv
Requiere: numpy, numba
"""

import argparse
import math
import os
import time

import numpy as np
from numba import njit, prange

GAMMA = 0.2375329
RHO_C = math.sqrt(3) / (32 * math.pi**2 * GAMMA**3)
K_LQC = math.sqrt(8 * math.pi * RHO_C)  # sqrt(a''/a) en el rebote, a_B = 1
M_STAR = 2.56e-6
LAM = math.sqrt(16 * math.pi / 3)
V0 = 3 * M_STAR**2 / (32 * math.pi)
PHI_B = 3.54  # valor del inflatón en el rebote (≈ 70 e-folds de inflación)

PASOS_POR_OSCILACION = 200.0


@njit(cache=True)
def V(p):
    e = math.exp(-LAM * p)
    return V0 * (1 - e) ** 2


@njit(cache=True)
def Vp(p):
    e = math.exp(-LAM * p)
    return 2 * V0 * LAM * e * (1 - e)


@njit(cache=True)
def Vpp(p):
    e = math.exp(-LAM * p)
    return 2 * V0 * LAM * LAM * e * (2 * e - 1)


@njit(cache=True)
def derivadas(y, k):
    # y = [ln a, H, φ, φ̇, Re Q, Im Q, Re Q̇, Im Q̇]
    lna, H, p, pd = y[0], y[1], y[2], y[3]
    rho = 0.5 * pd * pd + V(p)
    f = 24 * math.pi * pd * pd / rho
    s = 1.0 if pd * H >= 0 else -1.0
    U = Vpp(p) + f * V(p) + 2 * s * math.sqrt(f) * Vp(p)
    w2 = k * k * math.exp(-2 * lna) + U
    d = np.empty(8)
    d[0] = H
    d[1] = -4 * math.pi * pd * pd * (1 - 2 * rho / RHO_C)
    d[2] = pd
    d[3] = -3 * H * pd - Vp(p)
    d[4] = y[6]
    d[5] = y[7]
    d[6] = -3 * H * y[6] - w2 * y[4]
    d[7] = -3 * H * y[7] - w2 * y[5]
    return d


@njit(cache=True)
def proyectar(y):
    """Impone la ligadura de Friedmann efectiva lejos del rebote.

    Integrar solo Ḣ acumula un error que se comporta como curvatura
    espacial y crece; se corrige fijando |H| con la ligadura.
    """
    rho = 0.5 * y[3] ** 2 + V(y[2])
    if rho < 0.5 * RHO_C:
        h = math.sqrt(8 * math.pi / 3 * rho * (1 - rho / RHO_C))
        y[1] = h if y[1] > 0 else -h
    return y


@njit(cache=True)
def paso_rk4(y, dt, k):
    k1 = derivadas(y, k)
    k2 = derivadas(y + 0.5 * dt * k1, k)
    k3 = derivadas(y + 0.5 * dt * k2, k)
    k4 = derivadas(y + dt * k3, k)
    return proyectar(y + dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4))


@njit(cache=True)
def paso_temporal(y, t, k):
    H = abs(y[1])
    dt = min(0.002 + 0.001 * abs(t), 0.01 / max(H, 1e-30))
    w2 = k * k * math.exp(-2 * y[0])
    if w2 > 0:
        dt = min(dt, 2 * math.pi / math.sqrt(w2) / PASOS_POR_OSCILACION)
    return dt


@njit(cache=True)
def fondo_hasta(t_obj):
    """Integra solo el fondo desde el rebote (t = 0) hasta t_obj."""
    y = np.zeros(8)
    y[2] = PHI_B
    y[3] = -math.sqrt(2 * (RHO_C - V(PHI_B)))
    t = 0.0
    signo = 1.0 if t_obj > 0 else -1.0
    while signo * (t_obj - t) > 0:
        dt = min(paso_temporal(y, t, 0.0), signo * (t_obj - t))
        y = paso_rk4(y, signo * dt, 0.0)
        t += signo * dt
    return y


@njit(cache=True)
def inicio_inflacion():
    """Tiempo y estado cuando ε = 1.5 φ̇²/ρ cae por debajo de 1."""
    y = np.zeros(8)
    y[2] = PHI_B
    y[3] = -math.sqrt(2 * (RHO_C - V(PHI_B)))
    t = 0.0
    while True:
        dt = paso_temporal(y, t, 0.0)
        y = paso_rk4(y, dt, 0.0)
        t += dt
        rho = 0.5 * y[3] ** 2 + V(y[2])
        if t > 1 and 1.5 * y[3] ** 2 / rho < 1:
            return t, y


@njit(cache=True)
def vacio_adiabatico(y, k):
    """Condición inicial de vacío adiabático de orden 0 en la variable v = aQ."""
    a = math.exp(y[0])
    H, p, pd = y[1], y[2], y[3]
    rho = 0.5 * pd * pd + V(p)
    f = 24 * math.pi * pd * pd / rho
    s = 1.0 if pd * H >= 0 else -1.0
    U = Vpp(p) + f * V(p) + 2 * s * math.sqrt(f) * Vp(p)
    Hd = -4 * math.pi * pd * pd * (1 - 2 * rho / RHO_C)
    app_a = a * a * (Hd + 2 * H * H)  # a''/a en tiempo conforme
    w = math.sqrt(k * k + a * a * U - app_a)
    # v = 1/sqrt(2w), v' = -i w v ;  Q = v/a ;  Q̇ = (v' - aH v)/a²
    v = 1 / math.sqrt(2 * w)
    y[4] = v / a
    y[5] = 0.0
    y[6] = (-a * H * v) / (a * a)
    y[7] = (-w * v) / (a * a)
    return y


@njit(cache=True)
def evolucionar_modo(y0, t0, k):
    """Evoluciona un modo hasta que está congelado fuera del horizonte."""
    y = vacio_adiabatico(y0.copy(), k)
    t = t0
    while True:
        dt = paso_temporal(y, t, k)
        y = paso_rk4(y, dt, k)
        t += dt
        if t > 1:
            a = math.exp(y[0])
            if k / (a * abs(y[1])) < 1e-3:
                R2 = (y[4] ** 2 + y[5] ** 2) * (y[1] / y[3]) ** 2
                return k**3 * R2 / (2 * math.pi**2)


@njit(cache=True)
def bogoliubov_modo(y0, t0, k, t_m):
    """Evoluciona el modo hasta t_m (fase cinética tras el rebote, donde el
    modo vuelve a ser adiabático) y lo proyecta sobre las soluciones
    adiabáticas de frecuencia positiva y negativa:
        v = α u + β u*,   u = e^{-iωη}/sqrt(2ω)
        |α|² = |ωv + i v'|²/(2ω),   |β|² = |ωv - i v'|²/(2ω)
    Devuelve (|β|², |α|² - |β|² - 1); lo segundo debe ser ~0 (wronskiano).
    """
    y = vacio_adiabatico(y0.copy(), k)
    t = t0
    while t < t_m:
        dt = min(paso_temporal(y, t, k), t_m - t)
        y = paso_rk4(y, dt, k)
        t += dt
    a = math.exp(y[0])
    H, p, pd = y[1], y[2], y[3]
    rho = 0.5 * pd * pd + V(p)
    f = 24 * math.pi * pd * pd / rho
    s = 1.0 if pd * H >= 0 else -1.0
    U = Vpp(p) + f * V(p) + 2 * s * math.sqrt(f) * Vp(p)
    Hd = -4 * math.pi * pd * pd * (1 - 2 * rho / RHO_C)
    w = math.sqrt(k * k + a * a * U - a * a * (Hd + 2 * H * H))
    vr, vi = a * y[4], a * y[5]
    dvr, dvi = a * a * (H * y[4] + y[6]), a * a * (H * y[5] + y[7])
    # ωv ∓ i v'  (parte real, parte imaginaria)
    b2 = ((w * vr + dvi) ** 2 + (w * vi - dvr) ** 2) / (2 * w)
    a2 = ((w * vr - dvi) ** 2 + (w * vi + dvr) ** 2) / (2 * w)
    return b2, a2 - b2 - 1.0


@njit(parallel=True, cache=True)
def bogoliubov(ks, y0, t0, t_m):
    out = np.empty((ks.size, 2))
    for i in prange(ks.size):
        b2, w = bogoliubov_modo(y0, t0, ks[i], t_m)
        out[i, 0] = b2
        out[i, 1] = w
    return out


@njit(parallel=True, cache=True)
def espectro(ks, y0, t0):
    out = np.empty(ks.size)
    for i in prange(ks.size):
        out[i] = evolucionar_modo(y0, t0, ks[i])
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ti", type=float, default=-1e4, help="tiempo inicial (contracción); use --ti=-1e3")
    ap.add_argument("--tm", type=float, default=1e4, help="tiempo de proyección tras el rebote")
    ap.add_argument("--n", type=int, default=400, help="número de modos")
    ap.add_argument("--kmin", type=float, default=0.003, help="k mínimo en unidades de k_LQC")
    ap.add_argument("--kmax", type=float, default=30.0, help="k máximo en unidades de k_LQC")
    ap.add_argument("--exacto", action="store_true",
                    help="evolución completa de cada modo hasta salir del horizonte (lento)")
    ap.add_argument("--salida", default=None)
    args = ap.parse_args()

    t_on, y_on = inicio_inflacion()
    print(f"ρ_c = {RHO_C:.4f}   k_LQC = {K_LQC:.4f}   φ_B = {PHI_B}")
    print(f"Inicio de la inflación: t = {t_on:.3e}, ln a = {y_on[0]:.3f}, "
          f"k_I = a H = {math.exp(y_on[0]) * abs(y_on[1]):.3e}")

    y_i = fondo_hasta(args.ti)
    rho_i = 0.5 * y_i[3] ** 2 + V(y_i[2])
    print(f"Estado inicial en t_i = {args.ti:.1e}: ln a = {y_i[0]:.3f}, "
          f"ρ = {rho_i:.2e}, V/ρ = {V(y_i[2]) / rho_i:.1e}")

    xs = np.geomspace(args.kmin, args.kmax, args.n)
    ks = xs * K_LQC
    t0 = time.time()
    if args.exacto:
        salida = args.salida or os.path.join("resultados", "espectro_lqc_exacto.csv")
        P_lqc = espectro(ks, y_i, args.ti)
        P_bd = espectro(ks, y_on, t_on)
        F = P_lqc / P_bd
        columnas = np.column_stack([xs, P_lqc, P_bd, F])
        nombres = "k_sobre_kLQC, P_LQC, P_BD, F=P_LQC/P_BD"
    else:
        salida = args.salida or os.path.join("resultados", "espectro_lqc.csv")
        res = bogoliubov(ks, y_i, args.ti, args.tm)
        F = 1 + 2 * res[:, 0]
        columnas = np.column_stack([xs, res[:, 0], res[:, 1], F])
        nombres = "k_sobre_kLQC, beta2, error_wronskiano, F_promedio=1+2*beta2"
        print(f"Error máximo del wronskiano |α|²-|β|²-1: {np.abs(res[:, 1]).max():.1e}")
    print(f"{args.n} modos en {time.time() - t0:.1f} s")

    os.makedirs(os.path.dirname(salida), exist_ok=True)
    cabecera = (
        f"Espectro escalar LQC (métrica vestida), Starobinsky M={M_STAR}, phi_B={PHI_B}, "
        f"gamma={GAMMA}, t_i={args.ti}, t_m={args.tm}, k_LQC={K_LQC:.6f}\n" + nombres
    )
    np.savetxt(salida, columnas, delimiter=",", header=cabecera, fmt="%.8e")
    print(f"Guardado en {salida}")
    paso = max(1, args.n // 12)
    for x, f in zip(xs[::paso], F[::paso]):
        print(f"  k/k_LQC = {x:8.3f}   F = {f:10.4f}")


if __name__ == "__main__":
    main()
