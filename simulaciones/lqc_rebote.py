"""
Cosmología cuántica de bucles (LQC): el rebote que sustituye al Big Bang.

Universo plano y homogéneo lleno de un campo escalar sin masa (ρ = p).
Unidades de Planck: G = ħ = c = 1.

  Relatividad general:  H² = (8π/3) ρ
                        Ḣ  = -8π ρ                    -> singularidad (a -> 0)
  LQC efectiva:         H² = (8π/3) ρ (1 - ρ/ρ_c)
                        Ḣ  = -8π ρ (1 - 2ρ/ρ_c)        -> rebote en ρ = ρ_c

con ρ = ρ_c (a_B / a)^6 y  ρ_c = sqrt(3) / (32 π² γ³) ≈ 0.41 ρ_Planck.

Solución exacta de la versión efectiva:  a(t)^6 = a_B^6 (1 + 24 π ρ_c t²).
El script integra numéricamente las ecuaciones y la compara con la exacta.

Uso:  python3 simulaciones/lqc_rebote.py
Solo usa la biblioteca estándar.
"""

import math

GAMMA = 0.2375329
RHO_C = math.sqrt(3) / (32 * math.pi**2 * GAMMA**3)
A_B = 1.0  # factor de escala en el rebote (normalización libre)


def rho(lna):
    return RHO_C * math.exp(-6 * lna) * A_B**6


def derivadas(estado, cuantica):
    lna, H = estado
    r = rho(lna)
    if cuantica:
        dH = -8 * math.pi * r * (1 - 2 * r / RHO_C)
    else:
        dH = -8 * math.pi * r
    return (H, dH)


def rk4(estado, dt, cuantica):
    def suma(e, d, f):
        return (e[0] + f * d[0], e[1] + f * d[1])

    k1 = derivadas(estado, cuantica)
    k2 = derivadas(suma(estado, k1, dt / 2), cuantica)
    k3 = derivadas(suma(estado, k2, dt / 2), cuantica)
    k4 = derivadas(suma(estado, k3, dt), cuantica)
    return (
        estado[0] + dt / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]),
        estado[1] + dt / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]),
    )


def exacta(t):
    x = A_B**6 * (1 + 24 * math.pi * RHO_C * t**2)
    a = x ** (1 / 6)
    H = 8 * math.pi * RHO_C * t / (1 + 24 * math.pi * RHO_C * t**2)
    return math.log(a), H


def main():
    print(f"ρ_c = {RHO_C:.4f} ρ_Planck  (γ = {GAMMA})")
    t0, t1, dt = -3.0, 3.0, 1e-4
    lna0, H0 = exacta(t0)

    # LQC efectiva: atraviesa t = 0 sin problemas
    estado, t = (lna0, H0), t0
    rho_max, a_min, err_max = 0.0, float("inf"), 0.0
    filas = []
    marcas = {round(x, 1) for x in (-3, -1, -0.3, -0.1, 0, 0.1, 0.3, 1, 3)}
    while t < t1 - dt / 2:
        estado = rk4(estado, dt, cuantica=True)
        t += dt
        rho_max = max(rho_max, rho(estado[0]))
        a_min = min(a_min, math.exp(estado[0]))
        err_max = max(err_max, abs(estado[0] - exacta(t)[0]))
        if abs(t - round(t, 1)) < dt / 2 and round(t, 1) in marcas:
            filas.append((t, math.exp(estado[0]), rho(estado[0]) / RHO_C, estado[1] + 0.0))

    print("\nLQC efectiva (contracción -> rebote -> expansión)")
    print("     t        a(t)      ρ/ρ_c        H")
    for t, a, r, H in filas:
        print(f"  {t:+6.2f}   {a:8.4f}   {r:8.4f}   {(0.0 if abs(H) < 5e-5 else H):+9.4f}")
    print(f"  ρ máxima alcanzada = {rho_max / RHO_C:.6f} ρ_c   (a mínimo = {a_min:.6f})")
    print(f"  error máximo frente a la solución exacta en ln a: {err_max:.1e}")

    # Relatividad general: mismo estado inicial, colapsa
    estado, t = (lna0, -math.sqrt(8 * math.pi / 3 * rho(lna0))), t0
    while t < t1:
        estado = rk4(estado, dt, cuantica=False)
        t += dt
        if not math.isfinite(estado[1]) or rho(estado[0]) > 1e6 * RHO_C:
            break
    print("\nRelatividad general con el mismo estado inicial")
    print(f"  ρ supera 10⁶ ρ_c en t ≈ {t:+.3f}: singularidad del Big Bang/Big Crunch")
    print("\n-> En LQC la geometría cuántica actúa como repulsión cuando ρ ~ ρ_c.")


if __name__ == "__main__":
    main()
