"""
Órdenes de magnitud que sostienen el planteamiento (docs/planteamiento.md).

  1. Viaje al futuro ya medido: dilatación temporal en los relojes GPS y en
     relojes ópticos separados por 33 cm de altura.
  2. Energía negativa disponible en el laboratorio: efecto Casimir.
  3. Escala de Planck: dónde la teoría semiclásica deja de valer.

Uso:  python3 simulaciones/escalas_fisicas.py
Solo usa la biblioteca estándar.
"""

import math

c = 299_792_458.0          # m/s
G = 6.674_30e-11           # m^3 kg^-1 s^-2
hbar = 1.054_571_817e-34   # J s
M_tierra = 5.972e24        # kg
R_tierra = 6.371e6         # m
g = 9.806_65               # m/s^2
DIA = 86_400.0             # s


def gps():
    r = R_tierra + 20_200e3                  # radio orbital GPS
    v = math.sqrt(G * M_tierra / r)          # velocidad orbital circular
    # Fracción de ritmo (orden 1/c^2) respecto a un reloj en la superficie
    grav = G * M_tierra / c**2 * (1 / R_tierra - 1 / r)
    cine = -v**2 / (2 * c**2)
    print("1) Relojes GPS (aproximación de campo débil, orden 1/c^2)")
    print(f"   gravitatorio : {grav * DIA * 1e6:+.1f} µs/día (adelanta)")
    print(f"   cinemático   : {cine * DIA * 1e6:+.1f} µs/día (atrasa)")
    print(f"   neto         : {(grav + cine) * DIA * 1e6:+.1f} µs/día")
    print("   (Sin esta corrección el error de posición crecería ~10 km/día.)")


def chou_2010():
    h = 0.33
    frac = g * h / c**2
    print("\n1b) Relojes ópticos separados 33 cm (Chou et al. 2010)")
    print(f"   Δf/f = g h / c^2 = {frac:.2e}")
    print(f"   equivale a {frac * 79 * 365.25 * DIA * 1e9:.0f} ns en 79 años de vida")


def casimir():
    print("\n2) Densidad de energía de Casimir entre placas ideales")
    print("   ρ = -π² ħ c / (720 a⁴)")
    for a in (1e-6, 1e-8):
        rho = -math.pi**2 * hbar * c / (720 * a**4)
        print(f"   a = {a:.0e} m  ->  ρ = {rho:.2e} J/m³  ({rho / c**2:.2e} kg/m³)")
    # Escala de violación necesaria (orden de magnitud, Morris-Thorne):
    # para una garganta de radio b, |ρ| ~ c^4 / (8π G b^2)
    b = 1.0
    rho_nec = c**4 / (8 * math.pi * G * b**2)
    print(f"   Agujero de gusano de garganta b = 1 m necesita |ρ| ~ {rho_nec:.1e} J/m³")
    rho_cas = math.pi**2 * hbar * c / (720 * 1e-6**4)
    print(f"   -> unos {math.log10(rho_nec / rho_cas):.0f} órdenes de magnitud por encima de Casimir a 1 µm.")


def planck():
    lP = math.sqrt(hbar * G / c**3)
    tP = lP / c
    EP = math.sqrt(hbar * c**5 / G)
    print("\n3) Escala de Planck (donde la gravedad semiclásica deja de valer)")
    print(f"   longitud l_P = {lP:.3e} m")
    print(f"   tiempo   t_P = {tP:.3e} s")
    print(f"   energía  E_P = {EP / 1.602_176_634e-10:.3e} GeV")


if __name__ == "__main__":
    gps()
    chou_2010()
    casimir()
    planck()
