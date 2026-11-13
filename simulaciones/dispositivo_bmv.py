"""
Números de diseño para un dispositivo de laboratorio que ponga a prueba la
naturaleza cuántica de la gravedad (docs/plan-dispositivo.md).

  1. Las tres cosas que se llaman "gravedad artificial" y su tamaño real.
  2. Por qué el laboratorio no ve la longitud de Planck pero sí puede
     acercarse a la masa de Planck.
  3. Experimento de entrelazamiento por gravedad (Bose et al. 2017;
     Marletto y Vedral 2017): fase, entrelazamiento, fuerza de Casimir-Polder
     que hay que evitar y vacío necesario.
  4. Distancia entre el estado del arte y lo que hace falta.
  5. Régimen en el que una discretización del tiempo a escala de Planck
     sería visible (Christodoulou, Di Biagio y Martin-Dussaud 2022).

Todas las cifras son estimaciones de orden de magnitud con fórmulas estándar.

Uso:  python3 simulaciones/dispositivo_bmv.py
Solo usa la biblioteca estándar.
"""

import math

G = 6.674_30e-11
c = 299_792_458.0
hbar = 1.054_571_817e-34
kB = 1.380_649e-23
u = 1.660_539e-27

l_P = math.sqrt(hbar * G / c**3)
t_P = l_P / c
m_P = math.sqrt(hbar * c / G)


def seccion(titulo):
    print("\n" + "=" * 74)
    print(titulo)
    print("=" * 74)


def gravedad_artificial():
    seccion("1) Tres significados de 'gravedad artificial'")
    # (a) Rotación
    r = 10.0
    w = math.sqrt(9.81 / r)
    print(f"a) Rotación: 1 g con radio {r:.0f} m exige {w:.2f} rad/s ({w * 60 / (2 * math.pi):.1f} rpm).")
    print("   Es inercia (principio de equivalencia); no dice nada de gravedad cuántica.")
    # (b) Masa fuente
    m, d = 90e-6, 2.5e-3
    print(f"b) Masa fuente: una esfera de oro de {m * 1e6:.0f} mg a {d * 1e3:.1f} mm produce "
          f"a = {G * m / d**2:.1e} m/s²")
    print("   (medido por Westphal et al. 2021). Esto SÍ es crear un campo gravitatorio.")
    # (c) Arrastre de sistemas de referencia en el laboratorio
    M, R, omega = 1000.0, 0.5, 1000.0
    L = 0.5 * M * R**2 * omega
    lab = G * L / (c**2 * R**3)
    gpb = 39e-3 / 3600 * math.pi / 180 / (365.25 * 86400)
    print(f"c) Campo 'gravitomagnético' de un rotor de {M:.0f} kg, R = {R} m, {omega:.0f} rad/s:")
    print(f"   Ω ~ G L / (c² r³) = {lab:.1e} rad/s, {gpb / lab:.0e} veces menor que el de la")
    print(f"   Tierra medido por Gravity Probe B ({gpb:.1e} rad/s, con un 19 % de error).")
    print("   No existe ningún 'generador de gravedad' sin masa en ninguna teoría aceptada.")


def escalas():
    seccion("2) Longitud de Planck frente a masa de Planck")
    Lab = 1e-6
    print(f"   l_P = {l_P:.2e} m   t_P = {t_P:.2e} s   m_P = {m_P:.2e} kg ({m_P * 1e9:.1f} µg)")
    print(f"   Corrección por área discreta a escala de 1 µm: (l_P/L)² = {(l_P / Lab)**2:.0e}")
    print("   -> la discretización del espacio es inaccesible de forma directa.")
    print(f"   Pero m_P es mesoscópica: un grano de {m_P * 1e9:.0f} µg. La gravedad cuántica")
    print("   se vuelve visible cuando una masa cercana a m_P se comporta cuánticamente.")


def bmv(m=1e-14, dx=250e-6, d=450e-6, tau=2.5, rho=3500.0, eps=5.7, T=0.15, verbose=True):
    """Dos masas m, cada una en superposición de dos posiciones separadas dx,
    con centros a distancia d, durante un tiempo tau. Geometría colineal."""
    r = {"LL": d, "RR": d, "LR": d + dx, "RL": d - dx}
    phi = {k: G * m * m * tau / (hbar * v) for k, v in r.items()}
    delta = phi["LR"] + phi["RL"] - phi["LL"] - phi["RR"]
    concurrencia = abs(math.sin(delta / 2))
    R = (3 * m / (4 * math.pi * rho)) ** (1 / 3)
    rmin = r["RL"]
    V_g = G * m * m / rmin
    V_cp = 23 * hbar * c / (4 * math.pi) * ((eps - 1) / (eps + 2)) ** 2 * R**6 / rmin**7
    # Gas residual (helio): cada colisión delata la posición
    v_gas = math.sqrt(8 * kB * T / (math.pi * 4 * u))
    n_max = 0.1 / (tau * v_gas * math.pi * R**2)  # 0.1 colisiones en tau
    P_max = n_max * kB * T
    if verbose:
        print(f"   m = {m:.0e} kg (R = {R * 1e6:.2f} µm de diamante), Δx = {dx * 1e6:.0f} µm, "
              f"d = {d * 1e6:.0f} µm, τ = {tau} s")
        for k in ("RL", "LL", "LR"):
            print(f"     fase φ_{k} (r = {r[k] * 1e6:.0f} µm) = {phi[k]:.3f} rad")
        print(f"   Fase de entrelazamiento Δ = {delta:.3f} rad  ->  concurrencia = {concurrencia:.3f}")
        print(f"   Casimir-Polder / gravedad a {rmin * 1e6:.0f} µm = {V_cp / V_g:.3f}")
        print(f"   Vacío para < 0.1 colisiones en τ (He a {T} K): n < {n_max:.1e} m⁻³, "
              f"P < {P_max:.0e} Pa")
    return delta, concurrencia, V_cp / V_g


def experimento():
    seccion("3) Experimento de entrelazamiento por gravedad (parámetros de Bose et al.)")
    bmv()
    print("\n   Sensibilidad a la masa (misma geometría):")
    print("        m [kg]     Δ [rad]   concurrencia")
    for m in (1e-16, 1e-15, 1e-14, 3e-14):
        delta, conc, _ = bmv(m=m, verbose=False)
        print(f"       {m:.0e}   {delta:9.2e}   {conc:.2e}")
    print("   -> la señal crece como m²: por debajo de ~1e-15 kg es indetectable.")


def distancia():
    seccion("4) Distancia al estado del arte")
    m_record = 25_000 * u
    print(f"   Mayor masa en superposición espacial (Fein et al. 2019): ~25 kDa = {m_record:.1e} kg")
    print(f"   Masa necesaria: 1e-14 kg  ->  faltan {math.log10(1e-14 / m_record):.1f} órdenes de magnitud")
    print("   Separación de la superposición: cientos de nm con moléculas; hacen falta ~250 µm.")
    print("   Nanopartícula levitada en el estado fundamental (Delić et al. 2020): ~1e-18 kg,")
    print("   localizada (sin superposición espacial).")


def tiempo_discreto():
    seccion("5) ¿Algo específico de la gravedad cuántica de bucles? Tiempo discreto")
    print("   La fase de entrelazamiento es φ = m c² δτ / ħ, con δτ la diferencia de")
    print("   tiempo propio entre ramas. Si δτ solo pudiera valer múltiplos de t_P,")
    print("   φ saltaría en pasos de m c² t_P / ħ = m / m_P.")
    for m in (1e-14, 1e-11, m_P):
        print(f"     m = {m:.1e} kg  ->  paso de fase = {m / m_P:.1e} rad")
    print("   Con masas menores el escalón existe pero exige medir la fase con esa precisión.")
    print(f"   Escalones de ~1 rad requieren masas ~m_P = {m_P:.1e} kg en superposición:")
    print(f"   {math.log10(m_P / 1e-14):.1f} órdenes por encima del experimento de la sección 3.")
    print("   Aviso: que el tiempo propio sea discreto es una hipótesis inspirada en LQG,")
    print("   no una predicción firme de la teoría.")


if __name__ == "__main__":
    gravedad_artificial()
    escalas()
    experimento()
    distancia()
    tiempo_discreto()
