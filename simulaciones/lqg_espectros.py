"""
Gravedad cuántica de bucles (LQG): resultados cinemáticos básicos.

  1. Espectro discreto del operador área y el "hueco de área".
  2. Cómo el espectro se vuelve casi continuo al crecer el área
     (condición necesaria para recuperar la relatividad general).
  3. El parámetro de Barbero-Immirzi γ fijado por la entropía de
     Bekenstein-Hawking S = A / (4 ℓ_P²), según tres esquemas de conteo.

Todo en unidades de Planck (ℓ_P = 1). Ver docs/gravedad-cuantica-bucles.md.

Uso:  python3 simulaciones/lqg_espectros.py
Solo usa la biblioteca estándar.
"""

import math

GAMMA_REF = 0.2375329  # valor de Meissner / Domagala-Lewandowski


def area_un_pinchazo(j, gamma=GAMMA_REF):
    """A_j = 8 π γ ℓ_P² sqrt(j (j+1)) para una arista de espín j."""
    return 8 * math.pi * gamma * math.sqrt(j * (j + 1))


def espines(jmax):
    return [k / 2 for k in range(1, int(2 * jmax) + 1)]


def espectro():
    print("1) Espectro del área (un solo pinchazo), γ = %.7f" % GAMMA_REF)
    for j in espines(3):
        print(f"   j = {j:>3}:  A = {area_un_pinchazo(j):8.3f} ℓ_P²")
    gap = area_un_pinchazo(0.5)
    print(f"   Hueco de área (j = 1/2): {gap:.3f} ℓ_P² = {gap * 2.612e-70:.2e} m²")
    print("   -> No existe ninguna superficie con área entre 0 y este valor.")


def densidad_espectro(cota=24.0):
    """Cuenta los autovalores distintos de un área con varios pinchazos.

    Un área con pinchazos de espines j_1..j_n vale 8πγ Σ sqrt(j_i(j_i+1)).
    Se cuentan los valores distintos de la suma (en unidades de 8πγ) hasta
    'cota' y se mide su separación media en cada tramo.
    """
    pesos = [math.sqrt(j * (j + 1)) for j in espines(cota)]
    valores = {0.0}
    frontera = {0.0}
    while frontera:
        nueva = set()
        for v in frontera:
            for w in pesos:
                s = round(v + w, 9)
                if s <= cota and s not in valores:
                    nueva.add(s)
        valores |= nueva
        frontera = nueva
    orden = sorted(valores)
    print("\n2) Densidad del espectro con varios pinchazos (área en unidades de 8πγ ℓ_P²)")
    print("   tramo de área     autovalores distintos   separación media")
    paso = cota / 4
    for i in range(4):
        lo, hi = i * paso, (i + 1) * paso
        trozo = [v for v in orden if lo < v <= hi]
        sep = (hi - lo) / len(trozo) if trozo else float("nan")
        print(f"   ({lo:5.1f}, {hi:5.1f}]   {len(trozo):>12}            {sep:.2e}")
    print("   -> La separación cae muy deprisa: a escalas macroscópicas el área")
    print("      parece continua, como en la relatividad general.")


def biseccion(f, a, b, tol=1e-12):
    fa = f(a)
    for _ in range(200):
        m = (a + b) / 2
        fm = f(m)
        if (fa > 0) == (fm > 0):
            a, fa = m, fm
        else:
            b = m
        if b - a < tol:
            break
    return (a + b) / 2


def immirzi():
    """γ tal que la entropía contada sea A/4.

    Cada esquema de conteo da una ecuación del tipo
        1 = Σ_j g(j) exp(-2π γ sqrt(j(j+1)))
    con distinta multiplicidad g(j) por pinchazo.
    """
    J = espines(400)

    def ecuacion(g):
        return lambda gam: sum(g(j) * math.exp(-2 * math.pi * gam * math.sqrt(j * (j + 1))) for j in J) - 1

    esquemas = [
        ("Solo j = 1/2, m = ±1/2 (Ashtekar-Baez-Krasnov)", None),
        ("Todos los j, m = ±j (Domagala-Lewandowski, Meissner)", lambda j: 2.0),
        ("Todos los j, todos los m (Ghosh-Mitra)", lambda j: 2 * j + 1),
    ]
    print("\n3) Parámetro de Barbero-Immirzi fijado por S = A/4")
    for nombre, g in esquemas:
        if g is None:
            gam = math.log(2) / (math.pi * math.sqrt(3))
        else:
            gam = biseccion(ecuacion(g), 0.05, 1.0)
        print(f"   {nombre:<55} γ = {gam:.6f}")
    print("   -> El valor depende de cómo se cuentan los estados del horizonte.")
    print("      Fijar el esquema correcto es un problema abierto de la teoría.")


if __name__ == "__main__":
    espectro()
    densidad_espectro()
    immirzi()
