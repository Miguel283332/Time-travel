"""
Modelos cuánticos de curvas temporales cerradas (CTC): Deutsch (D-CTC) y
post-selección (P-CTC).

Esto NO simula un viaje en el tiempo real. Calcula lo que predicen los dos
modelos operacionales más usados para un sistema cuántico que interactúa con
una copia de sí mismo "llegada del futuro". Sirve para comparar sus
predicciones (sección 6 de docs/planteamiento.md).

Convención: el espacio de Hilbert es H_A ⊗ H_CTC (primero el sistema que
respeta la cronología, A, y después el qubit que recorre la CTC).

  D-CTC (Deutsch 1991):
      rho_CTC = Tr_A[ U (rho_A ⊗ rho_CTC) U† ]          (punto fijo)
      rho_A'  = Tr_CTC[ U (rho_A ⊗ rho_CTC) U† ]

  P-CTC (Lloyd et al. 2011, teleportación post-seleccionada):
      C       = Tr_CTC[ U ]  (traza parcial del unitario, operador sobre A)
      rho_A'  = C rho_A C† / Tr(C rho_A C†)   (si el denominador es 0, el
                                                proceso tiene probabilidad 0)

Uso:  python3 simulaciones/ctc_cuantico.py
Requiere: numpy
"""

import numpy as np

I2 = np.eye(2)
X = np.array([[0, 1], [1, 0]], dtype=complex)
P0 = np.diag([1, 0]).astype(complex)
P1 = np.diag([0, 1]).astype(complex)


def kron(*ops):
    out = np.eye(1, dtype=complex)
    for op in ops:
        out = np.kron(out, op)
    return out


# Puertas sobre A ⊗ CTC
CNOT_A_CTC = kron(P0, I2) + kron(P1, X)  # A controla, CTC es el blanco
CNOT_CTC_A = kron(I2, P0) + kron(X, P1)  # CTC controla, A es el blanco
SWAP = np.array([[1, 0, 0, 0], [0, 0, 1, 0], [0, 1, 0, 0], [0, 0, 0, 1]], dtype=complex)


def traza_sobre_A(M):
    return M.reshape(2, 2, 2, 2).trace(axis1=0, axis2=2)


def traza_sobre_CTC(M):
    return M.reshape(2, 2, 2, 2).trace(axis1=1, axis2=3)


def deutsch(U, rho_A, iteraciones=20000):
    """Punto fijo de Deutsch.

    Se toma el promedio de Cesàro de la órbita que parte del estado
    máximamente mezclado. Si el punto fijo no es único, esto selecciona el que
    se alcanza desde I/2 (la prescripción de máxima entropía de Deutsch en los
    casos de este archivo).
    """
    Ud = U.conj().T
    sigma = np.eye(2, dtype=complex) / 2
    acumulado = np.zeros((2, 2), dtype=complex)
    for _ in range(iteraciones):
        acumulado += sigma
        sigma = traza_sobre_A(U @ kron(rho_A, sigma) @ Ud)
    rho_ctc = acumulado / iteraciones
    residuo = np.abs(traza_sobre_A(U @ kron(rho_A, rho_ctc) @ Ud) - rho_ctc).max()
    salida = traza_sobre_CTC(U @ kron(rho_A, rho_ctc) @ Ud)
    return salida, rho_ctc, residuo


def postseleccion(U, rho_A):
    C = traza_sobre_CTC(U)  # operador 2x2 sobre A
    salida = C @ rho_A @ C.conj().T
    peso = np.real(np.trace(salida))
    if peso < 1e-12:
        return None, peso
    return salida / peso, peso


def fmt(rho):
    if rho is None:
        return "imposible (probabilidad 0)"
    r = np.real_if_close(np.round(rho, 3))
    return np.array2string(r, separator=", ").replace("\n", "")


ket0 = np.array([1, 0], dtype=complex)
ket_mas = np.array([1, 1], dtype=complex) / np.sqrt(2)
rho_0 = np.outer(ket0, ket0.conj())
rho_mas = np.outer(ket_mas, ket_mas.conj())


def caso(titulo, explicacion, U, rho_A):
    print("=" * 72)
    print(titulo)
    print(explicacion)
    d_out, d_ctc, res = deutsch(U, rho_A)
    p_out, peso = postseleccion(U, rho_A)
    print(f"  D-CTC  estado en la CTC : {fmt(d_ctc)}   (residuo {res:.1e})")
    print(f"  D-CTC  salida de A      : {fmt(d_out)}")
    print(f"  P-CTC  salida de A      : {fmt(p_out)}   (peso sin normalizar {peso:.3f})")


def main():
    caso(
        "1) Paradoja del abuelo",
        "   U = I ⊗ X: el qubit que viaja al pasado niega su propio valor.",
        kron(I2, X),
        rho_0,
    )
    print("  -> Deutsch evita la paradoja con una mezcla 50/50 en la CTC;")
    print("     la post-selección la declara un suceso de probabilidad nula.")

    caso(
        "2) Bucle de información ('paradoja del teorema no demostrado')",
        "   U = SWAP: A entrega su estado al pasado y recibe el que llega.",
        SWAP,
        rho_mas,
    )
    print("  -> Ambos modelos son autoconsistentes: lo que circula por la CTC")
    print("     es el propio estado de A, sin que nadie lo haya 'creado'.")

    U = CNOT_CTC_A @ CNOT_A_CTC
    print("=" * 72)
    print("3) No linealidad de D-CTC")
    print("   U = CNOT(CTC→A) · CNOT(A→CTC). Se comparan la salida de una")
    print("   mezcla con la mezcla de las salidas.")
    out_0, _, _ = deutsch(U, rho_0)
    out_mas, _, _ = deutsch(U, rho_mas)
    out_mezcla, _, _ = deutsch(U, (rho_0 + rho_mas) / 2)
    diferencia = np.abs(out_mezcla - (out_0 + out_mas) / 2).max()
    print(f"  D-CTC  salida(mezcla)          : {fmt(out_mezcla)}")
    print(f"  D-CTC  mezcla de las salidas   : {fmt((out_0 + out_mas) / 2)}")
    print(f"  diferencia máxima              : {diferencia:.3f}")
    print("  -> Una mecánica cuántica lineal daría diferencia 0. Esta no")
    print("     linealidad es la fuente de las capacidades 'imposibles' de")
    print("     D-CTC (PSPACE, distinguir estados no ortogonales).")
    print("=" * 72)


if __name__ == "__main__":
    main()
