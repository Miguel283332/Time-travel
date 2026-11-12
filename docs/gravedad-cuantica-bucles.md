# Demostrar la gravedad cuántica de bucles: planteamiento

> **Objetivo.** Llevar la gravedad cuántica de bucles (LQG) desde "teoría
> candidata" hasta "teoría confirmada". Este documento define qué significa
> eso, dónde está hoy la teoría, qué hemos reproducido ya y cuál es el camino
> más directo hacia una prueba observacional.

---

## 0. Qué significa "demostrar" una teoría física

Una teoría física no se demuestra como un teorema. Se **confirma** cuando
predice algo nuevo, ese algo se mide y las alternativas no lo predicen. Para
LQG esto se descompone en cinco niveles:

| Nivel | Pregunta | Tipo de trabajo |
|---|---|---|
| **N0. Premisa** | ¿Es cuántica la gravedad? | Experimental (común a todas las teorías) |
| **N1. Consistencia** | ¿Es la teoría matemáticamente coherente y sin ambigüedades? | Matemático |
| **N2. Límite clásico** | ¿Recupera la relatividad general a escalas grandes? | Matemático y numérico |
| **N3. Predicción nueva** | ¿Predice algo medible que la relatividad general no predice? | Teórico y observacional |
| **N4. Discriminación** | ¿Esa predicción la distingue de cuerdas, conjuntos causales, etc.? | Observacional |

Solo superar N3 y N4 cuenta como confirmación. N1 y N2 son requisitos previos.

---

## 1. La teoría en una página

- **Variables.** En lugar de la métrica, se usa una conexión $SU(2)$ de
  Ashtekar–Barbero y su momento conjugado, la tríada densitizada:
  $$A^i_a = \Gamma^i_a + \gamma K^i_a, \qquad \{A^i_a(x), E^b_j(y)\} = 8\pi G\gamma\,\delta^b_a\delta^i_j\,\delta^3(x-y).$$
  $\gamma$ es el **parámetro de Barbero–Immirzi**, un número libre de la teoría.
- **Estados.** Las **redes de espín**: grafos cuyas aristas llevan espines
  $j = \tfrac12, 1, \tfrac32,\dots$ y cuyos nodos llevan "entrelazadores".
  Cada arista es un "cuanto de área"; cada nodo, un "cuanto de volumen".
- **Geometría discreta.** El área de una superficie atravesada por aristas de
  espines $j_i$ es
  $$A = 8\pi\gamma\,\ell_P^2 \sum_i \sqrt{j_i(j_i+1)},$$
  con un valor mínimo distinto de cero (el **hueco de área**).
- **Dinámica.** Canónica: la ligadura hamiltoniana de Thiemann (1998).
  Covariante: las **espumas de espín** (modelo EPRL, 2008), que dan amplitudes
  de transición entre redes de espín.

---

## 2. Estado actual por niveles

| Nivel | Conseguido | Abierto |
|---|---|---|
| **N0** | Gravedad en mecánica cuántica (experimento COW 1975), aunque como campo clásico externo. | Entrelazamiento inducido por gravedad (Bose et al., Marletto–Vedral 2017): aún no alcanzado. |
| **N1** | Cinemática rigurosa; unicidad de la representación (teorema LOST, 2006). Área y volumen bien definidos. | Ambigüedades en la ligadura hamiltoniana; relación exacta entre las versiones canónica y covariante. |
| **N2** | Las amplitudes de vértice EPRL reproducen la acción de Regge (relatividad general discretizada) en el límite de espines grandes (Barrett et al. 2009). La cosmología de bucles recupera la ecuación de Friedmann a baja densidad. | El "problema de la planitud" en espumas de espín; límite continuo y renormalización de la teoría completa. |
| **N3** | Predicciones candidatas: rebote en vez de Big Bang con huellas en el fondo cósmico de microondas; transición agujero negro → agujero blanco; corrección logarítmica $-\tfrac32\ln A$ a la entropía. | Ninguna se ha medido todavía. |
| **N4** | — | Hace falta una predicción con números propios de LQG (por ejemplo, ligados a $\rho_c$ o a $\gamma$). |

---

## 3. Lo que ya hemos reproducido

### 3.1 Espectro de área y valor de γ — [`simulaciones/lqg_espectros.py`](../simulaciones/lqg_espectros.py)

- Hueco de área: $4\sqrt3\,\pi\gamma\,\ell_P^2 \approx 5{,}17\,\ell_P^2 \approx 1{,}35\times10^{-69}\ \text{m}^2$.
- Con varias aristas, la separación entre niveles del área cae de
  $7{,}5\times10^{-2}$ a $4\times10^{-5}$ (en unidades de $8\pi\gamma\ell_P^2$)
  al pasar de áreas de ~6 a ~24 unidades: a escala macroscópica el área es
  efectivamente continua (requisito de N2).
- Valor de $\gamma$ exigido por la entropía de Bekenstein–Hawking $S=A/4\ell_P^2$:

  | Esquema de conteo | γ |
  |---|---|
  | Ashtekar–Baez–Krasnov (solo $j=\tfrac12$) | 0,127384 |
  | Domagala–Lewandowski / Meissner | **0,237533** |
  | Ghosh–Mitra | 0,274067 |

  Coinciden con los valores publicados. Que dependan del esquema muestra una
  ambigüedad real de la teoría (N1).

### 3.2 Rebote cosmológico — [`simulaciones/lqc_rebote.py`](../simulaciones/lqc_rebote.py)

Ecuación de Friedmann efectiva de la cosmología cuántica de bucles:

$$H^2 = \frac{8\pi G}{3}\,\rho\left(1-\frac{\rho}{\rho_c}\right), \qquad \rho_c = \frac{\sqrt3}{32\pi^2\gamma^3}\,\rho_{\text{Pl}} \approx 0{,}41\,\rho_{\text{Pl}}.$$

- La integración numérica coincide con la solución exacta
  $a^6 = a_B^6(1+24\pi G\rho_c t^2)$ con error $\sim10^{-12}$.
- La densidad nunca supera $\rho_c$: el universo **rebota**. Con el mismo
  estado inicial, la relatividad general llega a la singularidad.
- Interpretación: la geometría cuántica actúa como **repulsión gravitatoria
  efectiva** a densidades de Planck. Es el vínculo con
  [`planteamiento.md`](planteamiento.md): es la única fuente conocida de algo
  parecido a "materia exótica" que sale de la propia geometría.

---

## 4. Camino hacia la confirmación

### Fase A — Dominio técnico (en curso)
- ✅ Espectro de área, densidad del espectro, γ.
- ✅ Rebote efectivo de la cosmología de bucles.
- ⬜ Espectro del operador volumen para un nodo de valencia 4.
- ⬜ Amplitud de vértice EPRL y comprobación numérica de su asintótica de Regge.

### Fase B — La predicción más accesible: el fondo cósmico de microondas

> **Estado:** pasos 1–3 hechos y primer contraste con Planck completado.
> Resultados en [`resultados-cmb.md`](resultados-cmb.md): con el estado
> inicial más simple, el rebote aumenta la potencia a gran escala y Planck lo
> acota a $k_{\rm LQC} < 4{,}3\times10^{-4}$ Mpc⁻¹. Con el vacío no oscilante
> el realce es menor (cota $1{,}07\times10^{-3}$ Mpc⁻¹) o desaparece; ningún
> estado probado suprime potencia. Siguiente: estado de Ashtekar–Gupt y
> enfoque híbrido.

Es la ruta con **datos públicos ya disponibles** (Planck, y pronto LiteBIRD y
CMB-S4):

1. Calcular el espectro primordial de perturbaciones con el rebote de LQC
   (enfoques de "métrica vestida", Agullo–Ashtekar–Nelson, o "híbrido",
   Fernández-Méndez, Mena Marugán y Olmedo).
2. Propagarlo hasta el espectro angular del CMB con un código estándar
   (CLASS o CAMB).
3. Ajustarlo a los datos de Planck en multipolos bajos, donde hay anomalías
   conocidas (falta de potencia a gran escala). Ashtekar, Gupt, Jeong y
   Sreenath (2020) argumentan que LQC puede aliviarlas.
4. Producir una **predicción** para observables aún no medidos (polarización
   tipo E y B a gran escala) que un experimento futuro pueda confirmar o
   refutar. Esto es lo que cuenta como N3/N4.

### Fase C — Segunda línea: agujeros negros
Transición agujero negro → agujero blanco (Haggard–Rovelli 2015): estimar el
tiempo de vida y las señales (rayos gamma, radio) de agujeros negros
primordiales que estarían "explotando" hoy.

### Fase D — Publicación y colaboración
- Publicar resultados en arXiv (sección gr-qc) y someterlos a revisión por pares.
- Grupos activos en LQG/LQC: Penn State (Ashtekar), Louisiana State (Agullo),
  CPT Marsella (Rovelli), Perimeter Institute, Erlangen (Thiemann), Western
  (Vidotto) y, en España, el **IEM-CSIC de Madrid** (Mena Marugán, cosmología
  cuántica de bucles híbrida).

---

## 5. Advertencias

- LQG tiene problemas abiertos reales (sección 2). Parte del trabajo puede
  consistir en mostrar que alguna versión de la teoría **falla**: eso también
  es avanzar.
- Confirmar LQG no da por sí solo gravedad artificial ni máquinas del tiempo.
  Lo que aporta al programa de [`planteamiento.md`](planteamiento.md) es el
  mecanismo de repulsión efectiva y el marco para decidir si las curvas
  temporales cerradas están permitidas.
- Las referencias están citadas de memoria; conviene verificarlas.

---

## Referencias

**Fundamentos**
- Ashtekar, A. (1986). New variables for classical and quantum gravity. *Phys. Rev. Lett.* 57, 2244.
- Barbero, J. F. (1995). Real Ashtekar variables for Lorentzian signature space-times. *Phys. Rev. D* 51, 5507.
- Immirzi, G. (1997). Real and complex connections for canonical gravity. *Class. Quantum Grav.* 14, L177.
- Rovelli, C. & Smolin, L. (1995). Discreteness of area and volume in quantum gravity. *Nucl. Phys. B* 442, 593.
- Ashtekar, A. & Lewandowski, J. (1997). Quantum theory of geometry I: Area operators. *Class. Quantum Grav.* 14, A55.
- Thiemann, T. (1998). Quantum spin dynamics (QSD). *Class. Quantum Grav.* 15, 839.
- Lewandowski, J., Okołów, A., Sahlmann, H. & Thiemann, T. (2006). Uniqueness of diffeomorphism invariant states on holonomy-flux algebras. *Commun. Math. Phys.* 267, 703.
- Engle, J., Pereira, R., Rovelli, C. & Livine, E. (2008). LQG vertex with finite Immirzi parameter. *Nucl. Phys. B* 799, 136.
- Barrett, J. W. et al. (2009). Asymptotic analysis of the EPRL four-simplex amplitude. *J. Math. Phys.* 50, 112504.

**Libros**
- Rovelli, C. & Vidotto, F. (2014). *Covariant Loop Quantum Gravity*. Cambridge University Press.
- Thiemann, T. (2007). *Modern Canonical Quantum General Relativity*. Cambridge University Press.

**Agujeros negros**
- Ashtekar, A., Baez, J., Corichi, A. & Krasnov, K. (1998). Quantum geometry and black hole entropy. *Phys. Rev. Lett.* 80, 904.
- Kaul, R. K. & Majumdar, P. (2000). Logarithmic correction to the Bekenstein-Hawking entropy. *Phys. Rev. Lett.* 84, 5255.
- Domagala, M. & Lewandowski, J. (2004). *Class. Quantum Grav.* 21, 5233.
- Meissner, K. A. (2004). *Class. Quantum Grav.* 21, 5245.
- Ghosh, A. & Mitra, P. (2005). *Phys. Lett. B* 616, 114.
- Rovelli, C. & Vidotto, F. (2014). Planck stars. *Int. J. Mod. Phys. D* 23, 1442026.
- Haggard, H. M. & Rovelli, C. (2015). Black hole fireworks. *Phys. Rev. D* 92, 104020.

**Cosmología cuántica de bucles**
- Ashtekar, A., Pawlowski, T. & Singh, P. (2006). Quantum nature of the big bang. *Phys. Rev. Lett.* 96, 141301; *Phys. Rev. D* 74, 084003.
- Ashtekar, A. & Singh, P. (2011). Loop quantum cosmology: a status report. *Class. Quantum Grav.* 28, 213001.
- Agullo, I., Ashtekar, A. & Nelson, W. (2012). Quantum gravity extension of the inflationary scenario. *Phys. Rev. Lett.* 109, 251301.
- Fernández-Méndez, M., Mena Marugán, G. A. & Olmedo, J. (2012). Hybrid quantization of an inflationary universe. *Phys. Rev. D* 86, 024003.
- Ashtekar, A., Gupt, B., Jeong, D. & Sreenath, V. (2020). Alleviating the tension in the cosmic microwave background using Planck-scale physics. *Phys. Rev. Lett.* 125, 051302.

- Martín de Blas, D. & Olmedo, J. (2016). Primordial power spectra for scalar perturbations in loop quantum cosmology. *JCAP* 06, 029.
- Ashtekar, A. & Gupt, B. (2017). Quantum gravity in the sky: interplay between fundamental theory and observations. *Class. Quantum Grav.* 34, 014002.

**Fenomenología y experimentos**
- Rovelli, C. & Speziale, S. (2003). Reconcile Planck-scale discreteness and the Lorentz-Fitzgerald contraction. *Phys. Rev. D* 67, 064019.
- Colella, R., Overhauser, A. W. & Werner, S. A. (1975). *Phys. Rev. Lett.* 34, 1472.
- Bose, S. et al. (2017). *Phys. Rev. Lett.* 119, 240401.
- Marletto, C. & Vedral, V. (2017). *Phys. Rev. Lett.* 119, 240402.
