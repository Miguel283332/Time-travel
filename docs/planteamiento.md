# Viajes en el tiempo y gravedad cuántica: planteamiento y base

> **Tesis.** La relatividad general clásica admite viajes al pasado (curvas
> temporales cerradas). La física semiclásica que podría prohibirlos deja de
> ser válida justo donde se forma la "máquina del tiempo". Por eso la
> pregunta de si los viajes al pasado son físicamente posibles la decide
> necesariamente una teoría de **gravedad cuántica**. El proyecto consiste en
> plantear esa pregunta de forma que una teoría así, junto con experimentos,
> pueda responderla.

---

## 0. Estado de la cuestión (con honestidad)

| Afirmación | Estado a fecha de este documento |
|---|---|
| Viajar **al futuro** más deprisa que otro observador (dilatación temporal) | **Demostrado** y medido a diario (GPS, relojes atómicos, muones). |
| La relatividad general tiene soluciones con viaje **al pasado** | **Demostrado** matemáticamente (Gödel 1949 y muchas más). |
| Esas soluciones son físicamente realizables | **No demostrado.** Requieren condiciones (energía negativa a gran escala, cilindros infinitos, etc.) no observadas. |
| La gravedad es cuántica | **No demostrado.** Hay experimentos en marcha (sección 8). |
| Existe una teoría de gravedad cuántica completa y confirmada | **No.** Hay candidatas (cuerdas/AdS-CFT, lazos, conjuntos causales, triangulaciones dinámicas causales, seguridad asintótica). |
| La gravedad cuántica permite viajar al pasado | **Pregunta abierta.** La opinión mayoritaria se inclina por la "protección de la cronología", pero no está demostrada. |

Consecuencia: el objetivo realista **no** es "demostrar que los viajes en el
tiempo existen" sino construir el programa que **decida** si existen. Un
resultado negativo (la gravedad cuántica los prohíbe) es tan valioso como uno
positivo.

---

## 1. Definiciones

Trabajamos en un espacio-tiempo $(M, g_{\mu\nu})$ con signatura $(-,+,+,+)$
que cumple las ecuaciones de Einstein

$$
G_{\mu\nu} + \Lambda g_{\mu\nu} = \frac{8\pi G}{c^4}\, T_{\mu\nu}.
$$

- **Curva temporal:** curva cuyo vector tangente $u^\mu$ cumple
  $g_{\mu\nu}u^\mu u^\nu < 0$ en todo punto (la trayectoria posible de un
  observador con masa).
- **Curva temporal cerrada (CTC):** curva temporal que vuelve a su punto de
  partida. Un observador que la recorre se encuentra con su propio pasado.
  **"Viaje al pasado" ≡ existencia de CTC.**
- **Región de violación de la cronología:** conjunto de puntos por los que
  pasa alguna CTC.
- **Horizonte de cronología:** frontera de esa región; es un horizonte de
  Cauchy $H^+(\Sigma)$ de alguna hipersuperficie inicial $\Sigma$. Se dice
  **generado de forma compacta** si sus generadores nulos, hacia el pasado,
  entran y permanecen en una región compacta (la máquina del tiempo "se
  fabrica" en una zona finita del espacio y no viene "del infinito").
- **Condiciones de energía** (para todo vector nulo $k^\mu$ y temporal $t^\mu$):
  - Nula (NEC): $T_{\mu\nu}k^\mu k^\nu \ge 0$.
  - Débil (WEC): $T_{\mu\nu}t^\mu t^\nu \ge 0$.
  - Nula promediada (ANEC): $\int_\gamma T_{\mu\nu}k^\mu k^\nu\, d\lambda \ge 0$
    a lo largo de geodésicas nulas completas $\gamma$.
  - ANEC acronal: la ANEC exigida solo en geodésicas nulas acronales (que no
    tienen dos puntos conectados por una curva temporal).

---

## 2. Base demostrada: el viaje al futuro

El tiempo propio de un reloj en un campo gravitatorio débil y en movimiento es

$$
d\tau \simeq \left(1 - \frac{GM}{rc^2} - \frac{v^2}{2c^2}\right) dt .
$$

| Experimento | Resultado |
|---|---|
| Hafele–Keating (1971) | Relojes de cesio en aviones comerciales, diferencias de decenas a centenares de ns, de acuerdo con la predicción. |
| Satélites GPS | +45,7 µs/día (gravedad) − 7,2 µs/día (velocidad) = **+38,5 µs/día**; sin corregirlo el error crecería ~10 km/día. |
| Chou et al. (2010) | Relojes ópticos separados **33 cm** de altura: $\Delta f/f \approx gh/c^2 \approx 3.6\times10^{-17}$. |
| Bothwell et al. (2022) | Corrimiento al rojo gravitatorio medido dentro de una sola muestra atómica, a escala milimétrica. |
| Muones cósmicos y de acelerador | Vidas medias dilatadas por el factor $\gamma$. |

Los números se reproducen con [`simulaciones/escalas_fisicas.py`](../simulaciones/escalas_fisicas.py).

El viaje al futuro es, por tanto, un hecho. Todo el problema está en el
**pasado**.

---

## 3. Base clásica: la relatividad general admite CTC

| Solución | Año | Mecanismo | Pega física |
|---|---|---|---|
| Universo de Gödel | 1949 | Universo de polvo en rotación global | Nuestro universo no rota así (límites del CMB). |
| Cilindro de van Stockum / Tipler | 1937 / 1974 | Cilindro masivo en rotación rápida | Necesita longitud infinita; con longitud finita, Hawking muestra que hace falta energía negativa. |
| Interior de Kerr ($r<0$) | 1963 | Agujero negro en rotación | Más allá del horizonte de Cauchy interior, que es inestable (inflación de masa). |
| Agujero de gusano de Morris–Thorne–Yurtsever | 1988 | Desplazar o poner en un pozo gravitatorio una de las bocas | Requiere violar la NEC (materia "exótica"). |
| Cuerdas cósmicas de Gott | 1991 | Dos cuerdas infinitas cruzándose a gran velocidad | Momento total taquiónico; no puede formarse en un universo abierto. |
| Burbuja de Alcubierre (+ Everett) | 1994 / 1996 | Dos burbujas superlumínicas forman una CTC | Requiere energía negativa enorme. |
| Máquina de Ori | 2005 | Núcleo de vacío, CTC formadas en región acotada | Necesita una configuración de materia muy particular; estabilidad no demostrada. |

### 3.1 El ejemplo canónico: agujero de gusano atravesable

$$
ds^2 = -e^{2\Phi(r)}c^2dt^2 + \frac{dr^2}{1 - b(r)/r} + r^2 d\Omega^2 .
$$

- $b(r)$ es la función de forma; en la garganta $b(r_0)=r_0$.
- Para que la garganta se "abra" (condición de *flare-out*, $b'(r_0)<1$), la
  ecuación de Raychaudhuri para haces nulos obliga a
  $T_{\mu\nu}k^\mu k^\nu < 0$ en la garganta: **violación de la NEC**.
- **Cómo se convierte en máquina del tiempo:** si una boca se mueve a velocidad
  relativista y vuelve (paradoja de los gemelos), o se coloca en un pozo
  gravitatorio, el tiempo de las dos bocas se desincroniza un $\Delta T$.
  Cuando $\Delta T$ supera el tiempo que tarda la luz en ir de una boca a la
  otra por fuera, aparecen CTC.

Orden de magnitud: una garganta de $b \sim 1\ \text{m}$ requiere densidades
$|\rho| \sim c^4/(8\pi G b^2) \sim 10^{42}\ \text{J/m}^3$, frente a
$\sim 10^{-4}\ \text{J/m}^3$ del efecto Casimir con placas a 1 µm: unos
**46 órdenes de magnitud** de diferencia (ver el script).

---

## 4. Los obstáculos: teoremas y el "puente" hacia la gravedad cuántica

### 4.1 Teoremas clásicos

- **Tipler (1976-77):** bajo la WEC y condiciones genéricas, crear CTC en una
  región finita produce singularidades.
- **Hawking (1992):** un horizonte de cronología generado de forma compacta
  **exige violar la condición de energía** (los generadores nulos tendrían que
  desenfocarse, lo que la NEC impide).
- **Conjetura de protección de la cronología (Hawking 1992):** "las leyes de la
  física no permiten la aparición de curvas temporales cerradas".

### 4.2 La mecánica cuántica abre la puerta... a medias

La teoría cuántica de campos **sí** permite energía negativa local (Casimir,
estados comprimidos, efecto Hawking). Pero la limita:

- **Desigualdades cuánticas (Ford–Roman):** cuanto más negativa la energía,
  menos tiempo puede durar.
- **ANEC demostrada** en teoría cuántica de campos en espacio plano
  (Faulkner et al. 2016; Hartman et al. 2017).
- **ANEC acronal (Graham–Olum 2007):** si se cumple en gravedad semiclásica,
  **prohíbe** agujeros de gusano atravesables y máquinas del tiempo.
  No está demostrada en espacio curvo general.

### 4.3 El puente: la teoría semiclásica se rompe justo en el horizonte

En gravedad semiclásica se resuelve
$G_{\mu\nu} = 8\pi G\, \langle \hat T_{\mu\nu} \rangle / c^4$.

- **Kim–Thorne (1991), Hawking (1992):** al acercarse al horizonte de
  cronología, las fluctuaciones del vacío recorren el bucle una y otra vez y
  $\langle \hat T_{\mu\nu}\rangle$ crece sin límite (efecto "realimentación").
- **Contraejemplos:** Li y Gott (1998) encontraron en el espacio de Misner un
  estado con $\langle \hat T_{\mu\nu}\rangle$ regular. La divergencia no es
  universal.
- **Teorema de Kay–Radzikowski–Wald (1997):** en un horizonte de cronología
  generado de forma compacta existen puntos donde **ninguna** función de dos
  puntos de tipo Hadamard está bien definida. Es decir, el propio objeto
  $\langle \hat T_{\mu\nu}\rangle$ deja de tener sentido allí.

**Conclusión de la sección (el núcleo del planteamiento):**

> Ni la relatividad general clásica (que permite CTC) ni la teoría
> semiclásica (que deja de estar definida en el horizonte de cronología)
> pueden dar el veredicto. Las distancias relevantes en el horizonte llegan a
> la escala de Planck ($\ell_P \approx 1.6\times10^{-35}$ m). **La respuesta
> sólo puede venir de la gravedad cuántica.**

---

## 5. Qué tiene que responder una teoría de gravedad cuántica

Sea $\mathcal{Q}$ una teoría candidata. Para decidir la cuestión debe
responder:

- **(Q1) Dominio.** ¿Incluye la teoría estados o geometrías con estructura
  causal no trivial? (¿suma la integral de caminos sobre geometrías con CTC?
  ¿hay cambio de topología?)
- **(Q2) Retroacción.** ¿Es finita la energía efectiva en el horizonte de
  cronología una vez incluidos los efectos a escala de Planck?
- **(Q3) Consistencia.** ¿Qué condición de autoconsistencia aplica: la de
  Novikov, la de Deutsch (D-CTC) o la de post-selección (P-CTC)? (sección 6)
- **(Q4) Materia exótica.** ¿Permite la completitud ultravioleta violaciones
  de la NEC suficientemente grandes y duraderas?
- **(Q5) Causalidad: ¿fundamental o emergente?** Si la causalidad es un
  axioma, la prohibición no se *demuestra*, se *supone*.

### 5.1 Qué dicen hoy las candidatas

| Enfoque | Posición respecto a las CTC |
|---|---|
| **Cuerdas / supergravedad** | Existen soluciones supersimétricas tipo Gödel (Gauntlett et al. 2003). Sin embargo, el análisis holográfico sugiere que la "pantalla holográfica" deja fuera las CTC (Boyda, Ganguli, Hořava, Varadarajan 2003). |
| **AdS/CFT y agujeros de gusano** | Gao–Jafferis–Wall (2017) construyen agujeros de gusano **atravesables** con un acoplamiento que viola la ANEC de forma controlada, pero el atajo es más lento que la comunicación por la frontera: **no hay violación de causalidad**. ER=EPR (Maldacena–Susskind 2013) relaciona agujeros de gusano y entrelazamiento. |
| **Gravedad cuántica de lazos / espumas de espín** | No hay resultado definitivo sobre CTC. Pregunta abierta. |
| **Conjuntos causales (Sorkin et al.)** | El espacio-tiempo es un orden parcial **acíclico**: las CTC quedan excluidas por axioma. Observar una CTC refutaría el enfoque. |
| **Triangulaciones dinámicas causales (CDT)** | Imponen una foliación global en el tiempo; excluyen CTC y cambio de topología por construcción. |
| **Integral de caminos euclídea** | Permite cambio de topología; las CTC lorentzianas no aparecen de forma natural. Abierto. |
| **Gravedad clásica-estocástica (Oppenheim)** | Alternativa a cuantizar la gravedad; si fuera correcta, el veredicto lo daría esa teoría híbrida. Hace predicciones comprobables (compromiso decoherencia-difusión). |

---

## 6. Modelos cuánticos de CTC (nivel operacional)

Aunque no sepamos la teoría completa, la mecánica cuántica ya ofrece modelos de
cómo se comportaría un sistema cuántico que interactúa con su versión
"llegada del futuro".

**Deutsch (D-CTC, 1991).** El sistema $A$ interactúa mediante $U$ con un
sistema que recorre la CTC. Se exige consistencia:

$$
\rho_{\text{CTC}} = \operatorname{Tr}_A\!\left[U(\rho_A\otimes\rho_{\text{CTC}})U^\dagger\right],
\qquad
\rho_A' = \operatorname{Tr}_{\text{CTC}}\!\left[U(\rho_A\otimes\rho_{\text{CTC}})U^\dagger\right].
$$

Siempre existe un punto fijo, de modo que la paradoja del abuelo se resuelve
con una mezcla (el abuelo muere "con probabilidad 1/2"). Precio: la evolución
es **no lineal**; permite resolver problemas PSPACE (Aaronson–Watrous 2009) y
distinguir estados no ortogonales.

**Post-selección (P-CTC, Lloyd et al. 2011).** Equivale a teleportación
cuántica post-seleccionada:

$$
\rho_A' = \frac{C\rho_A C^\dagger}{\operatorname{Tr}(C\rho_A C^\dagger)},
\qquad C = \operatorname{Tr}_{\text{CTC}}(U).
$$

Las paradojas tienen probabilidad cero: simplemente no ocurren. Potencia
computacional: PP.

**Novikov (Friedman et al. 1990).** Clásico: sólo ocurren las historias
globalmente autoconsistentes (problema de la bola de billar de Echeverria,
Klinkhammer y Thorne).

**Experimentos.** Lloyd et al. (2011) y Ringbauer et al. (2014) **simularon**
P-CTC y D-CTC con fotones. Son simulaciones, no CTC reales, pero demuestran
que los modelos son operacionalmente distintos y comprobables.

El script [`simulaciones/ctc_cuantico.py`](../simulaciones/ctc_cuantico.py)
calcula ambos modelos para tres circuitos: la paradoja del abuelo (Deutsch da
una mezcla 50/50; la post-selección da probabilidad 0), un bucle de
información y un ejemplo de no linealidad de D-CTC.

**Para el programa:** una teoría de gravedad cuántica que admita CTC debe
**derivar** cuál de estas reglas de consistencia (u otra) se cumple. Eso
convierte (Q3) en una predicción contrastable.

---

## 7. Estructura lógica de la demostración

Definimos las proposiciones:

- **A** — Una teoría $\mathcal{Q}$ de gravedad cuántica está confirmada
  experimentalmente.
- **B** — $\mathcal{Q}$ implica que existen evoluciones que, partiendo de datos
  iniciales **sin** CTC, generan CTC con energía efectiva finita en el
  horizonte de cronología.
- **C** — Las condiciones que exige B (materia, energía, configuración) existen
  en la naturaleza o pueden construirse.
- **T** — Los viajes al pasado son físicamente posibles.

Entonces:

$$
A \wedge B \wedge C \;\Rightarrow\; T
\qquad\qquad
A \wedge \neg B \;\Rightarrow\; \neg T\ \ \text{(protección de la cronología demostrada)}
$$

**Demostrar la gravedad cuántica (A) es necesario pero no suficiente.** Hacen
falta también B (un teorema dentro de la teoría) y C (física realizable). La
alternativa a todo esto es la observación directa de una CTC, que no tiene
precedentes.

Tres escenarios posibles al final del programa:

1. **Protección fuerte:** $\mathcal{Q}$ prohíbe las CTC (por ejemplo, por una
   versión cuántica de la ANEC acronal). Resultado: $\neg T$.
2. **Protección débil:** sólo hay CTC "virtuales" a escala de Planck
   (fluctuaciones de la espuma espacio-temporal), sin CTC macroscópicas
   utilizables. Resultado: $T$ en sentido técnico, sin viajes utilizables.
3. **CTC macroscópicas:** $\mathcal{Q}$ las permite y C se cumple. Resultado: $T$.

---

## 8. Hipótesis de trabajo y pruebas falsables

| # | Hipótesis | Prueba | Apoya si... | Refuta si... | Estado |
|---|---|---|---|---|---|
| H1 | La gravedad es un mediador cuántico | Entrelazamiento inducido por gravedad entre dos masas en superposición (Bose et al. 2017; Marletto–Vedral 2017) | Aparece entrelazamiento | No aparece con sensibilidad suficiente | Propuesto; aún no alcanzado |
| H2 | La estructura causal puede estar en superposición | "Interruptor cuántico gravitatorio": una masa en superposición de posiciones produce orden temporal indefinido (Zych, Costa, Pikovski, Brukner 2019) | Se viola la desigualdad de orden temporal | Se cumple | Teórico; el interruptor cuántico ya se ha realizado con fotones, sin gravedad |
| H3 | El tiempo propio se comporta cuánticamente | Interferometría de relojes en superposición de alturas (Margalit et al. 2015, reloj autointerferente) | Pérdida de visibilidad ∝ diferencia de tiempo propio | No hay pérdida | Demostraciones parciales |
| H4 | Se cumple la ANEC acronal en gravedad cuántica | Teórica (holografía, demostraciones en espacio curvo) | Demostración general | Contraejemplo controlado | Abierto |
| H5 | Hay invariancia de Lorentz a escala de Planck | Dispersión en fotones de estallidos de rayos gamma (Fermi, GRB 090510) | Sin dispersión | Dispersión dependiente de la energía | Sin violación detectada (más allá de $E_P$ a primer orden) |
| H6 | Existen agujeros de gusano astrofísicos | Microlente con masa negativa, ecos en ondas gravitacionales, sombras (EHT) | Señal inequívoca | No detección (sólo acota) | Sin detección |
| H7 | La regla de consistencia es D-CTC o P-CTC | Derivarla de $\mathcal{Q}$; contrastar con simulaciones fotónicas | Coincidencia con la teoría | Contradicción | Simulaciones hechas; derivación pendiente |

H1–H3 atacan **A**; H4 ataca **B**; H6 ataca **C**; H7 conecta los modelos
operacionales con la teoría fundamental.

---

## 9. Hoja de ruta

**Fase 1 — Base teórica clásica (este repositorio).**
Catálogo de soluciones con CTC; reproducción de los teoremas de Tipler y
Hawking; scripts de modelos D-CTC y P-CTC.

**Fase 2 — Semiclásica.**
Calcular $\langle \hat T_{\mu\nu}\rangle$ en modelos de juguete: espacio de
Misner, agujero de gusano en 1+1 dimensiones. Comprobar cuándo diverge y cuándo
no (Li–Gott). Estudiar la ANEC acronal en esos fondos.

**Fase 3 — Modelos de gravedad cuántica calculables.**
- Gravedad en 2+1 dimensiones y AdS₃/CFT₂ (agujero BTZ en rotación, que tiene
  regiones con CTC si no se recortan).
- Gravedad de Jackiw–Teitelboim (2D), donde la integral de caminos es
  calculable: ¿contribuyen geometrías con CTC?
- Agujeros de gusano atravesables tipo Gao–Jafferis–Wall: ¿puede la
  construcción producir CTC si se añade movimiento relativo de las bocas?

**Fase 4 — Seguimiento experimental.**
Seguir y, si es posible, colaborar con los experimentos de H1–H3 y H5–H6.

**Fase 5 — Criterios de éxito.**
- Demostrar o refutar (Q2)–(Q4) en al menos un modelo calculable.
- Identificar una predicción observable que distinga los escenarios 1, 2 y 3
  de la sección 7.

---

## 10. Conclusión

1. El viaje **al futuro** está demostrado.
2. El viaje **al pasado** está permitido por la relatividad general clásica,
   pero cada solución conocida exige condiciones no observadas.
3. La teoría semiclásica no puede decidir la cuestión porque deja de estar
   definida en el horizonte de cronología (Kay–Radzikowski–Wald).
4. Por tanto **la cuestión es de gravedad cuántica**, y demostrar la gravedad
   cuántica es un paso necesario, no suficiente.
5. El programa es falsable: puede acabar demostrando la **protección de la
   cronología**, que es un resultado igual de importante.

---

## Referencias

**Fundamentos**
- Hawking, S. W. & Ellis, G. F. R. (1973). *The Large Scale Structure of Space-Time*. Cambridge University Press.
- Visser, M. (1995). *Lorentzian Wormholes: From Einstein to Hawking*. AIP Press.
- Earman, J., Smeenk, C. & Wüthrich, C. (2009). Do the laws of physics forbid the operation of time machines? *Synthese* 169, 91.

**Soluciones con CTC**
- Gödel, K. (1949). *Rev. Mod. Phys.* 21, 447.
- Tipler, F. J. (1974). Rotating cylinders and the possibility of global causality violation. *Phys. Rev. D* 9, 2203.
- Morris, M. S. & Thorne, K. S. (1988). *Am. J. Phys.* 56, 395.
- Morris, M. S., Thorne, K. S. & Yurtsever, U. (1988). Wormholes, time machines, and the weak energy condition. *Phys. Rev. Lett.* 61, 1446.
- Gott, J. R. (1991). *Phys. Rev. Lett.* 66, 1126.
- Alcubierre, M. (1994). *Class. Quantum Grav.* 11, L73.
- Everett, A. E. (1996). Warp drive and causality. *Phys. Rev. D* 53, 7365.
- Ori, A. (2005). A class of time-machine solutions with a compact vacuum core. *Phys. Rev. Lett.* 95, 021101.

**Protección de la cronología y condiciones de energía**
- Tipler, F. J. (1977). Singularities and causality violation. *Ann. Phys.* 108, 1.
- Hawking, S. W. (1992). Chronology protection conjecture. *Phys. Rev. D* 46, 603.
- Kim, S.-W. & Thorne, K. S. (1991). Do vacuum fluctuations prevent the creation of closed timelike curves? *Phys. Rev. D* 43, 3929.
- Kay, B. S., Radzikowski, M. J. & Wald, R. M. (1997). *Commun. Math. Phys.* 183, 533.
- Li, L.-X. & Gott, J. R. (1998). Self-consistent vacuum for Misner space and the chronology protection conjecture. *Phys. Rev. Lett.* 80, 2980.
- Ford, L. H. & Roman, T. A. (1996). Quantum field theory constrains traversable wormhole geometries. *Phys. Rev. D* 53, 5496.
- Graham, N. & Olum, K. D. (2007). Achronal averaged null energy condition. *Phys. Rev. D* 76, 064001.
- Faulkner, T., Leigh, R. G., Parrikar, O. & Wang, H. (2016). Modular Hamiltonians for deformed half-spaces and the averaged null energy condition. *JHEP* 09, 038.
- Hartman, T., Kundu, S. & Tajdini, A. (2017). Averaged null energy condition from causality. *JHEP* 07, 066.

**Autoconsistencia y modelos cuánticos de CTC**
- Friedman, J. et al. (1990). Cauchy problem in spacetimes with closed timelike curves. *Phys. Rev. D* 42, 1915.
- Echeverria, F., Klinkhammer, G. & Thorne, K. S. (1991). *Phys. Rev. D* 44, 1077.
- Deutsch, D. (1991). Quantum mechanics near closed timelike lines. *Phys. Rev. D* 44, 3197.
- Aaronson, S. & Watrous, J. (2009). Closed timelike curves make quantum and classical computing equivalent. *Proc. R. Soc. A* 465, 631.
- Lloyd, S. et al. (2011). Closed timelike curves via postselection: theory and experimental test of consistency. *Phys. Rev. Lett.* 106, 040403.
- Ringbauer, M. et al. (2014). Experimental simulation of closed timelike curves. *Nat. Commun.* 5, 4145.

**Gravedad cuántica**
- Bombelli, L., Lee, J., Meyer, D. & Sorkin, R. D. (1987). Space-time as a causal set. *Phys. Rev. Lett.* 59, 521.
- Ambjørn, J., Jurkiewicz, J. & Loll, R. (2004). Emergence of a 4D world from causal quantum gravity. *Phys. Rev. Lett.* 93, 131301.
- Gauntlett, J. P. et al. (2003). All supersymmetric solutions of minimal supergravity in five dimensions. *Class. Quantum Grav.* 20, 4587.
- Boyda, E. K., Ganguli, S., Hořava, P. & Varadarajan, U. (2003). Holographic protection of chronology in universes of the Gödel type. *Phys. Rev. D* 67, 106003.
- Maldacena, J. & Susskind, L. (2013). Cool horizons for entangled black holes. *Fortsch. Phys.* 61, 781.
- Gao, P., Jafferis, D. L. & Wall, A. C. (2017). Traversable wormholes via a double trace deformation. *JHEP* 12, 151.
- Oppenheim, J. (2023). A postquantum theory of classical gravity? *Phys. Rev. X* 13, 041040.

**Experimentos**
- Hafele, J. C. & Keating, R. E. (1972). *Science* 177, 166.
- Colella, R., Overhauser, A. W. & Werner, S. A. (1975). Observation of gravitationally induced quantum interference. *Phys. Rev. Lett.* 34, 1472.
- Chou, C. W., Hume, D. B., Rosenband, T. & Wineland, D. J. (2010). Optical clocks and relativity. *Science* 329, 1630.
- Abdo, A. A. et al. (2009). A limit on the variation of the speed of light arising from quantum gravity effects. *Nature* 462, 331.
- Margalit, Y. et al. (2015). A self-interfering clock as a "which path" witness. *Science* 349, 1205.
- Bose, S. et al. (2017). Spin entanglement witness for quantum gravity. *Phys. Rev. Lett.* 119, 240401.
- Marletto, C. & Vedral, V. (2017). Gravitationally induced entanglement between two massive particles is sufficient evidence of quantum effects in gravity. *Phys. Rev. Lett.* 119, 240402.
- Zych, M., Costa, F., Pikovski, I. & Brukner, Č. (2019). Bell's theorem for temporal order. *Nat. Commun.* 10, 3772.
- Bothwell, T. et al. (2022). Resolving the gravitational redshift across a millimetre-scale atomic sample. *Nature* 602, 420.
