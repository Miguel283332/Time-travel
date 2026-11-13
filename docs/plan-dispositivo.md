# Plan: un dispositivo de laboratorio para poner a prueba la gravedad cuántica

> **Idea central.** No existe forma conocida de "fabricar gravedad" sin masa,
> y la gravedad cuántica de bucles (LQG) no la proporciona. Lo que sí se puede
> construir es un dispositivo en el que **una masa en superposición cuántica
> genera un campo gravitatorio** y se comprueba si ese campo también está en
> superposición. Ese es el experimento que puede decidir si la gravedad es
> cuántica, y es el único "generador de gravedad" que dice algo sobre LQG.

Los números de este documento salen de
[`simulaciones/dispositivo_bmv.py`](../simulaciones/dispositivo_bmv.py).

---

## 0. De dónde partimos

| Parte del proyecto | Qué aporta a este plan |
|---|---|
| [`planteamiento.md`](planteamiento.md) | La pregunta de fondo (viajes al pasado) depende de la gravedad cuántica. |
| [`gravedad-cuantica-bucles.md`](gravedad-cuantica-bucles.md) | Los cinco niveles para confirmar LQG (N0–N4). El nivel N0, "¿es cuántica la gravedad?", sigue sin respuesta experimental. |
| [`resultados-cmb.md`](resultados-cmb.md) | La vía cosmológica (N3): hecha en parte; con los estados probados, Planck solo acota el modelo. |
| **Este plan** | La vía de laboratorio: ataca **N0** con un dispositivo, y explora si puede llegar a N3. |

---

## 1. Qué puede ser y qué no un "dispositivo de gravedad artificial"

| Significado | ¿Existe? | Tamaño del efecto | ¿Dice algo de LQG? |
|---|---|---|---|
| **a) Rotación** (nave o centrifugadora) | Sí, física de 1900 | 1 g con radio 10 m a 9,5 rpm | No. Es inercia, no gravedad. |
| **b) Campo de una masa fuente** controlada | Sí | 90 mg de oro a 2,5 mm: $10^{-9}$ m/s² (medido en 2021) | **Sí, si la masa fuente está en superposición cuántica.** |
| **c) Arrastre de sistemas de referencia** con un rotor | En principio | Rotor de 1 t a 1000 rad/s: $7\times10^{-22}$ rad/s, ocho millones de veces menos que el efecto de la Tierra que midió Gravity Probe B | No de forma práctica. |
| **d) Generador de gravedad sin masa, antigravedad, pantallas** | **No** | — | Ninguna teoría aceptada lo predice, LQG incluida. Las afirmaciones experimentales (Podkletnov, 1992) no se han reproducido. |

Dos razones por las que LQG no da un generador del tipo (d):

1. LQG modifica la gravedad a la escala de Planck. A un micrómetro, la
   corrección relativa por la discretización del área es
   $(\ell_P/L)^2 \sim 10^{-58}$.
2. La repulsión efectiva que aparece en el rebote cosmológico
   ([`lqc_rebote.py`](../simulaciones/lqc_rebote.py)) solo actúa a densidades
   $\rho\sim\rho_c\approx 0{,}41\,\rho_{\rm Planck}\sim10^{96}$ kg/m³.

Por tanto el plan se construye sobre **(b)**.

---

## 2. Por qué esto sí es posible: la masa de Planck

| Escala de Planck | Valor | ¿Accesible? |
|---|---|---|
| Longitud | $1{,}6\times10^{-35}$ m | No: 20 órdenes por debajo de cualquier acelerador. |
| Tiempo | $5{,}4\times10^{-44}$ s | No directamente. |
| **Masa** | $2{,}2\times10^{-8}$ kg = **22 µg** | **Sí: es un grano de polvo.** |

La gravedad cuántica se manifiesta cuando un objeto con masa cercana a la de
Planck se comporta de forma cuántica. Eso no requiere energías enormes:
requiere **control cuántico de objetos cada vez más pesados**. Es una carrera
tecnológica en marcha.

---

## 3. El dispositivo

### 3.1 Principio (Bose et al. 2017; Marletto y Vedral 2017)

1. Dos micropartículas de masa $m$ levitadas en vacío, a distancia $d$.
2. Cada una se pone en superposición de dos posiciones separadas $\Delta x$.
3. Se dejan caer libremente un tiempo $\tau$. Solo interactúan por gravedad.
4. Cada combinación de posiciones acumula una fase distinta,
   $\phi_{ij} = Gm^2\tau/(\hbar\, r_{ij})$.
5. Se recombinan y se mide si las dos partículas han quedado **entrelazadas**.

**Por qué es decisivo:** un teorema de información cuántica dice que un
mediador clásico no puede crear entrelazamiento. Si aparece entrelazamiento y
el único canal es la gravedad, el campo gravitatorio no es clásico.
Christodoulou y Rovelli (2019) lo interpretan, desde LQG, como una
**superposición de geometrías del espacio-tiempo**.

### 3.2 Parámetros de referencia y lo que exigen

| Magnitud | Valor | Comentario |
|---|---|---|
| Masa de cada partícula | $10^{-14}$ kg | Microdiamante de 0,88 µm de radio. |
| Separación de la superposición $\Delta x$ | 250 µm | |
| Distancia entre centros $d$ | 450 µm | Mínima entre ramas: 200 µm. |
| Tiempo de interacción $\tau$ | 2,5 s | En caída libre. |
| Fase de entrelazamiento | 0,31 rad | Concurrencia 0,16: pequeña pero medible con estadística. |
| Fuerza de Casimir–Polder / gravedad | 0,02 | Por eso no pueden acercarse más de ~200 µm. |
| Vacío | $< 10^{-15}$ Pa a 0,15 K | Menos de 0,1 colisiones con gas durante la medida. |
| Dependencia con la masa | señal ∝ $m^2$ | Con $10^{-15}$ kg la fase cae a 0,003 rad. |

### 3.3 Subsistemas

| Subsistema | Función | Tecnología candidata |
|---|---|---|
| Fuente y carga de partículas | Microdiamantes con un espín (centro NV) | Diamante de alta pureza isotópica |
| Levitación | Sostener sin contacto ni luz que caliente | Trampa magnética (diamagnética) o de Paul |
| Enfriamiento | Llevar el movimiento cerca del estado fundamental | Realimentación óptica o eléctrica |
| Creación de la superposición | Separar las dos ramas 250 µm | Interferometría de Stern–Gerlach con gradiente magnético |
| Criostato y vacío extremo | Suprimir decoherencia | Criostato de dilución, bombeo criogénico |
| Apantallamiento | Eliminar fuerzas eléctricas y magnéticas entre partículas | Lámina conductora entre ambas; control de carga |
| Lectura | Medir correlaciones de espín | Resonancia magnética detectada ópticamente |
| Aislamiento | Vibraciones y gradientes de gravedad | Plataforma en caída libre o aislamiento activo |

---

## 4. Qué demostraría cada resultado

| Resultado | Conclusión | Consecuencia para LQG |
|---|---|---|
| **Hay entrelazamiento** | La gravedad no es un campo clásico. Quedan descartadas la gravedad semiclásica y las teorías clásico-estocásticas (Oppenheim, Diósi–Penrose). | **Supera N0.** Pero es compatible por igual con cuerdas y con cualquier teoría cuántica de la gravedad: a estas energías todas predicen lo mismo. **No confirma LQG.** |
| **No hay entrelazamiento** (con sensibilidad suficiente) | La gravedad es clásica o la mecánica cuántica falla para masas grandes. | **Refuta LQG** (y también cuerdas). |
| **Escalones en la fase** (sección 5) | El tiempo propio sería discreto a escala de Planck. | Sería una señal a favor de teorías con geometría discreta como LQG. **Es el único resultado que apuntaría a LQG en concreto.** |

El dispositivo puede **refutar** LQG o **superar su premisa**, pero no
confirmarla frente a sus rivales salvo en el régimen de la sección 5.

---

## 5. La extensión específica de LQG: tiempo discreto

La fase de entrelazamiento es $\phi = mc^2\,\delta\tau/\hbar$, donde $\delta\tau$
es la diferencia de tiempo propio entre ramas por dilatación gravitatoria. Si
$\delta\tau$ solo pudiera tomar múltiplos del tiempo de Planck, la fase crecería
a saltos de tamaño

$$\Delta\phi = \frac{mc^2\,t_P}{\hbar} = \frac{m}{m_P}.$$

| Masa en superposición | Tamaño del escalón |
|---|---|
| $10^{-14}$ kg (dispositivo de la sección 3) | $5\times10^{-7}$ rad |
| $10^{-11}$ kg | $5\times10^{-4}$ rad |
| $2{,}2\times10^{-8}$ kg (masa de Planck) | 1 rad |

Propuesto por Christodoulou, Di Biagio y Martin-Dussaud (2022). Dos
advertencias: que el tiempo propio sea discreto es una **hipótesis inspirada
en LQG**, no una predicción firme; y exige masas seis órdenes de magnitud
mayores que el experimento base.

---

## 6. Etapas

Los plazos y costes son **estimaciones mías de orden de magnitud**, no
presupuestos.

| Etapa | Objetivo | Estado del arte | Qué falta | Escala |
|---|---|---|---|---|
| **E0. Diseño y simulación** | Modelo completo del dispositivo: señal, ruido, decoherencia, tolerancias | Artículos de propuesta | Simulación propia con presupuesto de errores | Meses; solo cálculo. **Lo podemos hacer aquí.** |
| **E1. Gravedad de una masa pequeña** | Medir el campo de una fuente de miligramos | Hecho: 90 mg (Westphal et al. 2021); partícula levitada de 0,43 mg (Fuchs et al. 2024) | Bajar a microgramos | 1–3 años; laboratorio universitario |
| **E2. Levitación y enfriamiento** | Micropartícula de $10^{-14}$ kg cerca del estado fundamental | Hecho con $10^{-18}$ kg (Delić et al. 2020) | 4 órdenes en masa | 3–5 años; laboratorio especializado |
| **E3. Superposición de una masa grande** | Una partícula de $10^{-14}$ kg en dos sitios separados 250 µm | Récord: $4\times10^{-23}$ kg, cientos de nm (Fein et al. 2019) | **8,4 órdenes en masa y ~3 en distancia. Es el cuello de botella.** | 10–20 años |
| **E4. Entrelazamiento por gravedad** | El experimento de la sección 3 | No realizado | E2 + E3 + apantallamiento | Tras E3; consorcio internacional |
| **E5. Tiempo discreto** | Sección 5 | Solo propuesta teórica | 6 órdenes más de masa | Muy largo plazo |

Grupos que trabajan en esta línea: UCL (Bose), Viena (Aspelmeyer),
Southampton (Ulbricht), Warwick (Morley), Ben-Gurión (Folman), Groningen
(Mazumdar). Un proyecto así solo es viable en colaboración con alguno.

---

## 7. Lo que podemos hacer ya en este repositorio (E0)

1. **Simulador del experimento:** evolución de las dos partículas con
   decoherencia (gas, radiación térmica, ruido de la trampa) y cálculo del
   testigo de entrelazamiento medible.
2. **Presupuesto de errores:** cuántas repeticiones hacen falta para detectar
   una concurrencia de 0,16 con cada fuente de ruido.
3. **Optimización del diseño:** barrer masa, separación, tiempo y geometría
   (colineal frente a paralela) para maximizar la señal con el vacío y la
   temperatura alcanzables.
4. **Fuerzas parásitas:** Casimir–Polder con y sin lámina de apantallamiento,
   dipolos eléctricos y magnéticos residuales.
5. **Discriminación de modelos:** qué predicen para el mismo dispositivo la
   gravedad cuántica, la semiclásica y la clásico-estocástica.
6. **Tiempo discreto:** qué precisión de fase y qué masa permitirían ver los
   escalones de la sección 5.

El resultado de E0 sería un documento de diseño con números propios,
presentable a un grupo experimental.

---

## 8. Riesgos y advertencias

- **El plan no produce gravedad artificial utilizable.** El campo que genera
  el dispositivo es de $\sim10^{-17}$ m/s² a 200 µm; su valor es lo que revela
  sobre la naturaleza de la gravedad.
- **Un resultado positivo no confirma LQG** frente a otras teorías cuánticas
  (sección 4). La vía con números propios de LQG sigue siendo la cosmológica.
- **E3 puede no ser alcanzable.** Si existe un mecanismo que destruye las
  superposiciones de masas grandes (modelos de colapso), el experimento lo
  descubriría, lo que también sería un resultado de primer orden.
- **Hay debate teórico** sobre si el entrelazamiento por gravedad demuestra
  la cuantización del campo o solo su carácter no clásico.
- Las referencias están citadas de memoria; conviene verificarlas.

---

## Referencias

- Bose, S. et al. (2017). Spin entanglement witness for quantum gravity. *Phys. Rev. Lett.* 119, 240401.
- Marletto, C. & Vedral, V. (2017). Gravitationally induced entanglement between two massive particles is sufficient evidence of quantum effects in gravity. *Phys. Rev. Lett.* 119, 240402.
- Christodoulou, M. & Rovelli, C. (2019). On the possibility of laboratory evidence for quantum superposition of geometries. *Phys. Lett. B* 792, 64.
- Christodoulou, M., Di Biagio, A. & Martin-Dussaud, P. (2022). An experiment to test the discreteness of time. *Quantum* 6, 826.
- Westphal, T., Hepach, H., Pfaff, J. & Aspelmeyer, M. (2021). Measurement of gravitational coupling between millimetre-sized masses. *Nature* 591, 225.
- Fuchs, T. M. et al. (2024). Measuring gravity with milligram levitated masses. *Sci. Adv.* 10, eadk2949.
- Delić, U. et al. (2020). Cooling of a levitated nanoparticle to the motional quantum ground state. *Science* 367, 892.
- Fein, Y. Y. et al. (2019). Quantum superposition of molecules beyond 25 kDa. *Nat. Phys.* 15, 1242.
- Margalit, Y. et al. (2021). Realization of a complete Stern-Gerlach interferometer: toward a test of quantum gravity. *Sci. Adv.* 7, eabg2879.
- Romero-Isart, O. (2011). Quantum superposition of massive objects and collapse models. *Phys. Rev. A* 84, 052121.
- Everitt, C. W. F. et al. (2011). Gravity Probe B: final results. *Phys. Rev. Lett.* 106, 221101.
- Hathaway, G., Cleveland, B. & Bao, Y. (2003). Gravity modification experiment using a rotating superconducting disk and radio frequency fields. *Physica C* 385, 488.
- Oppenheim, J. (2023). A postquantum theory of classical gravity? *Phys. Rev. X* 13, 041040.
