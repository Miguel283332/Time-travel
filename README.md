# Time-travel

Planteamiento de un programa de investigación para decidir si los viajes al
pasado son físicamente posibles, usando la gravedad cuántica como árbitro.

**Idea central:** la relatividad general permite curvas temporales cerradas;
la física semiclásica deja de estar definida justo en el horizonte donde se
formarían; por tanto, el veredicto corresponde a la gravedad cuántica.
Demostrar la gravedad cuántica es necesario, pero no suficiente.

## Contenido

| Archivo | Qué es |
|---|---|
| [`docs/planteamiento.md`](docs/planteamiento.md) | Documento principal: definiciones, base clásica, teoremas, gravedad cuántica, estructura lógica, hipótesis falsables, hoja de ruta y referencias. |
| [`docs/gravedad-cuantica-bucles.md`](docs/gravedad-cuantica-bucles.md) | Programa para confirmar la gravedad cuántica de bucles: niveles de demostración, estado actual, resultados reproducidos y ruta hacia una predicción observable (fondo cósmico de microondas). |
| [`docs/resultados-cmb.md`](docs/resultados-cmb.md) | **Primer contraste con datos reales:** espectro primordial con rebote cuántico → fondo cósmico de microondas → Planck 2018. Resultado, validaciones y siguientes pasos. |
| [`simulaciones/lqc_perturbaciones.py`](simulaciones/lqc_perturbaciones.py) | Perturbaciones a través del rebote (métrica vestida, Starobinsky); espectro promediado con coeficientes de Bogoliubov o evolución exacta (`--exacto`). |
| [`simulaciones/cmb_lqc.py`](simulaciones/cmb_lqc.py) | CAMB + verosimilitud de bajo ℓ de Planck 2018: barrido en la escala del rebote, cota y predicción para la polarización. |
| [`simulaciones/comparar_estados.py`](simulaciones/comparar_estados.py) | Compara el efecto del rebote y la cota de Planck para distintos estados iniciales. |
| [`simulaciones/lqg_espectros.py`](simulaciones/lqg_espectros.py) | Espectro de área de LQG, densidad del espectro y parámetro de Barbero-Immirzi fijado por la entropía de agujeros negros. |
| [`simulaciones/lqc_rebote.py`](simulaciones/lqc_rebote.py) | Rebote de la cosmología cuántica de bucles frente a la singularidad de la relatividad general. |
| [`simulaciones/ctc_cuantico.py`](simulaciones/ctc_cuantico.py) | Modelos de Deutsch (D-CTC) y de post-selección (P-CTC): paradoja del abuelo, bucle de información y no linealidad. |
| [`simulaciones/escalas_fisicas.py`](simulaciones/escalas_fisicas.py) | Órdenes de magnitud: dilatación temporal (GPS, relojes a 33 cm), energía de Casimir frente a la que necesita un agujero de gusano, escala de Planck. |

## Ejecutar

```bash
pip install -r requirements.txt
python3 simulaciones/ctc_cuantico.py
python3 simulaciones/escalas_fisicas.py
python3 simulaciones/lqg_espectros.py
python3 simulaciones/lqc_rebote.py
python3 simulaciones/lqc_perturbaciones.py   # ~30 s con 4 núcleos
python3 simulaciones/cmb_lqc.py              # ~2 min; descarga los datos de Planck la primera vez
```
