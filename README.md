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
| [`simulaciones/ctc_cuantico.py`](simulaciones/ctc_cuantico.py) | Modelos de Deutsch (D-CTC) y de post-selección (P-CTC): paradoja del abuelo, bucle de información y no linealidad. |
| [`simulaciones/escalas_fisicas.py`](simulaciones/escalas_fisicas.py) | Órdenes de magnitud: dilatación temporal (GPS, relojes a 33 cm), energía de Casimir frente a la que necesita un agujero de gusano, escala de Planck. |

## Ejecutar

```bash
pip install numpy
python3 simulaciones/ctc_cuantico.py
python3 simulaciones/escalas_fisicas.py
```
