# Estado del proyecto · 2026-10-19

**Proyecto:** Kisspeptina y vías de fosforilación río abajo de KISS1R
**Responsable:** GB · Decisiones: `docs/decisions.md` · Detalle técnico: `docs/avances_2026-10-19.md`

---

## 1. Resumen

El análisis está **completo de punta a punta**: corrieron `00`–`08` y están todas las tablas de resultados. Faltan las tareas de anotación manual y las licencias, que pueden afinar la rama de literatura, y la redacción del manuscrito.

| Hipótesis | Resultado | Lectura |
|---|---|---|
| **H1** · convergencia entre comunidades | ✅ **Se sostiene** | Índice observado/esperado entre oncología y neuroendocrinología: 0,27 (2001–05) → 0,79 (2021–25) |
| **H2** · KISS1R recupera las quinasas canónicas | ❌ **No se sostiene** | El receptor ficticio HRH4 recupera tantas canónicas como KISS1R o más, en las dos series |
| **H3** · transferencia de quinasas oncológicas al hipotálamo | ⚪ **No evaluable** | La red hipotalámica estable sólo contiene MAPK14, que no lleva etiqueta oncológica ni neuroendocrina |

**Hipótesis priorizada:** **JNK1/2 (MAPK8, MAPK9)**. Son las únicas quinasas específicas: no salen con el receptor ficticio y desaparecen con siGPR54. Sólo tienen soporte en GSE25458 (endometrio humano).

---

## 2. Estado por script

| Script | Contenido | Estado |
|---|---|---|
| `00_setup.R` | Parámetros preregistrados, rutas, `.Renviron` | ✅ |
| `00_kinome_reference.R` | Kinoma UniProt (512 quinasas) | ✅ |
| `01_search_openalex.R` | Corpus OpenAlex (5.925 registros) | ✅ |
| `02_import_wos_scopus.R` | Scopus 4.213 y WoS 3.846, exclusiones fuera de tema | ✅ |
| `03_bibliometrics.R` | Comunidades e índice de conectividad (H1) | ✅ |
| `04_flags_context.R` | Marcas de fosforilación, plantillas de anotación | ✅ (falta la anotación manual) |
| `05_pubtator_genes.R` | Menciones de genes y ortología | ✅ |
| `06_cooccurrence.R` | Coocurrencia, `lit_kin` | ✅ (faltan CGC/CancerMine y dirección) |
| `07a_process_cel_GSE28383.R` | RMA de GSE28383 | ✅ |
| `07_geo_inference.R` | limma → decoupleR → CARNIVAL, controles, `data_kin` | ✅ |
| `08_integration.R` | Cuadrantes, Tabla 2, H2, H3, análisis temporal, figura 6 | ✅ |

---

## 3. Resultados principales

### Rama de literatura
- Corpus principal: Scopus 1996–2025. WoS y OpenAlex se usan como sensibilidad.
- Tres comunidades: oncología, neuroendocrina y molecular (receptor y farmacología).
- Las marcas de fosforilación se reparten de forma desigual: 21,5 % en oncología, 25 % en molecular y 7 % en neuroendocrino.
- 106 quinasas mencionadas. Ocho tienen literatura alta: MAP2K4, FGFR1, AMHR2, MAPK3, EPHB1, MAP2K7, MAPK1 y PRKCG.

### Rama de datos (exploratoria: sólo 2 series cumplen los criterios)
| Serie | Modelo | Señal |
|---|---|---|
| GSE25458 | Células Ishikawa, kp-10 × siRNA | Fuerte: 134 genes con FDR < 0,05; depende del receptor |
| GSE28383 | Hipotálamo de ratón, doble KO | Débil: 0 genes con FDR < 0,10; EGFR baja de forma consistente |

Configuración de CARNIVAL:
- solver HiGHS 1.15.1;
- poda de la red previa a 5 pasos;
- sólo se aceptan óptimos probados;
- 100 submuestreos al 80 %, umbral ≥ 0,7;
- receptor ficticio HRH4 con el mismo protocolo y semillas pareadas.

### Integración
- **Cuadrantes:** 1 confirmada (MAPK1), 7 hipótesis nuevas, 7 no reproducidas, 93 sin evidencia.
- **H2:**

| Serie | Canónicas con KISS1R | Canónicas con HRH4 |
|---|---|---|
| GSE25458 | 3 de 6 | 5 de 6 |
| GSE28383 | 1 de 6 | 2 de 6 |

  El signo esperado tampoco se cumple en GSE25458: AKT1, PIK3CA y MAPK14 salen bajas con el agonista.
- **Tabla 2:** MAPK9 y MAPK8 (específicas), luego PLK1, CDK2 y CSNK2A1 (no específicas).
- **Canónicas poco citadas en el corpus:** MAPK14 (la única en las dos series, y exclusiva de KISS1R en el hipotálamo) y GSK3B.
- **Análisis temporal:** la mediana del año de primera mención es 2001 para la confirmada (n = 1) y 2021 para las hipótesis nuevas (n = 5 con año conocido). Falta interpretarlo.

---

## 4. Manuscrito

| Sección | Estado | Ubicación |
|---|---|---|
| Introducción | Esqueleto listo; falta la redacción | [Claude Docs](https://claude.ai/code/artifact/428e50da-10a2-47d8-bd67-6b8f08b7f15c) |
| Métodos | Sin empezar; `docs/decisions.md` sirve de base | — |
| Resultados | Sin empezar; tablas en `output/tables/` | — |
| Discusión | Esqueleto listo; falta la redacción | [Claude Docs](https://claude.ai/code/artifact/27868ebc-788a-4666-b121-40e6b0579be6) |
| Figuras | Sólo la figura 6 (`fig6_quadrants.pdf`) | `output/figures/` |

---

## 5. Productos y reproducibilidad

- **Tablas clave** (`output/tables/`):
  - `kinase_quadrants.csv`
  - `table2_prioritized_hypotheses.csv`
  - `canonical_undercited.csv`
  - `h2_canonical_recovery.csv`
  - `h2_sign_check.csv`
  - `h3_transfer_test.csv`
  - `connectivity_by_period.csv`
  - `temporal_first_mention.csv`
  - `lit_kin.csv`
  - `data_kin.csv`
- **Caché de CARNIVAL:** `data/processed/carnival_cache/`, con un archivo por contraste más su versión `_fake_HRH4`. Si se vuelve a correr `07`, la red no se recalcula.
- **Solvers** (fuera del proyecto, en `../tools/`, con SHA-256 en `decisions.md`):
  - HiGHS 1.15.1, que se usa en el análisis;
  - CBC 2.10.12, que queda como respaldo.
- **Respaldo de la corrida anterior** (CBC con D = 4): `geo_inference_results_cbc_d4.rds` y `data_kin_cbc_d4.csv`.
- **Credenciales:** la clave de OpenAlex vive sólo en `.Renviron`, nunca en los scripts.

---

## 6. Pendientes

**A cargo de GB**
- [ ] Doble anotación de dirección (supresor o promotor) de 61 artículos, luego consolidación y kappa.
- [ ] Conseguir un segundo anotador.
- [ ] Descargar COSMIC CGC (requiere licencia) y revisar las licencias de CGC y CancerMine.
- [ ] Documentar la ecuación de búsqueda v1.2 (`KISSR*` → `kissr1–4`).
- [ ] Redactar la introducción y la discusión a partir de los esqueletos.

**Código y análisis (Claude, cuando se pida)**
- [ ] Script de consolidación y kappa de la anotación, cuando lleguen las dos planillas.
- [ ] Integrar CGC/CancerMine en `06` y volver a correr `06` → `08`.
- [ ] Opcional: estabilidad del control siGPR54, que hoy es una corrida única.
- [ ] Figuras 1–5 y tablas del manuscrito en formato final.
- [ ] Interpretar el análisis temporal.

**Riesgos que conviene vigilar**
- La anotación de dirección y CGC/CancerMine pueden reetiquetar comunidades de quinasas. Eso cambiaría `lit_kin` y los cuadrantes, aunque no H2 ni la especificidad, que salen sólo de los datos.
- Con dos series, `data_freq` sólo puede valer 0, 0,5 o 1. Cualquier conclusión del módulo de datos debe presentarse como exploratoria.
