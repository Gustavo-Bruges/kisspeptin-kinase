# Protocolo de reproducibilidad

**Objetivo:** que alguien con un clon del repositorio, las mismas versiones y los insumos de `docs/protocolo_datos.md` obtenga las mismas tablas de `output/`.

## 1. Entorno

| Componente | Versión usada (2026-10-19) |
|---|---|
| R | 4.6.1 (Windows, ucrt) |
| CARNIVAL | 2.22.0 |
| decoupleR | 2.17.0 |
| limma | 3.68.5 |
| oligo | 1.76.0 |
| progeny | 1.34.0 |
| OmnipathR | 4.0.0 (sólo para utilidades: las interacciones se bajan de la API web, ver `decisions.md`) |
| bibliometrix | 5.5.0 |
| igraph | 2.3.3 |
| openalexR | 3.1.2 |

**Pendiente:** `renv` todavía no está inicializado. Antes de la etiqueta `v1.0-analisis`:

```r
install.packages("renv")
renv::init(bare = TRUE)
renv::snapshot()          # writes renv.lock with every package version; commit it
```

Quien clone el repositorio después corre `renv::restore()`.

## 2. Solvers externos

Viven fuera del proyecto, en `../tools/`, y no se versionan. `R/07_geo_inference.R` los busca ahí o en la variable de entorno indicada.

| Solver | Versión | SHA-256 | Variable | Uso |
|---|---|---|---|---|
| HiGHS | 1.15.1 | `b6cb06d5…2808d7` (completo en `decisions.md`) | `HIGHS_PATH` | Análisis principal (`params$carnival_solver = "highs"`) |
| CBC | 2.10.12 | `6acf3e30…83b3a3` | `CBC_PATH` | Respaldo; óptimo sólo hasta D = 4 |

Para verificar el binario: `sha256sum ../tools/highs-1.15.1/bin/highs.exe` y comparar con `decisions.md`.

## 3. Orden de ejecución y tiempos

| Paso | Script | Tiempo | Requiere |
|---|---|---|---|
| 1 | `00_kinome_reference.R` | por medir | Internet (UniProt) |
| 2 | `01_search_openalex.R` | por medir | `OPENALEX_API_KEY`, `OPENALEX_MAILTO` |
| 3 | `02_import_wos_scopus.R` | por medir | Exportes de WoS y Scopus en `data/raw/` |
| 4 | `03_bibliometrics.R` | por medir | — |
| 5 | `04_flags_context.R` | por medir | — |
| 6 | `05_pubtator_genes.R` | por medir | Internet (PubTator3), archivos NCBI |
| 7 | `06_cooccurrence.R` | por medir | CGC y CancerMine (opcional por ahora) |
| 8 | `07a_process_cel_GSE28383.R` | por medir | CEL de GEO |
| 9 | `07_geo_inference.R` | **~20–30 h** sin caché; minutos con caché | HiGHS, 8 hilos |
| 10 | `08_integration.R` | por medir | — |

Sólo el tiempo de `07` está medido. Viene de la corrida de 2026-10-18/19: entre 168 y 400 s por corrida de CARNIVAL, con 100 submuestreos por contraste para KISS1R y otros 100 para HRH4. Para correr las series en paralelo, en dos terminales:

```sh
KISS_GSE=GSE25458 Rscript R/07_geo_inference.R
KISS_GSE=GSE28383 Rscript R/07_geo_inference.R
Rscript R/07_geo_inference.R        # without KISS_GSE: assembles data_kin from the cache
```

## 4. Fuentes de variación y cómo se controlan

| Fuente | Control |
|---|---|
| Semillas | `params$seed = 20261014`; los submuestreos usan `seed + i`, iguales para KISS1R y HRH4 |
| Solver MIP | HiGHS con `random_seed` fijo y `mip_rel_gap` de CARNIVAL; sólo se aceptan óptimos probados |
| Hilos | `params$carnival_threads = 8`. Con otro número de hilos el óptimo es el mismo, pero si hay empates puede cambiar la solución elegida: reporta los hilos usados |
| Versiones de bases | Cada `.rds` guarda la fecha y la versión de la consulta (`save_versioned()`); las cachés llevan la fecha en el nombre |
| OmniPath cambia con el tiempo | Caché fechada `omnipath_*_2026-10-18.tsv`; archivarla en Zenodo para la reproducción exacta |
| Caché de CARNIVAL | `data/processed/carnival_cache/<gse>_<contraste>_<solver>_d<D>_n<N>.rds`. Si cambias `solver`, `D` o `N`, el nombre cambia y se recalcula; cualquier otro cambio en `07` exige borrar la caché a mano |

## 5. Verificación de una reproducción

- [ ] `renv::status()` sin diferencias.
- [ ] El SHA-256 de HiGHS coincide.
- [ ] `output/tables/h2_canonical_recovery.csv`, `table2_prioritized_hypotheses.csv` y `data_kin.csv` son idénticos a los de la etiqueta (`git diff --stat v1.0-analisis -- output/`).
- [ ] En `07` resolvieron al menos el 90 % de las corridas de estabilidad (`n_solved`); si no, se emite un *warning*.
