# Kisspeptina y vías de fosforilación río abajo de KISS1R

Estudio integrado: bibliometría, minería de literatura (coocurrencia de genes) e inferencia causal sobre transcriptomas públicos.

**Pregunta.** ¿Qué vías de fosforilación río abajo de KISS1R respaldan los datos transcriptómicos públicos, y cuáles de ellas no ha explorado la literatura?

| Hipótesis | Contenido | Módulo |
|---|---|---|
| H1 | La literatura forma comunidades oncológica y neuroendocrina que estudian quinasas distintas y se citan poco entre sí | Literatura |
| H2 | La inferencia causal sobre series con KISS1R perturbado recupera las quinasas canónicas (ERK1/2, p38, PI3K/AKT, GSK3β) | Datos |
| H3 | Hay quinasas con alto soporte en datos y bajo soporte en literatura; las quinasas oncológicas aparecen en datos neuroendocrinos más de lo esperado | Cruce |

## Estructura

```
kisspeptin-kinase/
├── R/
│   ├── 00_setup.R                 paquetes, rutas, semillas, versiones
│   ├── 01_search_openalex.R       corpus OpenAlex (sensibilidad) con ecuación v1.1
│   ├── 02_import_wos_scopus.R     importación y deduplicación de exportes WoS/Scopus
│   ├── 03_bibliometrics.R         co-citación, comunidades, conectividad, burst, ARI entre bases
│   ├── 04_flags_context.R         marcas S (quinasas), E (efectores), contexto oncológico, plantilla de anotación
│   ├── 05_pubtator_genes.R        anotación PubTator3 y mapeo a símbolo humano oficial
│   ├── 06_cooccurrence.R          coocurrencia por oración, NPMI, backbone, estratos, matriz gen × contexto
│   ├── 07_geo_inference.R         limma, decoupleR (CollecTRI, PROGENy), CARNIVAL + controles
│   └── 08_integration.R           cuadrantes, prueba de H2, prueba de transferencia, análisis temporal
├── data/
│   ├── raw/{wos,scopus,openalex,pubtator,geo,curated}   insumos sin modificar (no versionar exportes con licencia)
│   ├── processed/                 objetos intermedios .rds
│   └── annotation/                anotación manual (dirección del efecto, validación)
├── docs/
│   ├── search_equations_v1.1.txt  ecuaciones para WoS, Scopus y OpenAlex
│   ├── preregistration_osf.md     borrador del preregistro (plantilla de datos secundarios)
│   ├── decisions.md               decisiones abiertas y cerradas, con fecha
│   ├── estado_proyecto_*.md       estado y resultados a la fecha
│   ├── protocolo_*.md             datos y licencias, reproducibilidad, citación
│   └── organigrama_kisspeptina.html
└── output/{figures,tables,supplementary}
```

## Orden de ejecución

1. Cerrar las decisiones de `docs/decisions.md` y registrar el preregistro en OSF.
2. Colocar los exportes de WoS y Scopus en `data/raw/wos/` y `data/raw/scopus/`.
3. Colocar en `data/raw/curated/`: `kinome_reference.csv`, `cancer_gene_census.csv`, `cancermine_collated.tsv`.
4. Ejecutar los scripts en orden numérico. Cada script lee de `data/` y escribe en `data/processed/` u `output/`.
5. La rama de datos (`07`) no debe leer resultados de `03`–`06` antes de `08` (regla de independencia).

## Requisitos

- R ≥ 4.3, `renv` para fijar versiones (`renv::init()` en la primera sesión).
- En `~/.Renviron`: `OPENALEX_API_KEY` (clave gratuita de OpenAlex) y `OPENALEX_MAILTO` (correo de contacto). No se escriben en los scripts.
- Solver para CARNIVAL: HiGHS 1.15.1 en `../tools/` o en `HIGHS_PATH` (ver `docs/protocolo_reproducibilidad.md`).

## Reproducibilidad

- Semilla global: `20261014`.
- Fecha de extracción y versión de cada base se guardan junto a cada objeto `.rds`.
- Los exportes crudos de WoS/Scopus no se publican (licencia); se publican identificadores (DOI, UT, EID), ecuaciones y código.

## Protocolos

| Documento | Contenido |
|---|---|
| [`CONTRIBUTING.md`](CONTRIBUTING.md) | Ramas, mensajes de commit, etiquetas de versión y cambios que afectan al preregistro |
| [`docs/protocolo_datos.md`](docs/protocolo_datos.md) | Qué se versiona, licencias de cada fuente, cómo reconstruir los datos y secretos |
| [`docs/protocolo_reproducibilidad.md`](docs/protocolo_reproducibilidad.md) | Versiones, solvers, orden de ejecución, tiempos y cachés |
| [`docs/protocolo_citacion.md`](docs/protocolo_citacion.md) | OSF, Zenodo, DOI y cómo citar |

Antes del primer commit en un clon nuevo: `git config core.hooksPath .githooks`.

## Licencia

Código bajo [MIT](LICENSE). Documentos, tablas, figuras y anotaciones bajo [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). Los datos de terceros mantienen su propia licencia (`docs/protocolo_datos.md`). Para citar, ver [`CITATION.cff`](CITATION.cff).

## Autor

Gustavo Bruges · Universidad Autónoma de Bucaramanga (UNAB) · ORCID 0000-0001-8566-1868
