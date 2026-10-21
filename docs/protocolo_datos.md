# Protocolo de datos, licencias y secretos

**Regla general:** se versionan el código, los identificadores, las anotaciones propias y los resultados derivados. No se versionan los textos con licencia, los archivos de más de 50 MB ni los secretos. Cada insumo excluido se puede reconstruir con un script o una descarga documentada.

## 1. Qué se versiona

| Ruta | Contenido | ¿Se versiona? | Motivo |
|---|---|---|---|
| `R/` | Scripts de análisis | ✅ | Licencia MIT |
| `docs/*.md`, `docs/*.txt`, `docs/*.html` | Preregistro, decisiones, ecuaciones, protocolos | ✅ | CC BY 4.0 |
| `docs/*.xlsx` | Hojas de cálculo de trabajo | ❌ | Se revisan una por una antes de publicarlas |
| `docs/manuscrito/`, `*.docx` | Borradores del manuscrito | ❌ | Privados hasta el envío; se publican como *preprint*, no en el repositorio |
| `output/tables/`, `output/supplementary/`, `output/figures/` | Resultados derivados (conteos, puntajes, identificadores) | ✅ | CC BY 4.0; no contienen texto de terceros |
| `data/annotation/community_labels_manual.csv`, `offtopic_exclusions.csv` | Decisiones manuales (identificador, título, motivo) | ✅ | Los títulos son metadatos bibliográficos |
| `data/annotation/direction_annotation_template.csv`, `direction_annotator*.xlsx` | Plantillas con resúmenes de Scopus | ❌ | La licencia de Scopus no permite redistribuir resúmenes |
| `data/annotation/direction_annotation_public.csv` | Anotación consolidada sin título ni resumen (PMID, DOI, códigos, anotador) | ✅ cuando exista | Se genera con el script de consolidación |
| `data/raw/curated/kinome_reference.csv`, `pkinfam_*.txt` | Kinoma UniProt | ✅ | UniProt, CC BY 4.0 |
| `data/raw/{wos,scopus}/` | Exportes de WoS y Scopus | ❌ | Licencia de la suscripción UNAB |
| `data/raw/{openalex,pubtator,geo}/` | Descargas de OpenAlex, PubTator3 y GEO | ❌ | Se pueden volver a descargar; ocupan mucho |
| `data/raw/curated/gene2pubmed_*`, `gene_orthologs_*` | NCBI | ❌ | Más de 100 MB (límite de GitHub) |
| `data/raw/curated/omnipath_*` | Caché fechada de OmniPath | ❌ | Las licencias varían según el recurso de origen |
| `data/raw/curated/cancer_gene_census*` | COSMIC CGC | ❌ **nunca** | Licencia personal de COSMIC |
| `data/processed/` | Objetos `.rds` intermedios y caché de CARNIVAL | ❌ | Derivados y grandes; se archivan en Zenodo (sección 4) |

## 2. Licencia de cada fuente

Estas licencias están anotadas de memoria. **Verifica cada una en la página de la fuente** antes de publicar y anota aquí la fecha de verificación.

| Fuente | Licencia | Verificada |
|---|---|---|
| Scopus (Elsevier), Web of Science (Clarivate) | Uso según la suscripción; sin redistribución de registros | [ ] |
| OpenAlex | CC0 | [ ] |
| PubTator3, GEO, gene2pubmed, gene_orthologs (NCBI) | Dominio público (obra del gobierno de EE. UU.); hay que respetar las condiciones de los autores de cada serie GEO | [ ] |
| UniProt | CC BY 4.0 | [ ] |
| OmniPath / CollecTRI / PROGENy | Distintas según el recurso; CollecTRI y PROGENy tienen licencia propia | [ ] |
| COSMIC Cancer Gene Census | Licencia académica no redistribuible | [ ] |
| CancerMine | CC0 | [ ] |

## 3. Cómo reconstruir los datos excluidos

| Insumo | Cómo obtenerlo |
|---|---|
| WoS y Scopus | Ejecutar las ecuaciones de `docs/search_equations_v1.1.txt` en la fecha registrada, o volver a exportar a partir de los identificadores (EID, UT, DOI) de `output/supplementary/S1_*` |
| OpenAlex | `R/01_search_openalex.R` (requiere `OPENALEX_API_KEY` y `OPENALEX_MAILTO`) |
| PubTator3 | `R/05_pubtator_genes.R` a partir de `corpus_pmids.txt` |
| GEO | `R/07a_process_cel_GSE28383.R` y `R/07_geo_inference.R` descargan con `GEOquery` |
| NCBI gene2pubmed y gene_orthologs | Descarga FTP de NCBI Gene, en la fecha del nombre del archivo |
| OmniPath | `R/07_geo_inference.R` consulta la API web; para reproducir exactamente, usar la caché fechada archivada (sección 4) |
| COSMIC CGC | Cada usuario la descarga con su propia licencia y la deja en `data/raw/curated/` |

## 4. Archivo de datos derivados

Cada versión etiquetada (`v1.0-…`) se archiva en Zenodo enlazado a GitHub. Hay dos depósitos:

- **El repositorio:** código, documentos y `output/`. Zenodo lo toma automáticamente de cada *release* de GitHub.
- **Los datos derivados**, en un depósito aparte: `data/processed/` sin los objetos que contienen texto de Scopus o WoS (`M_scopus.rds`, `M_wos.rds`, `corpus_flagged.rds`, `pubtator_text.rds`), la caché de CARNIVAL y las cachés de OmniPath cuya licencia lo permita.

## 5. Secretos y datos personales

- Las claves y el correo de contacto de las API viven **sólo** en `~/.Renviron` (`OPENALEX_API_KEY`, `OPENALEX_MAILTO`). Los scripts los leen con `Sys.getenv()` y nunca los imprimen.
- `.Renviron` y `.env` están en `.gitignore`.
- El *hook* `.githooks/pre-commit` bloquea los *commits* que incluyan:
  - claves o correos escritos en el código;
  - archivos de más de 50 MB;
  - archivos con contenido bajo licencia.
- Para activarlo, una sola vez por clon: `git config core.hooksPath .githooks`.
- Si una clave llega a publicarse, **revócala de inmediato** en el servicio. Borrarla del historial de git no basta, porque GitHub y sus copias pueden haberla guardado ya.

## 6. Antes de publicar algo nuevo

- [ ] ¿Contiene títulos con resumen, texto completo o registros exportados de Scopus o WoS? → no se versiona.
- [ ] ¿Pesa más de 50 MB? → Zenodo, no git.
- [ ] ¿Viene de COSMIC? → nunca.
- [ ] ¿Es un borrador del manuscrito? → queda en `docs/manuscrito/`, fuera de git.
- [ ] ¿Tiene claves, correos o rutas personales? → al `.Renviron`.
- [ ] ¿La fuente tiene una licencia sin verificar en la sección 2? → verificarla primero.
