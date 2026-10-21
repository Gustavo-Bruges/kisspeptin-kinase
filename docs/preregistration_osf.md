# Registro retrospectivo · Plantilla OSF *Secondary Data Preregistration*

**Título:** Literature and transcriptomic evidence on kinase signaling downstream of KISS1R: an integrated bibliometric, gene co-occurrence and causal network analysis

**Autor:** Gustavo Bruges (UNAB) · ORCID 0000-0001-8566-1868
**Fecha de registro:** ____ · **DOI OSF:** ____ · **Embargo:** ninguno (registro público)

> **Estado del registro.** El plan de análisis (secciones 1–9) se fijó el 2026-10-14. El análisis exploratorio completo se ejecutó entre el 2026-10-18 y el 2026-10-21, **antes de este registro**. Está archivado en la etiqueta git `v1.0-analisis` (https://github.com/Gustavo-Bruges/kisspeptin-kinase/releases/tag/v1.0-analisis). Este registro **no** debe leerse como un preregistro confirmatorio de H1–H3. Su propósito es fijar públicamente el plan original, declarar cada decisión tomada después de acceder a los datos (sección 10) y preregistrar de verdad las pruebas confirmatorias que aún no se han hecho (sección 11).

---

## 1. Pregunta e hipótesis

**Pregunta.** ¿Qué vías de fosforilación río abajo de KISS1R respaldan los datos transcriptómicos públicos, y cuáles de ellas no ha explorado la literatura?

- **H1.** La literatura sobre kisspeptina forma al menos dos comunidades (oncológica y neuroendocrina) que estudian conjuntos distintos de quinasas y presentan baja conectividad entre sí.
- **H2.** La inferencia causal sobre series con KISS1R perturbado recupera ERK1/2, p38, PI3K/AKT y GSK3β con mayor frecuencia que una perturbación ficticia.
- **H3.** Existen quinasas con alto soporte en datos y bajo soporte en literatura; las quinasas descritas en la comunidad oncológica aparecen en la red hipotalámica más de lo esperado por azar.

## 2. Datos

### 2.1 Corpus bibliográfico
- Bases: Scopus (principal), Web of Science y OpenAlex (sensibilidad).
- Ecuaciones: `docs/search_equations_v1.1.txt`. Periodo 1996–2025; artículos y revisiones.
- Fecha de extracción: 2026-10-18 (Scopus, WoS y OpenAlex).

### 2.2 Datos transcriptómicos
- Series: GSE25458 (agonista, signo +1) y GSE28383 (knockout, signo −1), más las adicionales que cumplan los criterios: perturbación directa de KISS1/KISS1R, controles en la misma serie, n ≥ 3 por grupo. No se encontró ninguna adicional (sección 10).

### 2.3 Recursos de conocimiento
- PubTator3 (anotación de genes), OmniPath (red previa con signo, caché del 2026-10-18), CollecTRI y PROGENy.
- Kinoma de referencia: UniProt `pkinfam` (descarga del 2026-10-18).
- COSMIC CGC y CancerMine: previstos como capa de anotación. **Todavía no integrados** a la fecha del registro.

## 3. Conocimiento previo de los datos (declaración obligatoria)

**A la fecha de este registro, el autor ya descargó y analizó todos los datos descritos en la sección 2.** En detalle:

- El 2026-10-14, antes de descargar ningún dato, quedaron fijados la ecuación de búsqueda v1.1, las hipótesis y los umbrales de la sección 5 (`params` en `R/00_setup.R`). Las fechas provienen del registro de decisiones (`docs/decisions.md`). El proyecto no estuvo bajo control de versiones hasta el 2026-10-20, así que esas fechas no se pueden verificar de forma independiente.
- Antes del 2026-10-14, el autor había leído las descripciones públicas de GSE25458 y GSE28383 y había ejecutado una consulta preliminar en OpenAlex (~5.000–6.000 registros), sin analizar su contenido. También había consultado síntesis de literatura sobre quinasas en cáncer asociadas a KISS1/KISS1R, que motivaron la ampliación del bloque K a la versión 1.1.
- Entre el 2026-10-18 y el 2026-10-21 se extrajeron los corpus, se descargaron y procesaron las series GEO y se ejecutaron todos los análisis (scripts `00`–`08`).
- **Los resultados de H1–H3 eran conocidos al registrar:**
  - H1 se sostiene;
  - H2 no se sostiene: el receptor ficticio recupera tantas quinasas canónicas como KISS1R o más;
  - H3 no es evaluable: la red hipotalámica estable contiene una sola quinasa.
  - La única hipótesis específica que surgió es JNK1/2 (MAPK8, MAPK9).
- Todas las decisiones tomadas después de acceder a los datos están en la sección 10, con su fecha y su motivo.

## 4. Variables y procedimientos

- **Comunidades:** acoplamiento bibliográfico + Louvain; etiquetado semiautomático por diccionario y confirmación manual.
- **Conectividad:** peso de acoplamiento entre comunidades / peso total, por quinquenio.
- **Subcorpus S:** marca local con el bloque K v1.1 sobre título y resumen.
- **Contexto oncológico:** diccionario de tipos y subtipos de cáncer; dirección del efecto por anotación manual doble (kappa de Cohen).
- **Genes:** PubTator3, símbolo humano oficial (org.Hs.eg.db + ortólogos).
- **Coocurrencia:** oración como unidad; NPMI; hipergeométrica con FDR; backbone SDSM.
- **lit_score:** log2 del enriquecimiento de cada quinasa en el corpus frente a su frecuencia en PubMed.
- **Soporte en datos:** fracción de series en cuya red causal (CARNIVAL sobre OmniPath sin aristas circulares, nodos estables en ≥ 70 % de 100 submuestreos al 80 %) aparece la quinasa.

## 5. Umbrales fijados (2026-10-14)

| Parámetro | Valor |
|---|---|
| Soporte mínimo por par de genes | 3 |
| FDR | 0,05 |
| Literatura alta | lit_score > percentil 75 |
| Datos altos | presente en ≥ 50 % de las series |
| Estabilidad CARNIVAL | ≥ 70 % de 100 corridas al 80 % |
| Permutaciones (transferencia) | 10.000 |
| Periodos | 1996–2008, 2009–2016, 2017–2025 |
| Semilla | 20261014 |

## 6. Análisis

1. **H1:** comparación de la proporción S/C y de los perfiles de quinasas entre comunidades; tendencia del índice de conectividad; concordancia entre bases (ARI).
2. **H2:** número de quinasas canónicas recuperadas en la red real frente a la red con receptor ficticio, por serie.
3. **H3:** cuadrantes literatura × datos; prueba de transferencia por permutación de etiquetas de comunidad.
4. **Temporal:** año de primera mención de las quinasas confirmadas frente a las nuevas.

## 7. Sensibilidad y robustez

- Unidad de coocurrencia = resumen; exclusión de revisiones.
- Bases alternativas (ARI).
- Red previa con y sin aristas circulares.
- Bloque K v1.0 frente a v1.1 (sesgo de origen oncológico del bloque).

## 8. Criterios de decisión y respaldo

- Si hay menos de 3 series perturbadas elegibles en la semana 4, el módulo de datos se presenta como exploratorio o se separa en un segundo artículo. **Se aplicó:** sólo 2 series fueron elegibles.
- Si no hay segundo anotador, la vista supresor/promotor se restringe a una muestra estratificada o se elimina.

## 9. Validación

- Precisión de PubTator3 y precisión relacional en 150 resúmenes aleatorios de S, con IC 95 %.
- Precisión de la ecuación en 200 registros aleatorios por base; sensibilidad frente a un conjunto de referencia de 25–30 artículos definido antes de la búsqueda.

---

## 10. Decisiones tomadas después de acceder a los datos

El detalle completo, con fecha y motivo, está en `docs/decisions.md`. La columna "¿Con resultados a la vista?" indica si la decisión se tomó conociendo el resultado de la hipótesis afectada.

| Fecha | Decisión | ¿Con resultados a la vista? |
|---|---|---|
| 2026-10-18 | Scopus como base principal: la suscripción WoS de UNAB empieza en 2006 | No |
| 2026-10-18 | Exclusión auditable de 22 registros de Scopus y 14 de WoS fuera de tema | No |
| 2026-10-18 | Tercera comunidad *molecular*, obtenida con un segundo nivel de Louvain dentro de oncología | Sí, para H1 (la estructura de comunidades ya se conocía) |
| 2026-10-18 | "Literatura alta" exige además ≥ 3 artículos (`lit_min_docs`) | Sí, para los cuadrantes |
| 2026-10-18 | Comunidad de cada quinasa = la de mayor lift | Sí, para H3 |
| 2026-10-18 | Sólo 2 series GEO elegibles; se aplica el criterio de respaldo (sección 8) | No |
| 2026-10-18 | GSE28383 normalizado desde los CEL; lote como covariable; sensibilidad sin el lote 2 | No |
| 2026-10-18 | Solver HiGHS, poda de la red previa a ≤ 5 pasos, sólo óptimos probados | Sí: se eligió D = 5 tras ver que con D = 4 ninguna canónica podía entrar en la red |
| 2026-10-18 | Receptor ficticio = HRH4, elegido con una regla fijada antes de correrlo | No, para H2 |
| 2026-10-19 | Receptor ficticio con el mismo protocolo de estabilidad que KISS1R (semillas pareadas) | Sí: se agregó tras ver la corrida única del ficticio, y cambió la clasificación de MAPK14 |
| 2026-10-19 | Regla de la Tabla 2: se excluyen las canónicas y se ordena por especificidad | Sí |
| 2026-10-19 | H3 se reporta como no evaluable, sin relajar el criterio de intersección | Sí |

## 11. Plan confirmatorio (preregistro de análisis aún no realizados)

Esta sección sí es un preregistro: ninguno de estos análisis se había hecho a la fecha del registro. **El pipeline queda congelado en `v1.0-analisis`:** mismo código y mismos `params`, sin cambios.

### C1 · Réplica computacional en series nuevas
- **Datos:** series de GEO o ArrayExpress que cumplan los criterios de la sección 2.2 y que no se hayan analizado en `v1.0-analisis`. La búsqueda se repite con las mismas consultas, entre la fecha de registro y ____.
- **Predicción:** en cada serie nueva con perturbación de KISS1R, MAPK8 y/o MAPK9 son nodos estables (≥ 70 %) con KISS1R y no lo son con HRH4, bajo el mismo protocolo pareado.
- **Criterio de éxito:** la predicción se cumple en al menos la mitad de las series nuevas elegibles. Con ninguna serie nueva, C1 se reporta como no realizado.
- **Para H2:** se reporta el mismo conteo de canónicas, real frente a ficticio. Una sola serie nueva no basta para revertir la conclusión de `v1.0-analisis`.

### C2 · Validación experimental de JNK1/2 (si se realiza)
- **Modelo:** ____ (por ejemplo, células Ishikawa u otra línea que exprese KISS1R).
- **Diseño:** siRNA control o siGPR54 × vehículo o kp-10; tiempos ____; n = ____ réplicas biológicas independientes por condición.
- **Medida primaria:** razón fosfo-JNK/JNK total (western blot o ELISA, ____).
- **Predicción:** kp-10 aumenta la razón fosfo-JNK/JNK frente al vehículo con siCtrl, y el efecto desaparece con siGPR54.
- **Análisis:** ANOVA de dos vías (siRNA × tratamiento). La prueba primaria es la interacción, con α = 0,05, bilateral. El tamaño muestral se justifica antes de empezar: ____.
- **Medida secundaria:** fosfo-p38/p38, con el mismo diseño. Sus resultados se reportan como secundarios.
