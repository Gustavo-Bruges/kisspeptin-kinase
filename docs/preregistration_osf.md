# Preregistro (borrador) · Plantilla OSF de análisis de datos secundarios

**Título:** Literature and transcriptomic evidence on kinase signaling downstream of KISS1R: an integrated bibliometric, gene co-occurrence and causal network analysis

**Autor:** Gustavo Bruges (UNAB) · ORCID 0000-0001-8566-1868
**Fecha prevista de registro:** ____ · **Embargo:** hasta la publicación

---

## 1. Pregunta e hipótesis

**Pregunta.** ¿Qué vías de fosforilación río abajo de KISS1R respaldan los datos transcriptómicos públicos, y cuáles de ellas no ha explorado la literatura?

- **H1.** La literatura sobre kisspeptina forma al menos dos comunidades (oncológica y neuroendocrina) que estudian conjuntos distintos de quinasas y presentan baja conectividad entre sí.
- **H2.** La inferencia causal sobre series con KISS1R perturbado recupera ERK1/2, p38, PI3K/AKT y GSK3β con mayor frecuencia que una perturbación ficticia.
- **H3.** Existen quinasas con alto soporte en datos y bajo soporte en literatura; las quinasas descritas en la comunidad oncológica aparecen en la red hipotalámica más de lo esperado por azar.

## 2. Datos

### 2.1 Corpus bibliográfico
- Bases: ____ (principal), ____ y OpenAlex (sensibilidad).
- Ecuaciones: `docs/search_equations_v1.1.txt`. Periodo 1996–2025; artículos y revisiones.
- Fecha de extracción: ____.

### 2.2 Datos transcriptómicos
- Series: GSE25458 (agonista, signo +1), GSE28383 (knockout, signo −1), más las adicionales que cumplan los criterios: perturbación directa de KISS1/KISS1R, controles en la misma serie, n ≥ 3 por grupo.

### 2.3 Recursos de conocimiento
- PubTator3 (anotación de genes), OmniPath (red previa con signo), CollecTRI y PROGENy, referencia del kinoma ____, COSMIC CGC ____ y CancerMine ____.

## 3. Conocimiento previo de los datos (declaración obligatoria)

Antes del registro, el autor:
- Leyó las descripciones públicas de GSE25458 y GSE28383, sin descargar ni analizar los datos de expresión.
- Ejecutó una consulta preliminar en OpenAlex que devolvió del orden de 30 páginas de resultados (~5.000–6.000 registros) sin analizar su contenido.
- Consultó síntesis de literatura sobre quinasas en cáncer asociadas a KISS1/KISS1R, que motivaron la ampliación del bloque K a la versión 1.1 (EIF2AK2, AXL, EGFR, TGFBR, ROCK, Src family).
- Revisó en PubMed que no existe un análisis bibliométrico dedicado a kisspeptina ni estudios de fosfoproteómica de KISS1R.

## 4. Variables y procedimientos

- **Comunidades:** acoplamiento bibliográfico + Louvain; etiquetado semiautomático por diccionario y confirmación manual.
- **Conectividad:** peso de acoplamiento entre comunidades / peso total, por quinquenio.
- **Subcorpus S:** marca local con el bloque K v1.1 sobre título y resumen.
- **Contexto oncológico:** diccionario de tipos y subtipos de cáncer; dirección del efecto por anotación manual doble (kappa de Cohen).
- **Genes:** PubTator3, símbolo humano oficial (org.Hs.eg.db + ortólogos).
- **Coocurrencia:** oración como unidad; NPMI; hipergeométrica con FDR; backbone SDSM.
- **lit_score:** log2 del enriquecimiento de cada quinasa en el corpus frente a su frecuencia en PubMed.
- **Soporte en datos:** fracción de series en cuya red causal (CARNIVAL sobre OmniPath sin aristas circulares, nodos estables en ≥ 70 % de 100 submuestreos al 80 %) aparece la quinasa.

## 5. Umbrales fijados

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

- Si hay menos de 3 series perturbadas elegibles en la semana 4, el módulo de datos se presenta como exploratorio o se separa en un segundo artículo.
- Si no hay segundo anotador, la vista supresor/promotor se restringe a una muestra estratificada o se elimina.

## 9. Validación

- Precisión de PubTator3 y precisión relacional en 150 resúmenes aleatorios de S, con IC 95 %.
- Precisión de la ecuación en 200 registros aleatorios por base; sensibilidad frente a un conjunto de referencia de 25–30 artículos definido antes de la búsqueda.
