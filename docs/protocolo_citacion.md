# Protocolo de citación y archivo

**Objetivo:** que cada versión citada en el manuscrito tenga un DOI permanente, ligado al preregistro y al autor.

## 1. Identificadores

| Identificador | Estado | Dónde se anota |
|---|---|---|
| ORCID 0000-0001-8566-1868 | ✅ | `CITATION.cff`, `README.md` |
| Preregistro OSF | ⏳ el borrador está en `docs/preregistration_osf.md`; falta registrarlo | `CITATION.cff` (`identifiers`), `README.md`, manuscrito |
| Repositorio GitHub | ⏳ por crear | `CITATION.cff` (`repository-code`) |
| DOI Zenodo del código | ⏳ se crea con la primera *release* | `CITATION.cff` (`doi`), insignia en el `README.md` |
| DOI Zenodo de los datos derivados | ⏳ depósito aparte (ver `protocolo_datos.md`, sección 4) | `README.md`, disponibilidad de datos del manuscrito |

## 2. Orden recomendado

1. **Registrar el preregistro en OSF antes de publicar resultados.** Como el repositorio será público y ya hay resultados, conviene declarar en el registro la fecha y la etiqueta git del análisis ya hecho. Es lo honesto, y los revisores lo verán igual en el historial.
2. Crear el repositorio en GitHub y hacer el primer *push* (ver `CONTRIBUTING.md`).
3. Vincular la cuenta de GitHub con Zenodo (zenodo.org → GitHub → activar el repositorio).
4. Publicar la *release* `v1.0-analisis`. Zenodo genera el DOI.
5. Completar `repository-code`, `doi` e `identifiers` en `CITATION.cff` y agregar la insignia del DOI al `README.md`.
6. Crear el depósito de datos derivados en Zenodo y enlazarlo como `isSupplementTo` del DOI del código.

## 3. En el manuscrito

- **Disponibilidad de código:** "El código está disponible en GitHub (URL) y archivado en Zenodo (DOI de la versión v1.1-envio), con licencia MIT."
- **Disponibilidad de datos:** "Los resultados derivados están en Zenodo (DOI, CC BY 4.0). Los registros de Scopus y WoS no pueden redistribuirse; se proporcionan sus identificadores y las ecuaciones de búsqueda."
- **Preregistro:** cita el DOI de OSF y remite a `docs/decisions.md` para las desviaciones.
- Cita siempre el DOI de la **versión** exacta, no el DOI conceptual, que apunta a la última.

## 4. Cuándo publicar una nueva *release*

- Al enviar el manuscrito (`v1.1-envio`).
- Después de cada ronda de revisión que cambie resultados (`v1.2-revision`, …).
- Nunca se reescribe una etiqueta ya publicada: si hace falta corregir algo, se crea una etiqueta nueva.
