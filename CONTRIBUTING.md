# Flujo de trabajo con git

El repositorio es **público**. Cualquier cosa que entre al historial debe tratarse como publicada para siempre, aunque después se borre. Revisa `docs/protocolo_datos.md` antes del primer *commit*.

## 1. Configuración inicial (una vez por clon)

```sh
git config core.hooksPath .githooks     # activa las comprobaciones de pre-commit
```

## 2. Ramas

| Rama | Uso |
|---|---|
| `main` | Siempre ejecutable de punta a punta. Sólo recibe *merges* desde ramas revisadas |
| `analysis/<tema>` | Cambios que alteran resultados (por ejemplo `analysis/sigpr54-stability`) |
| `docs/<tema>` | Documentos y manuscrito, sin cambios de resultados |
| `fix/<tema>` | Correcciones de errores |

Cada rama se integra con un *pull request*, aunque no haya más colaboradores. El PR deja constancia de qué cambió y por qué.

## 3. Mensajes de commit

Formato: `<tipo>(<script o área>): <qué cambió>`, en una línea de 72 caracteres como máximo, con el cuerpo opcional después de una línea en blanco.

| Tipo | Cuándo |
|---|---|
| `analysis` | Cambia un resultado o un parámetro (dejar la referencia a la fila de `docs/decisions.md`) |
| `fix` | Corrige un error sin cambiar la intención del análisis |
| `docs` | Documentos, protocolos o manuscrito |
| `data` | Anotación manual, listas curadas, exclusiones |
| `chore` | `.gitignore`, `renv.lock`, *hooks*, estructura |

Ejemplos:
- `analysis(07): fake receptor stability with paired seeds (decisions 2026-10-19)`
- `fix(08): compute series sign outside transmute()`
- `data(annotation): consolidated direction annotation, kappa 0.81`

## 4. Cambios que afectan al preregistro

Si un cambio modifica un parámetro de `params` en `R/00_setup.R` o una regla de análisis:

1. Agrega una fila en `docs/decisions.md` con la fecha, la decisión y el motivo.
2. Si contradice el preregistro, agrega además una fila en la tabla **Desviaciones posteriores al preregistro**.
3. En el mensaje del *commit* y en el PR, cita esa fila.

Un cambio así nunca va en un *commit* de tipo `fix` o `docs`.

## 5. Versiones (etiquetas)

| Etiqueta | Momento |
|---|---|
| `v0.1-preregistro` | Estado del código cuando se registra el preregistro en OSF |
| `v1.0-analisis` | Primer análisis completo (`00`–`08`) |
| `v1.1-envio` | Versión enviada a la revista o publicada como *preprint* |
| `v1.2-revision` | Después de responder a los revisores |

Cada etiqueta se publica como *release* de GitHub, y Zenodo la archiva con un DOI. En el manuscrito se cita el DOI de la versión exacta.

## 6. Antes de cada *push*

- [ ] El *hook* pasó sin `--no-verify`.
- [ ] `git status` no muestra archivos de `data/raw/` ni de `data/processed/`.
- [ ] Si cambió un resultado, se volvió a correr el script y `output/` está actualizado en el mismo *commit*.
- [ ] Si se instaló o actualizó un paquete, `renv::snapshot()` está hecho y `renv.lock` va en el *commit*.
