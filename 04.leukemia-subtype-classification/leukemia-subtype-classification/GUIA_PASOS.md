# Guía paso a paso (no subir a GitHub)

## 1. Preparar el ordenador (30 min)
1. Instala R y RStudio.
2. Abre esta carpeta en RStudio y ejecuta en la consola: `source("install_packages.R")`
3. Comprueba los datos: `library(ALL); data(ALL); dim(ALL); table(ALL$mol.biol)`

## 2. Ejecutar (5–10 min)
En la pestaña Terminal de RStudio, desde la carpeta del proyecto: `Rscript run_all.R`

Si algo falla, copia el error completo y pásamelo.

## 3. Mirar los resultados (lo importante, 1–2 h)
- `results/tables/classifier_performance.csv`: el control con etiquetas barajadas debe quedar cerca del azar. Si no, hay un problema.
- `results/tables/confusion_matrix.csv`: ¿qué subtipos se confunden? Mira sobre todo NEG.
- `results/tables/cluster_stability.csv` y `cluster_vs_subtype.csv`: ¿el clustering recupera los subtipos conocidos?
- Figuras en `results/figures/`.

## 4. Completar el README
Sustituye cada TODO por tus números reales y escribe con tus palabras una observación. No pongas nada que no hayas visto en tus resultados.

## 5. Probar un cambio tuyo
Cambia `n_genes` a 500 o 2000 en los dos scripts y `alpha` en `03_classify.R`, y anota qué pasa. Es lo que te convierte en autora y te da algo real que contar.

## 6. Subir a GitHub
Crea en github.com un repositorio público vacío llamado `leukemia-subtype-classification`. Después:

    git init
    git add .
    git reset GUIA_PASOS.md
    git commit -m "Add ALL subtype discovery and classification workflow"
    git branch -M main
    git remote add origin https://github.com/pchegom-cmd/leukemia-subtype-classification.git
    git push -u origin main

## Lo que te pueden preguntar
- **¿Por qué seleccionar genes dentro de cada pliegue?** Si se usan todas las muestras, las de test influyen en el modelo y la precisión sale inflada (data leakage).
- **¿Por qué balanced accuracy?** Con clases desiguales, un modelo que acierte solo la mayoritaria parecería bueno.
- **¿Por qué elastic net?** Hay muchos más genes que muestras; combina selección de variables (lasso) y estabilidad con genes correlacionados (ridge).
- **¿Por qué PAC?** Mide con qué frecuencia las parejas de muestras cambian de grupo al submuestrear: valores bajos indican clusters más estables.
- **¿Por qué UMAP solo para visualizar?** Distorsiona distancias globales y depende de sus parámetros; no define por sí solo los subtipos.
- **¿Cuál es la limitación?** Una sola cohorte, microarray, grupos pequeños y sin validación externa.
- **¿Hubo batch correction?** No: es una sola cohorte. Con varias, como en LUCID-AML, sería un paso central.

Sé sincera si te preguntan: es un proyecto de aprendizaje con datos públicos, hecho con ayuda de IA y que has ejecutado y modificado tú.
