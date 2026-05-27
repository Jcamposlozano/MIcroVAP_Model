# MicroNAV

Analítica MicroNAV

## Estructura (hexagonal-lite)
- `data/`: raw → interim → processed
- `src/MicroNAV/`: lógica (ingesta, features, modelos, evaluación, reporting)
- `notebooks/`: EDA y narrativa (idealmente llamando funciones de `src/`)
- `reports/`: figuras y salidas (HTML/PDF/tablas)
- `artifacts/`: modelos y métricas
- `configs/`: parámetros (`params.yaml`)

## Setup
```bash
poetry install
poetry run python -m MicroNAV.pipeline
```

## Notebooks
Recomendado: mantener notebooks livianos y mover lógica a `src/`.


## Github

### SSH

SSH: git@github.com:Jcamposlozano/MIcroVAP_Model.git

…or create a new repository on the command line
echo "# MIcroVAP_Model" >> README.md
git init
git add README.md
git commit -m "first commit"
git branch -M main
git remote add origin git@github.com:Jcamposlozano/MIcroVAP_Model.git
git push -u origin main

…or push an existing repository from the command line
git remote add origin git@github.com:Jcamposlozano/MIcroVAP_Model.git
git branch -M main
git push -u origin main


### HTTP

Http: https://github.com/Jcamposlozano/MIcroVAP_Model.git

…or create a new repository on the command line
echo "# MIcroVAP_Model" >> README.md
git init
git add README.md
git commit -m "first commit"
git branch -M main
git remote add origin https://github.com/Jcamposlozano/MIcroVAP_Model.git
git push -u origin main

…or push an existing repository from the command line
git remote add origin https://github.com/Jcamposlozano/MIcroVAP_Model.git
git branch -M main
git push -u origin main
