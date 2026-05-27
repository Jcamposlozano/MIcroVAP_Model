from MicroNAV.application.use_cases.load_microvap_dataset import Load_data
from MicroNAV.utils.logger import get_logger
from MicroNAV.application.use_cases.microvap_analysis import MicroVAPAnalysis

log = get_logger("MicroNAV")

RUN_OTUS = False
RUN_METABOLOMICS = True


def main():
    log.info("Pipeline iniciado.")

    # ── PASO 0: Carga y limpieza del dataset clínico ──────────────────────
    loader = Load_data("microvap_colombia_mothur.metadata_completa.xlsx")
    loader.clean_all_data()
    df = loader.read_file()

    # ── PASO 1: Análisis de microbioma 16S ───────────────────────────────
    # MicroVAPAnalysis resuelve las rutas desde data/raw/ automáticamente.
    # Solo se pasa ddpcr_spec_file si quieres sobrescribir el nombre por defecto.
    log.info("Iniciando análisis de microbioma 16S...")

    analysis = MicroVAPAnalysis()          # rutas automáticas desde data/raw/
    analysis.load_data()

    # Inyecta el DataFrame clínico ya limpio por el loader
    analysis.meta_df = df.copy()
    otu_rel = analysis.otu_good.div(analysis.otu_good.sum(axis=1), axis=0) * 100
    analysis.otu_df = (
        analysis.meta_df
        .set_index("specimen")
        .join(otu_rel, how="inner")
        .reset_index()
        .rename(columns={"index": "specimen"})
    )

    microbiome_results = analysis.run_all()
    log.info("Análisis de microbioma completado.")

    # ── PASO 2: Enriquecimiento del df clínico con índices de diversidad ──
    alpha_baseline = microbiome_results["alpha_diversity"]["baseline"]["raw"]
    df = df.merge(
        alpha_baseline[["specimen", "Shannon", "Simpson_D", "Chao1"]],
        on="specimen",
        how="left",
    )

    print(df.head())    
    '''
    # ── PASO 3: Pipeline clínico original ────────────────────────────────
    steps = [
        loader.update_patient,
        loader.update_admission,
        loader.lab_results,
        loader.update_specimen,
        loader.update_abd_results,
        loader.update_admission_derived,
        loader.update_admission_outcomes,
        loader.update_vent_settings,
        loader.update_vitals,
        loader.update_antibiotic_course,
        loader.update_micro_result,
        loader.update_trauma,
    ]

    if RUN_OTUS:
        steps.append(loader.updateSepecimentOTU)
    if RUN_METABOLOMICS:
        steps.append(loader.updateSpecimenMetabolomic)

    for step in steps:
        step(df)

    log.info("Pipeline finalizado.")
    return df, microbiome_results
    '''

if __name__ == "__main__":
    main()