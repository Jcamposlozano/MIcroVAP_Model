from MicroNAV.application.use_cases.load_microvap_dataset import Load_data
from MicroNAV.utils.logger import get_logger

log = get_logger("MicroNAV")

RUN_OTUS = False
RUN_METABOLOMICS = True

def main():
    log.info("Pipeline iniciado.")

    loader = Load_data("microvap_colombia_mothur.metadata_completa.xlsx")

    loader.clean_all_data()

    df = loader.read_file()

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


if __name__ == "__main__":
    main()