from MicroNAV.utils.logger import get_logger

log = get_logger("MicroNAV")

def main():
    log.info("Pipeline iniciado.")
    # TODO:
    # 1) Ingesta (raw -> interim)
    # 2) Preprocesamiento (interim -> processed)
    # 3) Features
    # 4) Entrenamiento / evaluación
    # 5) Export de métricas / reportes
    log.info("Pipeline finalizado (placeholder).")

if __name__ == "__main__":
    main()
