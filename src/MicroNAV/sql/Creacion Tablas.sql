DROP TABLE IF EXISTS specimen_otu;
DROP TABLE IF EXISTS antibiotic_course;
DROP TABLE IF EXISTS vitals;
DROP TABLE IF EXISTS lab_results;
DROP TABLE IF EXISTS abg_results;
DROP TABLE IF EXISTS vent_settings;
DROP TABLE IF EXISTS admission_outcomes;
DROP TABLE IF EXISTS admission_derived;
DROP TABLE IF EXISTS specimen;
DROP TABLE IF EXISTS dim_otu;
DROP TABLE IF EXISTS dim_antibiotic;
DROP TABLE IF EXISTS admission;
DROP TABLE IF EXISTS patient;

CREATE TABLE patient (
  patient_id INT NOT NULL,
  gender binary,
  height int,
  weight int,
  alcoholism binary,
  anemia binary,
  stroke binary,
  asthma binary,
  cancer binary,
  diabetes binary,
  coronary_heart_disease binary,
  mental_illness binary,
  chronic_kidney_disease binary,
  copd binary,
  heart_failure binary,
  hemodialysis binary,
  hypertension binary,
  ppi_chronic binary,
  obesity binary,
  sahos binary,
  smoking binary,
  statin binary,
  no_background binary,
  PRIMARY KEY (patient_id)
);

CREATE TABLE admission (
  admission_id INT NOT NULL,
  patient_id INT NOT NULL,
  admission_infectious binary,
  admission_general_surgery binary,
  admission_cardiovascular binary,
  admission_trauma binary,
  admission_nephrological binary,
  admission_neurological binary,
  admission_neurosurgical binary,
  admission_autoimmune binary,
  admission_toxic binary,
  admission_metabolic binary,
  other_diagnostic_admission binary,
  PRIMARY KEY (admission_id),
  FOREIGN KEY (patient_id) REFERENCES patient(patient_id)
);

CREATE TABLE specimen (
  specimen_id INT NOT NULL,
  admission_id INT NOT NULL,
  specimen text,
  timepoint int,
  tissue_type varchar(50),
  PRIMARY KEY (specimen_id),
  FOREIGN KEY (admission_id) REFERENCES admission(admission_id)
);

CREATE TABLE micro_result (
  micro_result_id INTEGER PRIMARY KEY AUTOINCREMENT,
  specimen_id INTEGER NOT NULL,
  marker_code TEXT NOT NULL,
  copies_ul_dna REAL,
  FOREIGN KEY (specimen_id) REFERENCES specimen(specimen_id),
  UNIQUE (specimen_id, marker_code)
);


CREATE TABLE dim_otu (
  otu_id INTEGER PRIMARY KEY AUTOINCREMENT,
  otu_code TEXT NOT NULL UNIQUE
);

CREATE TABLE specimen_otu (
  specimen_id INTEGER NOT NULL,
  otu_id INTEGER NOT NULL,
  relative_abundance REAL NOT NULL,  -- ej: 0.117370892
  PRIMARY KEY (specimen_id, otu_id),
  FOREIGN KEY (specimen_id) REFERENCES specimen(specimen_id),
  FOREIGN KEY (otu_id) REFERENCES dim_otu(otu_id)
);

CREATE TABLE dim_antibiotic (
  drug_id INTEGER PRIMARY KEY AUTOINCREMENT,
  drug_name TEXT NOT NULL UNIQUE   -- ej: 'pip-tazo', 'cefepime'
);

CREATE TABLE antibiotic_course (
  abx_course_id INTEGER PRIMARY KEY AUTOINCREMENT,
  admission_id INTEGER NOT NULL,
  drug_id INTEGER NOT NULL,
  start_antibiotic TEXT,   -- si viene como fecha: luego la estandarizas a ISO 'YYYY-MM-DD' o datetime
  end_antibiotic TEXT,
  FOREIGN KEY (admission_id) REFERENCES admission(admission_id),
  FOREIGN KEY (drug_id) REFERENCES dim_antibiotic(drug_id),
  UNIQUE (admission_id, drug_id, start_antibiotic, end_antibiotic)
);

CREATE INDEX idx_abx_admission ON antibiotic_course(admission_id);
CREATE INDEX idx_abx_drug ON antibiotic_course(drug_id);

CREATE TABLE vitals (
  vitals_id INTEGER PRIMARY KEY AUTOINCREMENT,
  admission_id INTEGER NOT NULL,
  measured_at TEXT,   -- idealmente 'YYYY-MM-DD HH:MM:SS'
  heart_rate REAL,
  respiratory_rate REAL,
  temperature REAL,
  sbp REAL,
  dbp REAL,
  mbp REAL,
  oxygen_saturation REAL,
  glasgow_score REAL,

  FOREIGN KEY (admission_id) REFERENCES admission(admission_id)
);

CREATE INDEX idx_vitals_admission ON vitals(admission_id);


CREATE TABLE lab_results (
  lab_id INTEGER PRIMARY KEY AUTOINCREMENT,
  admission_id INTEGER NOT NULL,
  measured_at TEXT,
  white_blood_cells REAL,
  neutrophils REAL,
  hemoglobin REAL,
  hematocrit REAL,
  platelets REAL,
  creatinine REAL,
  uremic_nitrogen REAL,
  glycemia REAL,
  sodium REAL,
  potassium REAL,
  total_bilirubin REAL,
  direct_bilirubin REAL,
  indirect_bilirubin REAL,
  alanine_gpt REAL,
  aspartate_got REAL,
  pcr REAL,
  procalcitonin REAL,
  pt REAL,
  ptt REAL,
  FOREIGN KEY (admission_id) REFERENCES admission(admission_id)
);

CREATE INDEX idx_lab_admission ON lab_results(admission_id);


CREATE TABLE abg_results (
  abg_id INTEGER PRIMARY KEY AUTOINCREMENT,

  admission_id INTEGER NOT NULL,
  measured_at TEXT,

  ph REAL,
  pao2 REAL,
  paco2 REAL,
  fio2 REAL,
  hco3 REAL,
  lactate REAL,
  pafi REAL,
  -- opcional si tu dataset ya la trae como variable
  hiperoxemia TEXT,
  FOREIGN KEY (admission_id) REFERENCES admission(admission_id)
);

CREATE INDEX idx_abg_admission ON abg_results(admission_id);

CREATE TABLE vent_settings (
  vent_id INTEGER PRIMARY KEY AUTOINCREMENT,
  admission_id INTEGER NOT NULL,
  measured_at TEXT,
  ventilatory_mode TEXT,
  tidal_volume REAL,
  peep REAL,
  plateau_pressure REAL,
  peak_pressure REAL,
  ibw_ideal REAL,
  vt_ml_ideal_weight REAL,
  classification_vt_ideal_weight TEXT,
  FOREIGN KEY (admission_id) REFERENCES admission(admission_id)
);

CREATE INDEX idx_vent_admission ON vent_settings(admission_id);

-- 1:1 con admission (PK = FK)
CREATE TABLE admission_outcomes (
  admission_id INTEGER PRIMARY KEY,
  pneumonia INTEGER,
  tracheostomy INTEGER,
  date_tracheostomy TEXT,
  hospital_mortality INTEGER,
  mortality_28d INTEGER,
  mortality_90d INTEGER,
  sdra INTEGER,
  FOREIGN KEY (admission_id) REFERENCES admission(admission_id)
);

-- Opcional 1:1 para derivadas
CREATE TABLE admission_derived (
  admission_id INTEGER PRIMARY KEY,
  days_stay_icu REAL,
  intubation_days REAL,
  length_of_stay REAL,
  days_antibiotic REAL,
  bmi REAL,
  pafi REAL,
  FOREIGN KEY (admission_id) REFERENCES admission(admission_id)
);


