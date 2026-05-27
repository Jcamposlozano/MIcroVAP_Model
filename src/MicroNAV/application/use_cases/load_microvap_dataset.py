from MicroNAV.adapters.outbound.db.sqlite_db import Db_datos

from pathlib import Path
import pandas as pd
import sys
from pathlib import Path
import os
import numpy as np
import re

class Load_data:
    
    def __init__(self,file):
        self.PROJECT_ROOT = Path.cwd()
        self.SRC_PATH = self.PROJECT_ROOT / "src"
        sys.path.append(str(self.SRC_PATH))  
        
        self.DB_PATH = self.PROJECT_ROOT / "data" / "microvap.db"
        self.DATA_RAW = self.PROJECT_ROOT / "data" / "raw"
        self.SQL_PATH = self.PROJECT_ROOT / "src" / "MicroNAV" / "sql"

        self.db_data = Db_datos(self.DB_PATH)
        self.file = file

    def clean_all_data(self):
        sql_file = self.SQL_PATH / "Creacion Tablas.sql"

        if not sql_file.exists():
            raise FileNotFoundError(f"No se encontró el archivo SQL: {sql_file}")

        with open(sql_file, "r", encoding="utf-8") as file:
            sql_script = file.read()

        self.db_data.run_query(sql_script)

        print("Base de datos limpiada y tablas creadas correctamente.")

    def build_measured_at(self, row):
        if pd.isna(row["admission_date_icu"]):
            return None
        
        base = pd.to_datetime(row["admission_date_icu"], errors="coerce")

        if pd.isna(base):
            return None

        tp = str(row["timepoint"]).lower()

        if "day" in tp:
            try:
                d = int(tp.replace("day", "").strip())
                return base + pd.Timedelta(days=d)
            except:
                return base

        if "baseline" in tp:
            return base

        return base

    def read_file(self):
        df = pd.read_excel(self.DATA_RAW / self.file)

        df.columns = df.columns.str.strip()

        df = df[df["reagent_tissue"] == "tissue"].copy()

        df["admission_date_icu"] = pd.to_datetime(
            df["admission_date_icu"],
            dayfirst=True,
            errors="coerce"
        ).dt.strftime("%Y-%m-%d")

        df = df.assign(
            admission_id=df["patient_id"].astype(str) + "_" + df["admission_date_icu"]
        )

        df = df.assign(
            measured_at=df.apply(self.build_measured_at, axis=1)
        )

        return df

    def update_patient(self, df):        
        variables_list= ['patient_id',
                        'gender',
                        'height',
                        'weight',
                        'alcoholism',
                        'anemia',
                        'stroke',
                        'asthma',
                        'cancer',
                        'diabetes',
                        'coronary_heart_disease',
                        'mental_illness',
                        'chronic_kidney_disease',
                        'copd',
                        'heart_failure',
                        'hemodialysis',
                        'hypertension',
                        'ppi_chronic',
                        'obesity',
                        'sahos',
                        'smoking',
                        'statin',
                        'no_background']

        df_patient = df[variables_list]
        df_grouped = (
                df_patient
                .groupby('patient_id', as_index=False)
                .agg({
                    'gender': 'first',
                    'height': 'max',
                    'weight': 'max',
                    'alcoholism': 'max',
                    'anemia': 'max',
                    'stroke': 'max',
                    'asthma': 'max',
                    'cancer': 'max',
                    'diabetes': 'max',
                    'coronary_heart_disease': 'max',
                    'mental_illness': 'max',
                    'chronic_kidney_disease': 'max',
                    'copd': 'max',
                    'heart_failure': 'max',
                    'hemodialysis': 'max',
                    'hypertension': 'max',
                    'ppi_chronic': 'max',
                    'obesity': 'max',
                    'sahos': 'max',
                    'smoking': 'max',
                    'statin': 'max',
                    'no_background': 'max'
                })
            )
        try:
            #print(df_grouped.head())
            n = self.db_data.load_df(df_grouped, "patient", clear_table=True)
            print(f"Patient Insertados: {n}")
        except Exception as e:
            print(f"Error: {e}")
            return None

    def update_admission(self,df):
        variables_list= [
                        'patient_id',
                        'admission_id',
                        'admission_infectious',
                        'admission_general_surgery',
                        'admission_cardiovascular',
                        'admission_trauma',
                        'admission_nephrological',
                        'admission_neurological',
                        'admission_neurosurgical',
                        'admission_autoimmune',
                        'admission_toxic',
                        'admission_metabolic',
                        'other_diagnostic_admission']

        df_admission = df[variables_list]

        df_grouped = (
            df_admission
            .groupby('admission_id', as_index=False)
            .agg({
                'patient_id': 'first',
                'admission_id': 'max',
                'admission_infectious': 'max',
                'admission_general_surgery': 'max',
                'admission_cardiovascular': 'max',
                'admission_trauma': 'max',
                'admission_nephrological': 'max',
                'admission_neurological': 'max',
                'admission_neurosurgical': 'max',
                'admission_autoimmune': 'max',
                'admission_toxic': 'max',
                'admission_metabolic': 'max',
                'other_diagnostic_admission': 'max'
            })
        )

        n = self.db_data.load_df(df_grouped, "admission", clear_table=True)
        print(f"Admission Insertados: {n}")

    def lab_results(self, df):
        variables_lab = [
            "admission_id",
            "measured_at",
            "white_blood_cells",
            "neutrophils",
            "hemoglobin",
            "hematocrit",
            "platelets",
            "creatitin",
            "uremic_nitrogen",
            "glycemia",
            "sodium",
            "potassium",
            "total_bilirubin",
            "direct_bilirubin",
            "indirect_bilirubin",
            "alanine_gpt",
            "aspartate_got",
            "pcr",
            "procalcitonin",
            "pt",
            "ptt"
        ]

        df_lab_results = df[variables_lab]

        df_lab_results = df_lab_results.rename(columns={
            "creatitin": "creatinine"
        })

        df_lab_results_group = (
            df_lab_results.groupby("admission_id", as_index=False)
            .agg({
                "measured_at": 'first',
                "white_blood_cells": 'first',
                "neutrophils": 'first',
                "hemoglobin": 'first',
                "hematocrit": 'first',
                "platelets": 'first',
                "creatinine": 'first',
                "uremic_nitrogen": 'first',
                "glycemia": 'first',
                "sodium": 'first',
                "potassium": 'first',
                "total_bilirubin": 'first',
                "direct_bilirubin": 'first',
                "indirect_bilirubin": 'first',
                "alanine_gpt": 'first',
                "aspartate_got": 'first',
                "pcr": 'first',
                "procalcitonin": 'first',
                "pt": 'first',
                "ptt": 'first'
            })
        )

        n2 = self.db_data.load_df(df_lab_results_group, "lab_results", clear_table=True)
        print("lab_results insertados:", n2)

    def update_specimen(self, df):
        variables_list= [
                    'specimen',
                    'admission_id',
                    'timepoint',
                    'tissue_type'
                    ]

        df_specimen = df[variables_list]

        df_specimen = df_specimen.rename(columns={
            'specimen': 'specimen_id'
        })

        df_grouped = (
            df_specimen
            .groupby('specimen_id', as_index=False)
            .agg({
                'admission_id': 'max',
                'timepoint': 'first',
                'tissue_type': 'first'
            })
        )

        n = self.db_data.load_df(df_grouped, "specimen", clear_table=True)
        print(f"Specimen Insertados: {n}")

    def normalize_otu(self, x):
        if pd.isna(x):
            return None
        s = str(x).strip().lower()
        m = re.search(r"otu\D*(\d+)", s) or re.search(r"(\d+)", s)
        if not m:
            return None
        return f"otu{int(m.group(1))}"

    def updateSepecimentOTU(self, df):
        # -------------------------
        # 1) Preparar df_long (specimen_id, otu_code, relative_abundance)
        # -------------------------
        df_work = df.copy().rename(columns={"specimen": "specimen_id"})
        df_work["specimen_id"] = df_work["specimen_id"].astype(str).str.strip()

        # Detecta columnas tipo: otu1, OTU_1, Otu00001, otu-0002, etc.
        otu_cols = [
            c for c in df_work.columns
            if re.match(r"(?i)^otu\D*\d+$", str(c).strip())
        ]
        if not otu_cols:
            raise ValueError("No encontré columnas OTU (ej: Otu00001, otu1, OTU_12, otu-003 ...)")

        df_long = df_work.melt(
            id_vars=["specimen_id"],
            value_vars=otu_cols,
            var_name="otu_code",
            value_name="relative_abundance",
        )

        df_long["relative_abundance"] = pd.to_numeric(df_long["relative_abundance"], errors="coerce")
        df_long = df_long.dropna(subset=["relative_abundance"])
        df_long = df_long[df_long["relative_abundance"] > 0].copy()

        df_long["otu_code"] = df_long["otu_code"].astype(str).str.strip()
        df_long["otu_code_norm"] = df_long["otu_code"].map(self.normalize_otu)

        # Si normalize_otu no pudo extraer número, lo botamos (o puedes levantar error)
        df_long = df_long.dropna(subset=["otu_code_norm"]).copy()


        # -------------------------
        # 2) Construir / actualizar catálogo dim_otu DESDE el dataset
        # -------------------------
        df_catalog = (
            df_long[["otu_code_norm"]]
            .drop_duplicates()
            .rename(columns={"otu_code_norm": "otu_code"})
        )

        # Si estás reconstruyendo el catálogo desde cero:
        # OJO: por FKs, primero limpiar specimen_otu, luego dim_otu.
        # Si no tienes db_data.execute, deja clear_table=True en dim_otu y specimen_otu y ya.
        # db_data.execute("DELETE FROM specimen_otu;")
        # db_data.execute("DELETE FROM dim_otu;")


        df_catalog['description'] = ''
        df_catalog['comment'] = ''

        n_cat = self.db_data.load_df(df_catalog, "dim_otu", clear_table=True)
        print("dim_otu insertados:", n_cat)

        # -------------------------
        # 3) Traer dim_otu ya con otu_id y preparar llave norm (por seguridad)
        # -------------------------
        df_dim_otu = self.db_data.extracData("SELECT otu_id, otu_code FROM dim_otu")
        df_dim_otu["otu_code_norm"] = df_dim_otu["otu_code"].map(self.normalize_otu)

        # Evita que un dim_otu sucio te duplique filas al hacer merge
        df_dim_otu = df_dim_otu.dropna(subset=["otu_code_norm"]).drop_duplicates(subset=["otu_code_norm"])


        # -------------------------
        # 4) Armar df_specimen_otu final (specimen_id, otu_id, relative_abundance)
        # -------------------------
        df_specimen_otu = df_long.merge(
            df_dim_otu[["otu_id", "otu_code_norm"]],
            on="otu_code_norm",
            how="left",
        )

        # Validación: si esto falla, dim_otu no tiene esas OTUs (o normalize_otu no coincide)
        if df_specimen_otu["otu_id"].isna().any():
            missing_norm = (
                df_specimen_otu.loc[df_specimen_otu["otu_id"].isna(), "otu_code_norm"]
                .dropna()
                .unique()
            )
            raise ValueError(f"Siguen faltando OTUs en dim_otu (norm): {missing_norm[:20]}")

        df_specimen_otu = df_specimen_otu[["specimen_id", "otu_id", "relative_abundance"]].copy()

        # Limpieza de tipos
        df_specimen_otu["specimen_id"] = df_specimen_otu["specimen_id"].astype(str).str.strip()
        df_specimen_otu["otu_id"] = df_specimen_otu["otu_id"].astype(int)
        df_specimen_otu["relative_abundance"] = pd.to_numeric(df_specimen_otu["relative_abundance"], errors="coerce")
        df_specimen_otu = df_specimen_otu.dropna(subset=["specimen_id", "otu_id", "relative_abundance"]).copy()

        # -------------------------
        # 4.1) FIX CLAVE: colapsar duplicados (specimen_id, otu_id)
        #     Esto evita: UNIQUE constraint failed
        # -------------------------
        df_specimen_otu = (
            df_specimen_otu
            .groupby(["specimen_id", "otu_id"], as_index=False)["relative_abundance"]
            .sum()
        )

        # Guard rail
        if df_specimen_otu.duplicated(["specimen_id", "otu_id"]).any():
            raise ValueError("Aún hay duplicados en (specimen_id, otu_id) después del groupby. Revisa el merge/normalización.")

    def update_abd_results(self, df):
        variables_list= [
            'admission_id',
            'measured_at', 
            'ph',
            'pao2',
            'paco2',
            'fio2',
            'hco3',
            'lactate',
            'pafi',
            'hiperoxemia'            
            ]

        df_abg_results = df[variables_list]

        df_abg_results_grouped = (
            df_abg_results
            .groupby('admission_id', as_index=False)
            .agg({
                'measured_at': 'first',
                'ph': 'first',
                'pao2': 'first',
                'paco2': 'first',
                'fio2': 'first',
                'hco3': 'first',
                'lactate': 'first',
                'pafi': 'first',
                'hiperoxemia': 'first' 
            })
        )
        n2 = self.db_data.load_df(df_abg_results_grouped, "abg_results", clear_table=True)
        print("Abg_results insertados:", n2)

    def update_admission_derived(self,df):
        variables_list= [
           "admission_id",
            "days_stay_icu",
            "intubation_days",
            "length_stay",
            "days_antibiotic",
            "bmi",
            "pafi",          
            ]
        df_admission_derived = df[variables_list]

        df_admission_derived_agg = (
                df_admission_derived
                .groupby('admission_id', as_index=False)
                .agg({
                    "days_stay_icu": "first",
                    "intubation_days": "first",
                    "length_stay": "first",
                    "days_antibiotic": "first",
                    "bmi": "first",
                    "pafi": "first"
                })
            )

        n2 = self.db_data.load_df(df_admission_derived_agg, "admission_derived", clear_table=True)
        print("Admission_derived insertados:", n2)

    def update_admission_outcomes(self, df):
        variables_list= [
                "admission_id",
                    "pneumonia",
                    "tracheostomy",
                    "date_tracheostomy",
                    "hospital_mortality",
                    "mortality_28d",
                    "mortality_90d",
                    "sdra"   
                    ]

        df_admission_outcomes = df[variables_list]

        df_admission_outcomes_agg = (
            df_admission_outcomes.groupby("admission_id", as_index=False)
            .agg({
                "pneumonia": 'first',
                "tracheostomy": 'first',
                "date_tracheostomy": 'first',
                "hospital_mortality": 'first',
                "mortality_28d": 'first',
                "mortality_90d": 'first',
                "sdra": 'first'   
            })
        )

        n2 = self.db_data.load_df(df_admission_outcomes_agg, "admission_outcomes", clear_table=True)
        print("admission_outcomes insertados:", n2)

    def update_vent_settings(self, df):
        variables_vent = [
            "admission_id",
            "measured_at",
            "ventilatory_mode",
            "tidal_volumen",   # luego renombramos
            "peep",
            "plateau_pressure",
            "peak_pressure",
            "ibw_ideal",
            "vt_ml_ideal_weight",
            "classification_vt_ideal_weight"
        ]

        df_variables_vent = df[variables_vent]

        df_variables_vent["tidal_volume"] = df["tidal_volumen"]

        df_variables_vent_group = (
            df_variables_vent.groupby("admission_id", as_index=False)
            .agg({
                "ventilatory_mode": 'first',
                "measured_at": 'first',
                "tidal_volume": 'first',
                "peep": 'first',
                "plateau_pressure": 'first',
                "peak_pressure": 'first',
                "ibw_ideal": 'first',
                "vt_ml_ideal_weight": 'first',
                "classification_vt_ideal_weight": 'first',

            })
        )

        n2 = self.db_data.load_df(df_variables_vent_group, "vent_settings", clear_table=True)
        print("vent_settings insertados:", n2)

    def update_vitals(self, df):
        variables_vitals = [
            "admission_id",
            "measured_at",
            "heart_rate",
            "respiratory_rate",
            "temperature",
            "sbp",
            "dbp",
            "mbp",
            "saturation",
            "glasgow"  
        ]

        df_vitals = df[variables_vitals]

        df_vitals_group = (
            df_vitals.groupby("admission_id", as_index=False)
            .agg({
                "heart_rate": 'first',
                "measured_at": 'first',
                "respiratory_rate": 'first',
                "temperature": 'first',
                "sbp": 'first',
                "dbp": 'first',
                "mbp": 'first',
                "saturation": 'first',
                "glasgow": 'first',
            })
        )

        df_vitals_group["oxygen_saturation"] = df["saturation"]
        df_vitals_group["glasgow_score"] = df["glasgow"]

        n2 = self.db_data.load_df(df_vitals_group, "vitals", clear_table=True)
        print("vitals insertados:", n2)

    def normalize_antibiotic(self, x):
        if pd.isna(x):
            return None
        
        s = str(x).strip().lower()
        
        if s in ["", "na", "nan", "none", "null"]:
            return None
        
        # normalizaciones básicas
        s = s.replace("/", "-")
        s = re.sub(r"\s+", "-", s)
        s = re.sub(r"-+", "-", s)
        
        return s

    def update_antibiotic_course(self, df):
        # -------------------------
        # 1) Preparar dataframe base
        # -------------------------
        abx_cols = [
            "type_antibiotic_1",
            "type_antibiotic_2",
            "type_antibiotic_3",
        ]

        base_cols = [
            "admission_id",
            "start_antibiotic",
            "end_antibiotic",
        ] + abx_cols

        df_abx_work = df[base_cols].copy()

        df_abx_work["admission_id"] = df_abx_work["admission_id"].astype(str).str.strip()

        # -------------------------
        # 2) Wide -> Long
        # -------------------------
        df_abx_long = df_abx_work.melt(
            id_vars=["admission_id", "start_antibiotic", "end_antibiotic"],
            value_vars=abx_cols,
            var_name="antibiotic_slot",
            value_name="drug_name"
        )

        # -------------------------
        # 3) Limpiar antibióticos
        # -------------------------
        df_abx_long["drug_name"] = df_abx_long["drug_name"].map(self.normalize_antibiotic)

        df_abx_long = df_abx_long.dropna(subset=["drug_name"]).copy()

        # -------------------------
        # 4) Crear catálogo dim_antibiotic
        # -------------------------
        df_dim_antibiotic = (
            df_abx_long[["drug_name"]]
            .drop_duplicates()
            .sort_values("drug_name")
            .reset_index(drop=True)
        )

        # Cargar catálogo
        n_cat = self.db_data.load_df(
            df_dim_antibiotic,
            "dim_antibiotic",
            clear_table=True
        )

        print("dim_antibiotic insertados:", n_cat)

        # -------------------------
        # 5) Traer catálogo con drug_id
        # -------------------------
        df_dim_db = self.db_data.extracData("""
            SELECT drug_id, drug_name
            FROM dim_antibiotic
        """)

        df_dim_db["drug_name"] = df_dim_db["drug_name"].map(self.normalize_antibiotic)

        # -------------------------
        # 6) Mapear drug_name -> drug_id
        # -------------------------
        df_antibiotic_course = df_abx_long.merge(
            df_dim_db,
            on="drug_name",
            how="left"
        )

        if df_antibiotic_course["drug_id"].isna().any():
            missing = df_antibiotic_course.loc[
                df_antibiotic_course["drug_id"].isna(),
                "drug_name"
            ].unique()
            raise ValueError(f"Antibióticos no encontrados en dim_antibiotic: {missing[:20]}")

        # -------------------------
        # 7) Preparar columnas finales
        # -------------------------
        df_antibiotic_course = df_antibiotic_course[
            [
                "admission_id",
                "drug_id",
                "start_antibiotic",
                "end_antibiotic",
            ]
        ].copy()

        df_antibiotic_course["drug_id"] = df_antibiotic_course["drug_id"].astype(int)

        # Fechas como TEXT ISO si se puede
        df_antibiotic_course["start_antibiotic"] = pd.to_datetime(
            df_antibiotic_course["start_antibiotic"],
            errors="coerce"
        ).dt.strftime("%Y-%m-%d")

        df_antibiotic_course["end_antibiotic"] = pd.to_datetime(
            df_antibiotic_course["end_antibiotic"],
            errors="coerce"
        ).dt.strftime("%Y-%m-%d")

        # -------------------------
        # 8) Eliminar duplicados para evitar UNIQUE constraint
        # -------------------------
        df_antibiotic_course = df_antibiotic_course.drop_duplicates(
            subset=[
                "admission_id",
                "drug_id",
                "start_antibiotic",
                "end_antibiotic",
            ]
        ).copy()

        # -------------------------
        # 10) Cargar antibiotic_course
        # -------------------------
        n_course = self.db_data.load_df(
            df_antibiotic_course,
            "antibiotic_course",
            clear_table=True
        )

        print("antibiotic_course insertados:", n_course)

    def update_micro_result(self, df):
        df_micro = df[
            [
                "specimen",
                "copies_ul_dna",
                "aureus_presence"
            ]
        ].copy()

        df_micro = df_micro.rename(columns={
            "specimen": "specimen_id",
            "copies_ul_dna": "copies_ul_dna"
        })

        df_micro["specimen_id"] = df_micro["specimen_id"].astype(str).str.strip()

        df_micro["aureus_presence"] = (
            df_micro["aureus_presence"]
            .astype(str)
            .str.strip()
            .str.lower()
        )

        # Modelo sparse:
        # solo insertamos si aureus_presence == yes
        df_micro = df_micro[df_micro["aureus_presence"] == "yes"].copy()

        df_micro["marker_code"] = "aureus"

        df_micro["copies_ul_dna"] = pd.to_numeric(
            df_micro["copies_ul_dna"],
            errors="coerce"
        )

        df_micro_result = df_micro[
            [
                "specimen_id",
                "marker_code",
                "copies_ul_dna"
            ]
        ].copy()

        # Evitar duplicados por UNIQUE(specimen_id, marker_code)
        df_micro_result = df_micro_result.drop_duplicates(
            subset=["specimen_id", "marker_code"]
        )

        n_course = self.db_data.load_df(
            df_micro_result,
            "micro_result",
            clear_table=True
        )
        print("Micro result insertados:", n_course)

    def update_trauma(self, df):
        variables_vent = [
            "admission_id",
            "injury_severity_score_iss", 	
            "trauma",
            "multiple_trauma", 	
            "craniofacial_trauma",
            "traumatic_brain_injury",
            "other",
            "skull",
            "thorax",
            "abdomen",
            "upper_limbs",
            "lower_limbs"
        ]

        df_trauma = df[variables_vent]

        df_trauma_group = (
            df_trauma.groupby("admission_id", as_index=False)
            .agg({
                "injury_severity_score_iss": 'first', 	
                "trauma": 'first',
                "multiple_trauma": 'first', 	
                "craniofacial_trauma": 'first',
                "traumatic_brain_injury": 'first',
                "other": 'first',
                "skull": 'first',
                "thorax": 'first',
                "abdomen": 'first',
                "upper_limbs": 'first',
                "lower_limbs": 'first'
            })
        )

        n2 = self.db_data.load_df(df_trauma_group, "trauma", clear_table=True)
        print("Trauma insertados:", n2)

    def updateSpecimenMetabolomic(self, df):
        """
        Carga datos metabolómicos en:
        - dim_metabolic
        - metabolomic

        Espera un DataFrame ancho donde:
        - exista la columna specimen
        - las columnas metabolómicas sean todas las columnas posteriores o seleccionadas
        mediante una lista de exclusión.
        """

        df_work = df.copy().rename(columns={"specimen": "specimen_id"})
        df_work["specimen_id"] = df_work["specimen_id"].astype(str).str.strip()

        # Columnas que NO son metabolitos
        exclude_cols = {
            "specimen_id",
            "specimen",
            "name",
            "samples_complete",
            "samples",
            "reagent_tissue",
            "copies_ul_dna",
            "patient_id",
            "admission_date_icu",
            "admission_id",
            "measured_at",
            "timepoint",
            "timepoint_character",
            "specimen_type",
            "tissue_type",
        }

        metabolite_cols = [
            c for c in df_work.columns
            if str(c).strip() not in exclude_cols
        ]

        if not metabolite_cols:
            raise ValueError("No encontré columnas metabolómicas para procesar.")

        # Pasar de formato ancho a largo
        df_long = df_work.melt(
            id_vars=["specimen_id"],
            value_vars=metabolite_cols,
            var_name="metabolito",
            value_name="valor",
        )

        df_long["metabolito"] = (
            df_long["metabolito"]
            .astype(str)
            .str.strip()
        )

        # Limpiar valores tipo na, NA, vacío, etc.
        df_long["valor"] = df_long["valor"].replace(
            ["na", "NA", "Na", "nan", "NaN", "", " "],
            np.nan
        )

        df_long["valor"] = pd.to_numeric(
            df_long["valor"],
            errors="coerce"
        )

        df_long = df_long.dropna(
            subset=["specimen_id", "metabolito", "valor"]
        ).copy()

        # Opcional: eliminar valores cero
        # Si cero significa ausencia, deja esto activo.
        # Si cero es un valor válido, comenta estas dos líneas.
        df_long = df_long[df_long["valor"] > 0].copy()

        # -------------------------
        # 1) Crear catálogo dim_metabolic
        # -------------------------
        df_catalog = (
            df_long[["metabolito"]]
            .drop_duplicates()
            .sort_values("metabolito")
            .reset_index(drop=True)
        )

        n_cat = self.db_data.load_df(
            df_catalog,
            "dim_metabolic",
            clear_table=True
        )

        print("dim_metabolic insertados:", n_cat)

        # -------------------------
        # 2) Leer catálogo con IDs
        # -------------------------
        df_dim_metabolic = self.db_data.extracData(
            "SELECT metabolito_id, metabolito FROM dim_metabolic"
        )

        df_dim_metabolic["metabolito"] = (
            df_dim_metabolic["metabolito"]
            .astype(str)
            .str.strip()
        )

        df_dim_metabolic = df_dim_metabolic.drop_duplicates(
            subset=["metabolito"]
        )

        # -------------------------
        # 3) Merge para obtener metabolito_id
        # -------------------------
        df_metabolomic = df_long.merge(
            df_dim_metabolic[["metabolito_id", "metabolito"]],
            on="metabolito",
            how="left"
        )

        if df_metabolomic["metabolito_id"].isna().any():
            missing = (
                df_metabolomic
                .loc[df_metabolomic["metabolito_id"].isna(), "metabolito"]
                .dropna()
                .unique()
            )
            raise ValueError(
                f"Faltan metabolitos en dim_metabolic: {missing[:20]}"
            )

        df_metabolomic = df_metabolomic[
            ["specimen_id", "metabolito_id", "valor"]
        ].copy()

        df_metabolomic["specimen_id"] = (
            df_metabolomic["specimen_id"]
            .astype(str)
            .str.strip()
        )

        df_metabolomic["metabolito_id"] = (
            df_metabolomic["metabolito_id"]
            .astype(int)
        )

        df_metabolomic["valor"] = pd.to_numeric(
            df_metabolomic["valor"],
            errors="coerce"
        )

        df_metabolomic = df_metabolomic.dropna(
            subset=["specimen_id", "metabolito_id", "valor"]
        ).copy()

        # -------------------------
        # 4) Colapsar duplicados
        # -------------------------
        df_metabolomic = (
            df_metabolomic
            .groupby(["specimen_id", "metabolito_id"], as_index=False)["valor"]
            .sum()
        )

        if df_metabolomic.duplicated(["specimen_id", "metabolito_id"]).any():
            raise ValueError(
                "Aún hay duplicados en (specimen_id, metabolito_id)."
            )

        # -------------------------
        # 5) Cargar tabla final
        # -------------------------
        n_met = self.db_data.load_df(
            df_metabolomic,
            "metabolomic",
            clear_table=True
        )

        print("metabolomic insertados:", n_met)

        return df_metabolomic