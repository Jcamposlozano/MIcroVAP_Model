import sqlite3
import pandas as pd
#from oauth2client.service_account import ServiceAccountCredentials

import os
import sys
import traceback


def get_resource_path(relative_path):
    """Obtener la ruta del archivo en entorno de desarrollo o en el ejecutable."""
    try:
        # Para el ejecutable de PyInstaller
        if getattr(sys, 'frozen', False):
            # Si el programa está ejecutándose desde el ejecutable
            base_path = sys._MEIPASS
        else:
            # Si se está ejecutando desde el entorno de desarrollo
            base_path = os.path.abspath(".")

        return os.path.join(base_path, relative_path)
    except Exception as e:
        print(f"Error al obtener la ruta del recurso: {e}")
        return None

class Db_datos:
    def __init__(self, db_path: str):
        self.DB_PATH = db_path
        print(self.DB_PATH)
        #print(f"Ruta de la base de datos: {self.DB_PATH}")

    def test(self):
        return "ok desde aqui"
    
    def run_query(self, query='', params=None):
        """Ejecuta una consulta en la base de datos SQLite."""
        try:
            con = sqlite3.connect(self.DB_PATH)
            cur = con.cursor()
            if query.strip().upper().startswith(('SELECT', 'WITH')):
                cur.execute(query, params or [])
                data = cur.fetchall()
            else:
                cur.execute(query, params or [])
                con.commit()
                data = None
            cur.close()
            con.close()
            return data
        except Exception as e:
            print(f"Error: {e}")
            return None

    def extracData(self, sql):
        """Extrae datos de la base de datos y los retorna como un DataFrame."""
        try:
            con = sqlite3.connect(self.DB_PATH)
            df_respuesta = pd.read_sql(sql, con)
            con.close()
            return df_respuesta
        except Exception as e:
            print(f"Error extracting data: {e}")
            return pd.DataFrame()

    def load_df(
        self,
        df: pd.DataFrame,
        table: str,
        *,
        clear_table: bool = False,
        if_exists: str = "append",
        chunksize: int = 5000,
        method: str = "multi",
        validate_columns: bool = True,
    ) -> int:
        if df is None:
            raise ValueError("df is None")
        if df.empty:
            return 0
        if if_exists not in ("append", "replace"):
            raise ValueError("if_exists debe ser 'append' o 'replace'")

        df = df.copy()
        df.columns = df.columns.str.strip()

        with sqlite3.connect(self.DB_PATH) as con:
            con.execute("PRAGMA foreign_keys = ON;")

            if validate_columns and if_exists == "append":
                cur = con.execute(f"PRAGMA table_info({table});")
                rows = cur.fetchall()
                # (cid, name, type, notnull, dflt_value, pk)
                table_cols = [r[1] for r in rows]

                # Detectar PK autogenerada (rowid alias): pk==1 y type == 'INTEGER'
                # Nota: en SQLite SOLO "INTEGER PRIMARY KEY" se comporta como rowid/autogen.
                pk_autogen_cols = [
                    r[1] for r in rows
                    if r[5] == 1 and (r[2] or "").strip().upper() == "INTEGER"
                ]

                # Columnas esperadas para INSERT: todas menos pk autogenerada
                insertable_cols = [c for c in table_cols if c not in pk_autogen_cols]

                # Recorta DF a columnas insertables (si trae drug_id, lo botamos)
                df = df[[c for c in df.columns if c in insertable_cols]]

                missing = [c for c in insertable_cols if c not in df.columns]
                if missing:
                    raise ValueError(
                        f"Faltan columnas en df para insertar en '{table}': {missing}"
                    )

            if clear_table:
                con.execute(f"DELETE FROM {table};")

            before = con.total_changes

            df.to_sql(
                table,
                con,
                if_exists=if_exists,
                index=False,
                method=method,
                chunksize=chunksize,
            )

            inserted = con.total_changes - before

        return inserted
    
    def pushData(self, df_data, table, clear_table, dictionary):
        """
        Inserta datos en una tabla SQLite desde un DataFrame.
        Si `clear_table` es True, limpia la tabla antes de insertar.
        """
        try:
            if clear_table:
                self.run_query(f"DELETE FROM {table}")
            
            columns_sql = ", ".join([col['column_sql'] for col in dictionary])
            placeholders = ", ".join("?" for _ in dictionary)
            query = f"INSERT INTO {table} ({columns_sql}) VALUES ({placeholders})"
            data = [
                tuple(df_data[col['column_df']].iloc[i] for col in dictionary)
                for i in range(len(df_data))
            ]      
            con = sqlite3.connect(self.DB_PATH)
            cur = con.cursor()
            cur.executemany(query, data)
            con.commit()
            cur.close()
            con.close()
            return len(data)
        except Exception as e:
            print(f"Error: {e}")
            print(traceback.format_exc())
            return 0
        