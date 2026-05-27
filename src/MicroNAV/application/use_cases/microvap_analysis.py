"""
MicroVAP Lung Microbiome 16S Analysis
Migrated from R (Dickson Lab / Gisell Bustos Moya) to Python

Estructura de carpetas esperada:
    <project_root>/
        data/
            raw/
                microvap_colombia_mothur.shared
                microvap_colombia_mothur.cons.taxonomy
                microvap_colombia_mothur.metadata_completa.txt
                ColombiaShip_Combined.xlsx          ← ddPCR (opcional)

Uso mínimo:
    analysis = MicroVAPAnalysis()
    results  = analysis.run_all()

Uso con rutas explícitas (override):
    analysis = MicroVAPAnalysis(
        shared_file    = "otra/ruta/archivo.shared",
        taxonomy_file  = "otra/ruta/archivo.taxonomy",
        metadata_file  = "otra/ruta/metadata.txt",
        ddpcr_spec_file= "otra/ruta/ddpcr.xlsx",   # opcional
    )
    results = analysis.run_all()
"""

import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier

warnings.filterwarnings("ignore")


# ─────────────────────────────────────────────────────────────────────────────
# Funciones auxiliares (sin cambios respecto a versión anterior)
# ─────────────────────────────────────────────────────────────────────────────

def _bray_curtis_matrix(mat: np.ndarray) -> np.ndarray:
    """Matriz simétrica de disimilitud Bray-Curtis (n × n)."""
    return squareform(pdist(mat, metric="braycurtis"))


def _hellinger(mat: np.ndarray) -> np.ndarray:
    """Transformación Hellinger: sqrt(x / row_sum)."""
    row_sums = np.where(mat.sum(axis=1, keepdims=True) == 0, 1,
                        mat.sum(axis=1, keepdims=True))
    return np.sqrt(mat / row_sums)


def _shannon(mat: np.ndarray) -> np.ndarray:
    """Índice de diversidad de Shannon por muestra (fila)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        p    = np.where(mat > 0, mat / mat.sum(axis=1, keepdims=True), 0)
        logp = np.where(p > 0, np.log(p), 0)
    return -(p * logp).sum(axis=1)


def _simpson(mat: np.ndarray) -> np.ndarray:
    """Índice de Simpson D por muestra (fila)."""
    p = mat / mat.sum(axis=1, keepdims=True)
    return (p ** 2).sum(axis=1)


def _chao1(mat: np.ndarray) -> np.ndarray:
    """Estimador de riqueza Chao1 por muestra (fila)."""
    chao = []
    for row in mat:
        s_obs = (row > 0).sum()
        f1    = (row == 1).sum()
        f2    = (row == 2).sum()
        chao.append(s_obs + (f1 ** 2) / (2 * f2) if f2 > 0
                    else s_obs + f1 * (f1 - 1) / 2)
    return np.array(chao)


def _pcoa(dist_matrix: np.ndarray, k: int = 2):
    """Análisis de Coordenadas Principales (equivalente a cmdscale de R)."""
    n  = dist_matrix.shape[0]
    D2 = dist_matrix ** 2
    J  = np.eye(n) - np.ones((n, n)) / n
    B  = -0.5 * J @ D2 @ J
    eigenvalues, eigenvectors = np.linalg.eigh(B)
    idx          = np.argsort(eigenvalues)[::-1]
    eigenvalues  = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    pos_mask     = eigenvalues > 0
    pct_explained = eigenvalues[pos_mask] / eigenvalues[pos_mask].sum() * 100
    coords = eigenvectors[:, :k] * np.sqrt(np.maximum(eigenvalues[:k], 0))
    return coords, eigenvalues, pct_explained


def _permanova(dist_matrix: np.ndarray, groups: np.ndarray,
               n_permutations: int = 999, seed: int = 42):
    """PERMANOVA simplificado (equivalente a adonis2 de R). Retorna R2 y p-valor."""
    rng           = np.random.default_rng(seed)
    n             = len(groups)
    unique_groups = np.unique(groups)

    def _f_stat(dist_mat, grp):
        ss_total  = (dist_mat ** 2).sum() / (2 * n)
        ss_within = sum(
            (dist_mat[np.ix_(np.where(grp == g)[0],
                              np.where(grp == g)[0])] ** 2).sum()
            / (2 * (grp == g).sum())
            for g in unique_groups
        )
        ss_between  = ss_total - ss_within
        df_between  = len(unique_groups) - 1
        df_within   = n - len(unique_groups)
        return ss_between, ss_total, (ss_between / df_between) / (ss_within / df_within)

    ss_b, ss_t, f_obs = _f_stat(dist_matrix, groups)
    r2    = ss_b / ss_t
    count = sum(
        1 for _ in range(n_permutations)
        if _f_stat(dist_matrix, rng.permutation(groups))[2] >= f_obs
    )
    return r2, (count + 1) / (n_permutations + 1)


# ─────────────────────────────────────────────────────────────────────────────
# Clase principal
# ─────────────────────────────────────────────────────────────────────────────

class MicroVAPAnalysis:
    """
    Pipeline completo de análisis de microbioma 16S para el estudio MicroVAP.

    Parámetros
    ----------
    shared_file      : ruta al archivo .shared de mothur
                       (por defecto: data/raw/microvap_colombia_mothur.shared)
    taxonomy_file    : ruta al archivo .cons.taxonomy
                       (por defecto: data/raw/microvap_colombia_mothur.cons.taxonomy)
    metadata_file    : ruta al metadata .txt (tab-separado)
                       (por defecto: data/raw/microvap_colombia_mothur.metadata_completa.txt)
    ddpcr_spec_file  : ruta al xlsx de ddPCR especímenes — opcional
                       (por defecto: data/raw/ColombiaShip_Combined.xlsx si existe)
    abundance_threshold : umbral mínimo de abundancia relativa (default 0.1 %)
    top_n_otus       : número de OTUs top para análisis downstream
    n_permutations   : permutaciones para PERMANOVA
    n_rf_runs        : número de ejecuciones de Random Forest
    random_seed      : semilla para reproducibilidad
    """

    # Grupos de comparación definidos en el script R original
    SAMPLE_GROUPS = {
        "baseline": ["baseline_no_vap", "baseline_vap"],
        "follow":   ["follow_no_vap",   "follow_vap"],
        "no_vap":   ["baseline_no_vap", "follow_no_vap"],
        "vap":      ["baseline_vap",    "follow_vap"],
    }
    COLORS = {
        "baseline_no_vap": "orange2",
        "baseline_vap":    "skyblue2",
        "follow_no_vap":   "cyan3",
        "follow_vap":      "plum3",
        "NTC":             "rosybrown2",
    }

    # Nombres de archivo predeterminados dentro de data/raw
    _DEFAULT_FILES = {
        "shared_file":     "microvap_colombia_mothur.shared",
        "taxonomy_file":   "microvap_colombia_mothur.cons.taxonomy",
        "metadata_file":   "microvap_colombia_mothur.metadata_completa.txt",
        "ddpcr_spec_file": "ColombiaShip_Combined.xlsx",
    }

    def __init__(
        self,
        shared_file:     str | Path | None = None,
        taxonomy_file:   str | Path | None = None,
        metadata_file:   str | Path | None = None,
        ddpcr_spec_file: str | Path | None = None,
        abundance_threshold: float = 0.1,
        top_n_otus:    int = 20,
        n_permutations: int = 999,
        n_rf_runs:     int = 100,
        random_seed:   int = 42,
    ):
        # ── Resolución de rutas ────────────────────────────────────────────
        self.PROJECT_ROOT = Path.cwd()
        self.SRC_PATH     = self.PROJECT_ROOT / "src"
        self.DATA_RAW     = self.PROJECT_ROOT / "data" / "raw"

        # Agrega src al sys.path si existe (patrón del proyecto MicroNAV)
        if self.SRC_PATH.exists():
            sys.path.append(str(self.SRC_PATH))

        # Resuelve cada archivo: usa el argumento si se pasó,
        # de lo contrario busca el nombre predeterminado en data/raw
        self.shared_file    = self._resolve(shared_file,     "shared_file")
        self.taxonomy_file  = self._resolve(taxonomy_file,   "taxonomy_file")
        self.metadata_file  = self._resolve(metadata_file,   "metadata_file")
        self.ddpcr_spec_file = self._resolve(ddpcr_spec_file, "ddpcr_spec_file",
                                              required=False)

        # ── Parámetros de análisis ─────────────────────────────────────────
        self.abundance_threshold = abundance_threshold
        self.top_n_otus          = top_n_otus
        self.n_permutations      = n_permutations
        self.n_rf_runs           = n_rf_runs
        self.random_seed         = random_seed

        # ── Atributos que se poblan con load_data() ────────────────────────
        self.otu_good:          pd.DataFrame | None = None
        self.otu_good_taxonomy: pd.DataFrame | None = None
        self.meta_df:           pd.DataFrame | None = None
        self.otu_df:            pd.DataFrame | None = None
        self.abs_df:            pd.DataFrame | None = None

    # ─────────────────────────────────────────────────────────────────────
    # Resolución de rutas
    # ─────────────────────────────────────────────────────────────────────

    def _resolve(self, provided: str | Path | None, key: str,
                 required: bool = True) -> Path | None:
        """
        Resuelve una ruta de archivo con la siguiente precedencia:
            1. Argumento explícito del usuario  → se usa tal cual (Path absoluto o relativo a cwd)
            2. Nombre predeterminado en data/raw → self.DATA_RAW / _DEFAULT_FILES[key]
            3. Si required=False y no existe → retorna None (sin lanzar excepción)
        """
        # 1. El usuario proporcionó una ruta
        if provided is not None:
            path = Path(provided)
            if not path.is_absolute():
                path = self.PROJECT_ROOT / path
            if not path.exists():
                raise FileNotFoundError(
                    f"Archivo '{key}' no encontrado en: {path}"
                )
            return path

        # 2. Buscar en data/raw con el nombre predeterminado
        default_path = self.DATA_RAW / self._DEFAULT_FILES[key]
        if default_path.exists():
            return default_path

        # 3. Opcional: retornar None sin error
        if not required:
            return None

        raise FileNotFoundError(
            f"Archivo '{key}' no encontrado.\n"
            f"  Ruta buscada: {default_path}\n"
            f"  Coloca el archivo ahí o pasa la ruta explícitamente al constructor."
        )

    # ─────────────────────────────────────────────────────────────────────
    # 1. CARGA DE DATOS
    # ─────────────────────────────────────────────────────────────────────

    def _clean_specimen(self, s: str) -> str:
        """Limpieza de nombres de muestra equivalente al script R."""
        s = re.sub(r"_S\d+_.*", "", s)
        s = re.sub(r"[_X]+", "_", s)
        return s.rstrip("_")

    def load_otu_data(self) -> "MicroVAPAnalysis":
        """Carga y umbraliza la matriz OTU (sección 2.2 del script R)."""
        encoding = self._detect_encoding(self.shared_file)
        raw      = pd.read_csv(self.shared_file, sep="\t", encoding=encoding)
        otu_trim = raw.iloc[:, 2:].copy()
        otu_trim.index = raw.iloc[:, 1].values

        otu_matrix = otu_trim.values.astype(float)
        row_sums   = otu_matrix.sum(axis=1, keepdims=True)
        otu_pct    = np.where(row_sums > 0, otu_matrix / row_sums * 100, 0)
        otu_matrix[otu_pct < self.abundance_threshold] = 0

        keep_cols  = otu_matrix.sum(axis=0) > 0
        otu_matrix = otu_matrix[:, keep_cols]
        col_names  = otu_trim.columns[keep_cols].tolist()
        row_names  = [self._clean_specimen(str(r)) for r in otu_trim.index]

        self.otu_good = pd.DataFrame(otu_matrix, index=row_names, columns=col_names)
        return self

    def load_taxonomy(self) -> "MicroVAPAnalysis":
        """Carga el archivo de taxonomía (sección 2.2 del script R)."""
        encoding = self._detect_encoding(self.taxonomy_file)
        tax1 = pd.read_csv(self.taxonomy_file, sep="\t", index_col=0, encoding=encoding)
        tax1.columns = ["Size", "Taxonomy"]

        levels = ["Domain", "Phylum", "Class", "Order", "Family", "Genus", "_extra"]
        tax2   = tax1["Taxonomy"].str.split(";", expand=True)
        tax2.columns = levels[: tax2.shape[1]]
        if "_extra" in tax2.columns:
            tax2 = tax2.drop(columns=["_extra"])
        for col in tax2.columns:
            tax2[col] = tax2[col].str.replace(r"\(.*?\)", "", regex=True).str.strip()

        tax_df = pd.concat(
            [pd.DataFrame({"OTU": tax1.index, "Size": tax1["Size"]}).reset_index(drop=True),
             tax2.reset_index(drop=True)],
            axis=1,
        )
        # Usar índice numérico para evitar ambigüedad OTU-como-índice y OTU-como-columna
        tax_df = tax_df.reset_index(drop=True)

        if self.otu_good is not None:
            keep = tax_df["OTU"].isin(self.otu_good.columns)
            self.otu_good_taxonomy = tax_df[keep].reset_index(drop=True).copy()
        else:
            self.otu_good_taxonomy = tax_df.reset_index(drop=True).copy()
        return self

    @staticmethod
    def _detect_encoding(path: Path) -> str:
        """
        Detecta el encoding de un archivo de texto leyendo los primeros bytes.
        Cubre los casos más comunes en archivos generados en Windows:
            UTF-16 LE/BE (BOM 0xFF 0xFE / 0xFE 0xFF)
            UTF-8 con BOM (0xEF 0xBB 0xBF)
            latin-1 / cp1252 (fallback para archivos Windows sin BOM)
        """
        with open(path, "rb") as f:
            bom = f.read(4)
        if bom[:2] in (b"\xff\xfe", b"\xfe\xff"):
            return "utf-16"
        if bom[:3] == b"\xef\xbb\xbf":
            return "utf-8-sig"
        # Intento utf-8 puro; si falla se usa latin-1 (nunca lanza error)
        try:
            with open(path, encoding="utf-8") as f:
                f.read(1024)
            return "utf-8"
        except UnicodeDecodeError:
            return "latin-1"

    def load_metadata(self) -> "MicroVAPAnalysis":
        """Carga el metadata y lo combina con la abundancia relativa (sección 2.3)."""
        encoding = self._detect_encoding(self.metadata_file)
        meta = pd.read_csv(self.metadata_file, sep="\t", encoding=encoding)
        meta["specimen"] = meta["specimen"].astype(str).apply(self._clean_specimen)
        self.meta_df = meta

        otu_rel = self.otu_good.div(self.otu_good.sum(axis=1), axis=0) * 100
        otu_rel.index = [self._clean_specimen(str(i)) for i in otu_rel.index]

        otu_df = (
            meta.set_index("specimen")
            .join(otu_rel, how="inner")
            .reset_index()
            .rename(columns={"index": "specimen"})
        )
        self.otu_df = otu_df
        return self

    def load_ddpcr(self) -> "MicroVAPAnalysis":
        """
        Carga datos de abundancia absoluta ddPCR — solo especímenes (sección 3.1).
        ddpcr_ctrl_file no se utiliza en esta versión del pipeline.
        """
        if self.ddpcr_spec_file is None:
            return self

        spec = pd.read_excel(self.ddpcr_spec_file, sheet_name=2)[
            ["Samples", "Copies/ul DNA"]
        ].rename(columns={"Samples": "sample_name", "Copies/ul DNA": "DNA_copies_per_ul"})

        def _assign_type(name: str) -> str:
            for label in ["follow_no_vap", "baseline_no_vap", "baseline_vap", "follow_vap"]:
                if label in name:
                    return label
            m = re.match(r"^[A-Za-z]+", name)
            return m.group() if m else "unknown"

        spec["Type"] = spec["sample_name"].apply(_assign_type)
        type_order   = ["NTC", "baseline_no_vap", "follow_no_vap", "baseline_vap", "follow_vap"]
        spec["Type"] = pd.Categorical(spec["Type"], categories=type_order, ordered=True)
        self.abs_df  = spec
        return self

    def load_data(self) -> "MicroVAPAnalysis":
        """Carga todos los archivos de entrada en el orden correcto."""
        self.load_otu_data()
        self.load_taxonomy()
        self.load_metadata()
        self.load_ddpcr()
        return self

    # ─────────────────────────────────────────────────────────────────────
    # 2. SUBCONJUNTOS
    # ─────────────────────────────────────────────────────────────────────

    def _subset(self, group_key: str) -> pd.DataFrame:
        """Retorna filas de otu_df para el grupo solicitado."""
        labels = self.SAMPLE_GROUPS[group_key]
        return self.otu_df[self.otu_df["samples"].isin(labels)].copy()

    def _otu_matrix(self, subset_df: pd.DataFrame):
        """Extrae columnas OTU como array numpy (sin varianza cero)."""
        otu_cols = [c for c in subset_df.columns if c.startswith("Otu")]
        mat      = subset_df[otu_cols].values.astype(float)
        keep     = mat.var(axis=0) > 0
        return mat[:, keep], np.array(otu_cols)[keep]

    # ─────────────────────────────────────────────────────────────────────
    # 3. ABUNDANCIA ABSOLUTA (ddPCR)
    # ─────────────────────────────────────────────────────────────────────

    def absolute_abundance_stats(self, group_types: list[str] | None = None) -> pd.DataFrame:
        """Mediana, Q1, Q3 de DNA_copies_per_ul por Type."""
        if self.abs_df is None:
            raise ValueError("Datos ddPCR no cargados. Llama load_ddpcr() primero.")
        df = self.abs_df.copy()
        if group_types:
            df = df[df["Type"].isin(group_types)]
        return (
            df.groupby("Type", observed=True)["DNA_copies_per_ul"]
            .agg(median_DNA_copies="median",
                 Q1=lambda x: x.quantile(0.25),
                 Q3=lambda x: x.quantile(0.75),
                 n="count")
            .reset_index()
        )

    def wilcoxon_absolute(self, type_a: str, type_b: str) -> dict:
        """Test de Wilcoxon entre dos tipos ddPCR."""
        if self.abs_df is None:
            raise ValueError("Datos ddPCR no cargados.")
        a    = self.abs_df.loc[self.abs_df["Type"] == type_a, "DNA_copies_per_ul"]
        b    = self.abs_df.loc[self.abs_df["Type"] == type_b, "DNA_copies_per_ul"]
        stat, p = stats.mannwhitneyu(a, b, alternative="two-sided")
        return {"type_a": type_a, "type_b": type_b, "statistic": stat, "p_value": p}

    def kruskal_absolute(self, group_types: list[str]) -> dict:
        """Test de Kruskal-Wallis entre múltiples tipos ddPCR."""
        if self.abs_df is None:
            raise ValueError("Datos ddPCR no cargados.")
        groups  = [self.abs_df.loc[self.abs_df["Type"] == t, "DNA_copies_per_ul"].values
                   for t in group_types]
        stat, p = stats.kruskal(*groups)
        return {"statistic": stat, "p_value": p, "groups": group_types}

    # ─────────────────────────────────────────────────────────────────────
    # 4. DIVERSIDAD BETA — PCoA
    # ─────────────────────────────────────────────────────────────────────

    def pcoa_braycurtis(self, group_key: str, k: int = 2) -> dict:
        """PCoA con disimilitud Bray-Curtis + PERMANOVA."""
        sub      = self._subset(group_key)
        mat, _   = self._otu_matrix(sub)
        dist_mat = _bray_curtis_matrix(mat)
        coords, eigvals, pct = _pcoa(dist_mat, k=k)
        groups   = sub["samples"].values
        r2, pval = _permanova(dist_mat, groups, self.n_permutations, self.random_seed)

        coord_df = pd.DataFrame(coords, index=sub["specimen"].values,
                                columns=[f"PCoA{i+1}" for i in range(k)])
        coord_df["SampleType"] = groups
        return {"coords": coord_df, "pct_explained": pct, "eigenvalues": eigvals,
                "group_labels": sub["samples"], "permanova": {"R2": r2, "p_value": pval}}

    def pcoa_hellinger(self, group_key: str, k: int = 2) -> dict:
        """PCoA con distancia Hellinger + Euclidiana + PERMANOVA."""
        sub      = self._subset(group_key)
        mat, _   = self._otu_matrix(sub)
        hel      = _hellinger(mat)
        dist_mat = squareform(pdist(hel, metric="euclidean"))
        coords, eigvals, pct = _pcoa(dist_mat, k=k)
        groups   = sub["samples"].values
        r2, pval = _permanova(dist_mat, groups, self.n_permutations, self.random_seed)

        coord_df = pd.DataFrame(coords, index=sub["specimen"].values,
                                columns=[f"PCoA{i+1}" for i in range(k)])
        coord_df["SampleType"] = groups
        return {"coords": coord_df, "pct_explained": pct, "eigenvalues": eigvals,
                "group_labels": sub["samples"], "permanova": {"R2": r2, "p_value": pval}}

    def pca_hellinger(self, group_key: str, n_components: int = 2) -> dict:
        """PCA sobre matriz Hellinger (equivalente a rda de R)."""
        sub    = self._subset(group_key)
        mat, _ = self._otu_matrix(sub)
        hel    = _hellinger(mat)
        pca    = PCA(n_components=n_components)
        coords = pca.fit_transform(hel)
        pct    = pca.explained_variance_ratio_ * 100

        coord_df = pd.DataFrame(coords, index=sub["specimen"].values,
                                columns=[f"PC{i+1}" for i in range(n_components)])
        coord_df["SampleType"] = sub["samples"].values
        return {"coords": coord_df, "pct_explained": pct, "group_labels": sub["samples"]}

    def kruskal_pcoa(self, pcoa_result: dict) -> dict:
        """Kruskal-Wallis sobre PCoA1 y PCoA2 por SampleType."""
        df   = pcoa_result["coords"]
        axes = [c for c in df.columns if c.startswith("PCoA")]
        return {
            ax: dict(zip(["statistic", "p_value"],
                         stats.kruskal(*[g[ax].values for _, g in df.groupby("SampleType")])))
            for ax in axes
        }

    # ─────────────────────────────────────────────────────────────────────
    # 5. DIVERSIDAD ALFA
    # ─────────────────────────────────────────────────────────────────────

    def alpha_diversity(self, group_key: str) -> pd.DataFrame:
        """Shannon, Simpson D, (1-D) y Chao1 por muestra."""
        sub      = self._subset(group_key)
        otu_cols = [c for c in sub.columns if c.startswith("Otu")]
        mat      = sub[otu_cols].values.astype(float)
        simp_d   = _simpson(mat)
        return pd.DataFrame({
            "specimen":          sub["specimen"].values,
            "samples":           sub["samples"].values,
            "Shannon":           _shannon(mat),
            "Simpson_D":         simp_d,
            "Simpson_1_minus_D": 1 - simp_d,
            "Chao1":             _chao1(mat),
            "Total_Reads":       mat.sum(axis=1),
            "Dominance":         mat.max(axis=1),
        })

    def alpha_diversity_stats(self, group_key: str) -> dict:
        """Resumen estadístico y tests de Wilcoxon para cada índice alfa."""
        alpha_df = self.alpha_diversity(group_key)
        indices  = ["Shannon", "Simpson_D", "Simpson_1_minus_D", "Chao1"]
        summary  = alpha_df.groupby("samples")[indices].agg(["mean", "sem"]).reset_index()

        wilcoxon = {}
        labels   = alpha_df["samples"].unique()
        if len(labels) == 2:
            for idx in indices:
                a = alpha_df.loc[alpha_df["samples"] == labels[0], idx].dropna()
                b = alpha_df.loc[alpha_df["samples"] == labels[1], idx].dropna()
                stat, p = stats.mannwhitneyu(a, b, alternative="two-sided")
                wilcoxon[idx] = {"statistic": stat, "p_value": p}

        return {"summary": summary, "wilcoxon": wilcoxon}

    # ─────────────────────────────────────────────────────────────────────
    # 6. ABUNDANCIA RELATIVA
    # ─────────────────────────────────────────────────────────────────────

    def top_otus(self, group_key: str, n: int | None = None) -> list[str]:
        """Nombres de las top-n OTUs más abundantes para un grupo."""
        n        = n or self.top_n_otus
        sub      = self._subset(group_key)
        otu_cols = [c for c in sub.columns if c.startswith("Otu")]
        return sub[otu_cols].mean().nlargest(n).index.tolist()

    def relative_abundance_df(self, group_key: str, n: int | None = None) -> pd.DataFrame:
        """DataFrame en formato largo con taxonomía (Phylum) unida."""
        n         = n or self.top_n_otus
        sub       = self._subset(group_key)
        top       = self.top_otus(group_key, n=n)
        meta_cols = [c for c in sub.columns if not c.startswith("Otu")]
        long      = sub[meta_cols + top].melt(
            id_vars=meta_cols, var_name="OTU", value_name="Percentage"
        )
        if self.otu_good_taxonomy is not None:
            tax_ref = self.otu_good_taxonomy[["OTU", "Phylum"]].reset_index(drop=True)
            long = long.merge(tax_ref, on="OTU", how="left")
        return long

    def relative_abundance_summary(self, group_key: str, by: str = "pneumonia") -> pd.DataFrame:
        """Media ± SEM de abundancia relativa agregada por grupo."""
        long = self.relative_abundance_df(group_key)
        return (
            long.groupby([by, "OTU", "Phylum"])["Percentage"]
            .agg(Mean_perc="mean", SEM=lambda x: np.sqrt(x.var() / len(x)))
            .reset_index()
        )

    def relative_abundance_wilcoxon(self, group_key: str, level: str = "OTU") -> pd.DataFrame:
        """Test de Wilcoxon por OTU o Phylum."""
        long   = self.relative_abundance_df(group_key)
        labels = long["samples"].unique()
        if len(labels) != 2:
            raise ValueError(f"Se esperaban 2 grupos de muestras, se encontraron {len(labels)}.")

        results = []
        for entity, grp in long.groupby(level):
            a = grp.loc[grp["samples"] == labels[0], "Percentage"]
            b = grp.loc[grp["samples"] == labels[1], "Percentage"]
            if len(a) < 2 or len(b) < 2:
                continue
            stat, p = stats.mannwhitneyu(a, b, alternative="two-sided")
            results.append({level: entity, "p_value": p, "statistic": stat})

        return pd.DataFrame(results).sort_values("p_value").reset_index(drop=True)

    # ─────────────────────────────────────────────────────────────────────
    # 7. BRAY-CURTIS INTRA-GRUPO
    # ─────────────────────────────────────────────────────────────────────

    def bray_curtis_within(self) -> pd.DataFrame:
        """Disimilitud Bray-Curtis media dentro de cada grupo de muestras."""
        all_rows = []
        for group_key in self.SAMPLE_GROUPS:
            sub          = self._subset(group_key)
            otu_cols     = [c for c in sub.columns if c.startswith("Otu")]
            dist_mat     = _bray_curtis_matrix(sub[otu_cols].values.astype(float))
            specimens    = sub["specimen"].values
            sample_types = sub["samples"].values

            for i in range(len(specimens)):
                within = [dist_mat[i, j] for j in range(len(specimens))
                          if i != j and sample_types[i] == sample_types[j]]
                if within:
                    all_rows.append({
                        "specimen":  specimens[i],
                        "pneumonia": sub["pneumonia"].values[i],
                        "samples":   sample_types[i],
                        "Mean_BC":   np.mean(within),
                        "SEM_BC":    np.std(within) / np.sqrt(len(within)),
                    })
        return pd.DataFrame(all_rows)

    # ─────────────────────────────────────────────────────────────────────
    # 8. RANDOM FOREST
    # ─────────────────────────────────────────────────────────────────────

    def _family_aggregate(self, subset_df: pd.DataFrame) -> pd.DataFrame:
        """Agrega datos OTU al nivel de Familia (equivalente a fam_df_* en R)."""
        otu_cols = [c for c in subset_df.columns if c.startswith("Otu")]
        long     = subset_df[["specimen"] + otu_cols].melt(
            id_vars="specimen", var_name="OTU", value_name="Percentage"
        )
        if self.otu_good_taxonomy is not None:
            tax_ref = self.otu_good_taxonomy[["OTU", "Family"]].reset_index(drop=True)
            long = long.merge(tax_ref, on="OTU", how="left")
        wide = (
            long.groupby(["specimen", "Family"])["Percentage"].sum()
            .reset_index()
            .pivot(index="specimen", columns="Family", values="Percentage")
            .fillna(0)
        )
        wide.columns.name = None
        meta_cols = [c for c in subset_df.columns if not c.startswith("Otu")]
        return subset_df[meta_cols].set_index("specimen").join(wide).reset_index()

    def random_forest(self, group_key: str, target: str = "pneumonia") -> dict:
        """
        Ejecuta n_rf_runs modelos de Random Forest y retorna la importancia de features.

        Retorna
        -------
        dict con:
            importance_df : DataFrame con Feature y MeanDecreaseAccuracy media
            top_features  : top-20 features por importancia
            oob_errors    : lista de tasas de error OOB por ejecución
        """
        sub    = self._subset(group_key)
        fam_df = self._family_aggregate(sub)

        family_cols = [c for c in fam_df.columns
                       if fam_df[c].dtype in [np.float64, np.int64] and c != target]
        X = fam_df[family_cols].fillna(0)
        y = fam_df[target].astype(str)

        importances_list, oob_errors = [], []
        rng = np.random.default_rng(self.random_seed)

        for _ in range(self.n_rf_runs):
            rf = RandomForestClassifier(
                n_estimators=500, oob_score=True,
                random_state=int(rng.integers(0, 2**31)), n_jobs=-1,
            )
            rf.fit(X, y)
            importances_list.append(pd.Series(rf.feature_importances_, index=family_cols))
            oob_errors.append(1 - rf.oob_score_)

        mean_imp = (
            pd.DataFrame(importances_list).mean()
            .sort_values(ascending=False)
            .reset_index()
        )
        mean_imp.columns = ["Feature", "MeanDecreaseAccuracy"]
        return {"importance_df": mean_imp, "top_features": mean_imp.head(20),
                "oob_errors": oob_errors}

    # ─────────────────────────────────────────────────────────────────────
    # 9. RUN ALL
    # ─────────────────────────────────────────────────────────────────────

    def run_all(self) -> dict:
        """
        Ejecuta el pipeline completo y retorna un dict anidado con todos los resultados.

        Estructura
        ----------
        {
            "otu_good"           : DataFrame,
            "otu_good_taxonomy"  : DataFrame,
            "meta_df"            : DataFrame,
            "otu_df"             : DataFrame,
            "absolute_abundance" : {summary, kruskal, pairwise},   # solo si ddpcr_spec_file existe
            "pcoa_braycurtis"    : {baseline, follow, no_vap, vap},
            "pcoa_hellinger"     : {baseline, follow, no_vap, vap},
            "pca_hellinger"      : {baseline, follow, no_vap, vap},
            "alpha_diversity"    : {baseline, follow, no_vap, vap},
            "rel_abundance"      : {baseline, follow, no_vap, vap},
            "bray_curtis"        : DataFrame,
            "random_forest"      : {baseline, follow, no_vap, vap},
        }
        """
        self.load_data()
        results: dict = {
            "otu_good":          self.otu_good,
            "otu_good_taxonomy": self.otu_good_taxonomy,
            "meta_df":           self.meta_df,
            "otu_df":            self.otu_df,
        }

        # Abundancia absoluta (solo si el archivo ddPCR existe)
        if self.abs_df is not None:
            all_types = ["NTC", "baseline_no_vap", "follow_no_vap", "baseline_vap", "follow_vap"]
            results["absolute_abundance"] = {
                "summary":  self.absolute_abundance_stats(all_types),
                "kruskal":  self.kruskal_absolute(all_types),
                "pairwise": [self.wilcoxon_absolute(a, b) for a, b in [
                    ("baseline_no_vap", "baseline_vap"),
                    ("follow_no_vap",   "follow_vap"),
                    ("baseline_no_vap", "follow_no_vap"),
                    ("baseline_vap",    "follow_vap"),
                ]],
            }

        # Diversidad beta
        for key in ("pcoa_braycurtis", "pcoa_hellinger", "pca_hellinger"):
            results[key] = {}
        for gk in self.SAMPLE_GROUPS:
            results["pcoa_braycurtis"][gk] = self.pcoa_braycurtis(gk)
            results["pcoa_hellinger"][gk]  = self.pcoa_hellinger(gk)
            results["pca_hellinger"][gk]   = self.pca_hellinger(gk)

        # Diversidad alfa
        results["alpha_diversity"] = {
            gk: {"raw": self.alpha_diversity(gk), "stats": self.alpha_diversity_stats(gk)}
            for gk in self.SAMPLE_GROUPS
        }

        # Abundancia relativa
        results["rel_abundance"] = {
            gk: {
                "long":     self.relative_abundance_df(gk),
                "summary":  self.relative_abundance_summary(gk),
                "wilcoxon": self.relative_abundance_wilcoxon(gk),
            }
            for gk in self.SAMPLE_GROUPS
        }

        # Bray-Curtis intra-grupo
        results["bray_curtis"] = self.bray_curtis_within()

        # Random Forest
        results["random_forest"] = {
            gk: self.random_forest(gk, target=target)
            for gk, target in [("baseline", "pneumonia"), ("follow", "pneumonia"),
                                ("no_vap", "samples"),    ("vap", "samples")]
        }

        return results


# ─────────────────────────────────────────────────────────────────────────────
# Función de acceso rápido
# ─────────────────────────────────────────────────────────────────────────────

def build_analysis(
    shared_file:     str | Path | None = None,
    taxonomy_file:   str | Path | None = None,
    metadata_file:   str | Path | None = None,
    ddpcr_spec_file: str | Path | None = None,
    **kwargs,
) -> dict:
    """
    Wrapper de una línea. Retorna el dict completo de run_all().

    Si no se pasan rutas, busca automáticamente en data/raw/.

    Ejemplo
    -------
    results = build_analysis()                        # rutas automáticas
    results = build_analysis(shared_file="mi/ruta")  # override parcial

    df_alpha = results["alpha_diversity"]["baseline"]["raw"]
    df_pcoa  = results["pcoa_braycurtis"]["no_vap"]["coords"]
    """
    return MicroVAPAnalysis(
        shared_file=shared_file, taxonomy_file=taxonomy_file,
        metadata_file=metadata_file, ddpcr_spec_file=ddpcr_spec_file,
        **kwargs,
    ).run_all()
