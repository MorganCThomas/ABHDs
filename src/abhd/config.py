"""Single place for paths and analysis parameters. Change cutoffs here, then `make` re-runs what depends on them."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]

# --- paths -------------------------------------------------------------------
RAW = ROOT / "data" / "raw"
DERIVED = ROOT / "data" / "derived"
FIGURES = ROOT / "figures"

ENTRIES_JSON = RAW / "uniprot" / "abhd_entries.json"
SEQUENCES_FASTA = RAW / "uniprot" / "abhd_sequences.fasta"
STRUCTURES = RAW / "structures"
FUNFAM_DIR = RAW / "christina"                      # read-only input from Christina
FUNFAM_META = FUNFAM_DIR / "funfam_metadata.csv"
FUNFAM_HH = FUNFAM_DIR / "human_abhd_funfam_hhalign_matrix.csv"   # HHalign E-values between FunFam HMMs

FF_UNIPROT_MAP = DERIVED / "ff_uniprot_map.json"
PROTBERT_DIR = DERIVED / "protbert"
PROTBERT_NPZ = DERIVED / "protbert.npz"
ELBOW_CSV = DERIVED / "nucleophile_elbow.csv"
STRUCT_SUMMARY = DERIVED / "structure_summary.csv"
SOAP_NPZ = DERIVED / "soap.npz"
PROTBERT_COS = DERIVED / "protbert_cosine.csv"
SOAP_COS = DERIVED / "soap_cosine.csv"

# --- UniProt query (stage 1) -------------------------------------------------
UNIPROT_QUERY = 'family:"ab hydrolase superfamily" AND organism_id:9606'

# --- nucleophile elbow (stage 4) ---------------------------------------------
# Elbow = nucleophile at the beta5 -> alphaC turn: strained backbone (phi>0, psi~-120),
# preceded by a beta-strand residue and followed by helix. Calibrated on known nucleophiles
# (ABHD2 S207, ABHD14B S111, LPL S159, MGLL S122: phi 52-55, psi -117..-126).
NUC_RESIDUES = {"SER": ["OG"], "CYS": ["SG"], "THR": ["OG1"], "ASP": ["OD1", "OD2"]}
ELBOW_PHI = (30, 90)
ELBOW_PSI = (-170, -80)
BETA_PHI_MAX = -50                 # residue i-1: phi in [-180,-50] and (psi >= 80 or psi <= -160)
BETA_PSI_MIN = 80
BETA_PSI_WRAP = -160
HELIX_PHI = (-110, -30)            # residue i+1: helical
HELIX_PSI = (-80, 10)
PLDDT_MIN = 70.0                   # AlphaFold confidence at the nucleophile CA

# Catalytic triad geometry (heavy-atom distances, Angstrom). Deliberately permissive for AF models.
D_NUC_HIS = 4.5                    # nucleophile O/S -> His ring N
D_HIS_ACID = 4.0                   # other His ring N -> Asp/Glu carboxylate O

# The motif regex used in the original notebook (kept for comparison only).
LEGACY_REGEX = r"^(?:G.|.G)[SCT](?:G.|.G)|(?:.P)[S](?:..)$"

# --- SOAP (stage 5) -----------------------------------------------------------
SOAP_RCUT = 8.0
SOAP_NMAX = 8
SOAP_LMAX = 8

# --- ProtBert (stage 3) -------------------------------------------------------
PROTBERT_MODEL = "Rostlab/prot_bert"
PROTBERT_MAX_LEN = 512             # tokens; longer sequences are truncated
PROTBERT_BATCH = 8
