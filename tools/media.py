"""Growth media definitions for the MIT1002 model.

Each medium is a dict mapping an exchange reaction ID to the bound used for it,
in the form COBRApy expects for ``model.medium``. Import the registry rather
than reading a file::

    from tools.media import MEDIA

    model.medium = MEDIA["minimal_glucose"]

These definitions used to be serialised to ``media_definitions.pkl`` and loaded
with ``pickle.load`` from ~17 places. That was replaced with a plain importable
module because the pickle was a cache of literal dicts -- there was no
computation to save, so it bought nothing while costing three things:

* It was a tracked binary, so a media change could not be reviewed in a diff.
  That matters because media determine which growth phenotypes pass, so an
  unreviewable media change silently changes test outcomes.
* It could go stale. Editing the definitions without regenerating the pickle
  left every consumer silently using the old media.
* ``pickle.load`` is tied to Python and library versions and executes arbitrary
  code on load.

The exchange reactions listed here must exist in the model, so not every
metabolite in the experimental medium can be represented. Where the exchange
reaction is missing the line is commented out rather than deleted, so the
omission stays visible. All exchange reactions use ModelSEED nomenclature.

Primary sources for the recipes are in ``data/media/sources/``. To regenerate
the KBase and CarveMe tables derived from these definitions, run
``python scripts/export_media_tables.py``.
"""


# This file is for defining the different minimal media with no carbon sources
# for the growth tests. These media are defined as dictionaries, where the keys
# are the exchange reactions for the metabolites in the media, and the values
# are the lower bound for the exchange reaction. The lower bound is set to 1000
# for most of the metabolites, other than oxygen.
# The exchange reactions listed here must be present in the model, so not
# every metabolite in the experimental media can actuually be in the media.
# Where the exchange reaction is missing, I've commented out that line.
# All exchange reactions are defined using the standard modelSEED nomenclature

# minimal_media
# The minimal media I initally used for simualtions/gap filling
# Not necessarily based on anything used in the lab
minimal_media = {
    "EX_cpd00067_e0": 1000,  # H+_e0
    "EX_cpd00058_e0": 1000,  # Cu2+_e0
    "EX_cpd00007_e0": 20,  # O2_e0
    "EX_cpd00971_e0": 1000,  # Na+_e0
    "EX_cpd00063_e0": 1000,  # Ca2+_e0
    "EX_cpd00048_e0": 1000,  # Sulfate_e0
    "EX_cpd10516_e0": 1000,  # fe3_e0
    "EX_cpd00254_e0": 1000,  # Mg_e0
    "EX_cpd00009_e0": 1000,  # Phosphate_e0
    "EX_cpd00205_e0": 1000,  # K+_e0
    "EX_cpd00013_e0": 1000,  # NH3_e0
    "EX_cpd00099_e0": 1000,  # Cl-_e0
    "EX_cpd00030_e0": 1000,  # Mn2+_e0
    "EX_cpd00075_e0": 1000,  # Nitrite_e0
    "EX_cpd00001_e0": 1000,  # H2O_e0
    "EX_cpd00034_e0": 1000,  # Zn2+_e0
    "EX_cpd00149_e0": 1000,  # Co2+_e0
}

# Minimal media with glucose and acetate
minimal_glucose = minimal_media.copy()
minimal_glucose["EX_cpd00027_e0"] = 10  # Glucose_e0

minimal_acetate = minimal_media.copy()
minimal_acetate["EX_cpd00029_e0"] = 10  # Acetate_e0

# A version of the minimal media with all of the vitamins in the biomass
# Useful for gap filling to avoid gap filling vitmain biosynthesis pathways
minimal_vitamins = minimal_media.copy()
vitamins = {
    "EX_cpd00010_e0": 100,  # 'CoA [c0]',
    "EX_cpd11493_e0": 100,  # 'ACP [c0]',
    "EX_cpd00015_e0": 100,  # 'FAD [c0]',
    "EX_cpd00016_e0": 100,  # Pyridoxal phosphate [c0] (Vitamin B6)
    "EX_cpd00220_e0": 100,  # Riboflavin [c0] (Vitamin B2)
    "EX_cpd00017_e0": 100,  # 'S-Adenosyl-L-methionine [c0]',
    "EX_cpd00201_e0": 100,  # '10-Formyltetrahydrofolate [c0]',
    "EX_cpd00087_e0": 100,  # Tetrahydrofolate (Folate/Vitamin B9)
    "EX_cpd00345_e0": 100,  # 5-Methyltetrahydrofolate (Active folic acid)
    "EX_cpd00028_e0": 100,  # 'Heme [c0]',
    "EX_cpd00557_e0": 100,  # 'Siroheme [c0]',
    "EX_cpd00264_e0": 100,  # 'Spermidine [c0]',
    "EX_cpd00118_e0": 100,  # 'Putrescine [c0]',
    "EX_cpd00056_e0": 100,  # TPP (Thiamin pyrophosphate/Vitamin B1)
    "EX_cpd15560_e0": 100,  # 'Ubiquinone-8 [c0]',
    "EX_cpd15540_e0": 100,  # 'Phosphatidylglycerol dioctadecanoyl [c0]',
    "EX_cpd15533_e0": 100,  # 'phosphatidylethanolamine dioctadecanoyl [c0]',
    "EX_cpd03736_e0": 100,  # 'Lauroyl-KDO2-lipid IV(A) [c0]',
    "EX_cpd02229_e0": 100,  # 'Bactoprenyl diphosphate [c0]',
    "EX_cpd15665_e0": 100,  # 'Peptidoglycan polymer (n subunits) [c0]',
}
minimal_vitamins.update(vitamins)


# mbm_media
# Minimal Basal Medium used by Zac in the Moran lab for the first round
# of growth tests for A. mac MI1002
# Does not contain an organic carbon source. They added organics to 12 mM
# C-equivalents (e.g. 2 mM glucose = 12 mM C-equivalent)
mbm_media = {
    "EX_cpd00007_e0": 20,  # O2_e0
    # Artificial Sea Water (ASW) Solution
    "EX_cpd00971_e0": 1000,  # Na+_e0 (in NaCl, Na2SO4, NaF, and NaHCO3)
    "EX_cpd00099_e0": 1000,  # Cl-_e0 (in NaCl, MgCl2, CaCl2, and SrCl2)
    "EX_cpd00048_e0": 1000,  # Sulfate (O4S) (in Na2SO4)
    "EX_cpd00205_e0": 1000,  # K+ (in KCl, and KBr)
    "EX_cpd00966_e0": 1000,  # Bromide (Br-) (in KBr)
    "EX_cpd09225_e0": 1000,  # Boric acid (H3BO3) (in H3BO3)
    "EX_cpd00552_e0": 1000,  # Fluoride (F-) (in NaF)
    "EX_cpd00242_e0": 1000,  # Biocarbonate (HCO3-) (in NaHCO3)
    "EX_cpd00254_e0": 1000,  # Mg_e0 (in MgCl2)
    "EX_cpd00063_e0": 1000,  # Ca2+_e0 (in CaCl2)
    "EX_cpd09695_e0": 1000,  # Strontium (Sr2+) (in SrCl2)
    # FeEDTA (Is a chelating agent)
    # TODO: Check if I chould use uncharge Fe instead of Fe3+
    "EX_cpd10516_e0": 1000,  # fe3_e0 (in FeEDTA)
    "EX_cpd00240_e0": 1000,  # EDTA (in FeEDTA)
    # Basal Medium
    "EX_cpd28238_e0": 1000,  # tris-hydrochloride (tris-HCl)
    "EX_cpd00013_e0": 1000,  # Ammonia (in NH4Cl)
    # Cl (from NH4Cl) already included
    # K (from K2HPO4) already included
    "EX_cpd00009_e0": 1000,  # Phosphate (HO4P) (in K2HPO4)
    "EX_cpd00001_e0": 1000,  # H2O (in H2O)
    "EX_cpd00067_e0": 1000,  # H+_e0 (in H2O)
    # Vitamin Supplement
    # Water already included
    "EX_cpd00104_e0": 1000,  # Biotin (Vitamin H)
    "EX_cpd00393_e0": 1000,  # Folate (Folic acid)
    "EX_cpd00263_e0": 1000,  # Pyridoxine (Pyridoxol)
    "EX_cpd00220_e0": 1000,  # Riboflavin
    "EX_cpd00305_e0": 1000,  # Thiamine
    "EX_cpd00133_e0": 1000,  # Nicotinic acid (Niacinamide)
    "EX_cpd00644_e0": 1000,  # Pantothenic acid (Pantothenate)
    "EX_cpd01826_e0": 1000,  # Cyanocobalamin (Dicopac)
    "EX_cpd00443_e0": 1000,  # p-Aminobenzoic acid (ABEE)
    # Not in the media definition, but needed for growth
    "EX_cpd00058_e0": 1000,  # Cu2+_e0  NOT IN MBM MEDIA
    "EX_cpd00030_e0": 1000,  # Mn2+_e0  NOT IN MBM MEDIA
    "EX_cpd00034_e0": 1000,  # Zn2+_e0  NOT IN MBM MEDIA
    "EX_cpd00149_e0": 1000,  # Co2+_e0  NOT IN MBM MEDIA
}

# l1_media
# L1 Minimal Media
# A general purpose marine medium for growing coastal algae
# An enriched seawater medium, with everything added to filtered natural seawater
# FIXME: Does that mean there are other carbon/nitrogen sources in the media?
l1_media = {
    "EX_cpd00007_e0": 20,  # O2_e0
    "EX_cpd00001_e0": 1000,  # H2O
    "EX_cpd00067_e0": 1000,  # H+_e0
    # L1 salts
    "EX_cpd00971_e0": 1000,  # Na+_e0 (in NaNO3, NaH2PO4, NaSiO3, Na2EDTA, NaMoO4, Na3VO4)
    "EX_cpd00209_e0": 1000,  # NO3- (in NaNO3)
    "EX_cpd00009_e0": 1000,  # Phosphate (HO4P) (in NaH2PO4)
    "EX_cpd20826_e0": 1000,  # Silica (O2Si) (in NaSiO3)
    # Trace element solution
    "EX_cpd00240_e0": 1000,  # EDTA (in Na2EDTA)
    "EX_cpd10516_e0": 1000,  # fe3_e0 (in FeCl3)
    "EX_cpd00099_e0": 1000,  # Cl- (in FeCl3, MnCl2, CoCl2)
    "EX_cpd00030_e0": 1000,  # Mn2+ (in MnCl2)
    "EX_cpd00034_e0": 1000,  # Zn2+ (in ZnSO4)
    "EX_cpd00048_e0": 1000,  # Sulfate (O4S) (in ZnSO4, CuSO4, NiSO4)
    "EX_cpd00149_e0": 1000,  # Co2+ (in CoCl2)
    "EX_cpd00058_e0": 1000,  # Cu2+_e0 (in CuSO4)
    "EX_cpd11574_e0": 1000,  # Molybdate (MoO4) (in NaMoO4)
    "EX_cpd03387_e0": 1000,  # Selenite (O3Se) (in H2SeO3)
    "EX_cpd00244_e0": 1000,  # Ni2+ (in NiSO4)
    "EX_cpd08438_e0": 1000,  # Ortho-vanadate (H2O4V) (in Na3VO4)
    "EX_cpd00205_e0": 1000,  # K+ (in K2CrO4)
    "EX_cpd11595_e0": 1000,  # Chromate (H2CrO4) (in K2CrO4)
    # Vitamin solution
    "EX_cpd00305_e0": 1000,  # Thiamine HCl (Vitamin B1)
    "EX_cpd00104_e0": 1000,  # Biotin (Vitamin H)
    "EX_cpd01826_e0": 1000,  # Cyanocobalamin (Vitamin B12)
    # Not in L1, but needed to grow
    # "EX_cpd00013_e0": 1000,  # NH3_e0 (Needed for growth on glucose and acetate)
    "EX_cpd00063_e0": 1000,  # Ca2+_e0 (Needed for growth on alanine)
    "EX_cpd00254_e0": 1000,  # Mg_e0 (Needed for growth on alanine)
}

bashir_c_free = {
    "EX_cpd00007_e0": 20,  # O2_e0
    "EX_cpd00001_e0": 1000,  # H2O
    "EX_cpd00067_e0": 1000,  # H+_e0
    ##################
    # Nitrogen Source
    ##################
    # NH4Cl
    "EX_cpd00013_e0": 1000,  # Ammonia
    "EX_cpd00099_e0": 1000,  # Cl-
    ####################
    # Phosphorus Source
    ####################
    # NaHPo4.H2O
    # Assuming I do not need to add H or H2O
    "EX_cpd00009_e0": 1000,  # Phosphate (HO4P) (in NaH2PO4)
    "EX_cpd00971_e0": 1000,  # Na+_e0
    ########################
    # Trace Metal Mix (TMM)
    ########################
    # ZnSO4*7H2O
    "EX_cpd00034_e0": 1000,  # Zn2+
    "EX_cpd00048_e0": 1000,  # Sulfate (O4S)
    # CoCl2*6H2O
    # Cl- already included in NH4Cl
    "EX_cpd00149_e0": 1000,  # Co2+
    # MnCl2*4H2O
    # Cl- already included in NH4Cl
    "EX_cpd00030_e0": 1000,  # Mn2+
    # Na2MoO4*2H2O
    # Na+ already included in phosphorus source
    "EX_cpd11574_e0": 1000,  # Molybdate (MoO4)
    # Na2SeO3
    # Na+ already included in phosphorus source
    "EX_cpd03387_e0": 1000,  # Selenite (O3Se)
    # NiCl2*6H2O
    # Cl- already included in NH4C
    "EX_cpd00244_e0": 1000,  # Ni2+
    ##################################
    # Artificial Sea Water (ASW) Base
    ##################################
    # NaCl
    # Na+ already included in phosphorus source
    # Cl- already included in NH4Cl
    # KCl
    # Cl- already included in NH4Cl
    "EX_cpd00205_e0": 1000,  # K+
    # CaCl2
    # Cl- already included in NH4Cl
    "EX_cpd00063_e0": 1000,  # Ca2+
    # MgCl2*6H2O
    # Cl- already included in NH4Cl
    "EX_cpd00254_e0": 1000,  # Mg2+
    # MgSO4*7H2O
    # Mg2+ already included in MgCl2
    # Sulfate already included in ZnSO4
}

# Define Franzi's media
# Based on 2024-02-19_Kratzl_Marine_Pro_Medium.xlsx
marine_broth_wo_yeast_and_peptone = {
    # Testing if Cu rescues TPP prucution
    "EX_cpd00058_e0": 1000,  # Cu2+_e0
    # The actual definition
    "EX_cpd00007_e0": 20,  # O2_e0
    "EX_cpd00001_e0": 1000,  # H2O
    "EX_cpd00067_e0": 1000,  # H+_e0
    ################
    # Salt Solution
    ################
    # NaCl
    "EX_cpd00971_e0": 1000,  # Na+
    "EX_cpd00099_e0": 1000,  # Cl-
    # MgSO4*7H2O
    "EX_cpd00254_e0": 1000,  # Mg2+
    "EX_cpd00048_e0": 1000,  # Sulfate (O4S)
    # MgCl2
    # Mg2+ already included in MgSO4*7H2O
    # Cl- already included in NaCl
    # KCl
    # Cl- already included in NaCl
    "EX_cpd00205_e0": 1000,  # K+
    # CaCl2
    # Cl- already included in NaCl
    "EX_cpd00063_e0": 1000,  # Ca2+
    # Boric acid
    "EX_cpd09225_e0": 1000,  # Boric acid (H3BO3)
    # NaHCO3
    # Na+ already included in NaCl
    "EX_cpd00242_e0": 1000,  # Biocarbonate (HCO3-)
    # Na2PO4
    # Na+ already included in NaCl
    "EX_cpd00009_e0": 1000,  # Phosphate (HO4P)
    # FeCl3*6H2O
    # Cl- already included in NaCl
    "EX_cpd10516_e0": 1000,  # fe3_e0
    ####################
    # Nitrogen Solution
    ####################
    # NH4Cl
    # Cl- already included in NaCl
    "EX_cpd00013_e0": 1000,  # Ammonia
    # KNO3
    # K+ already included in KCl
    "EX_cpd00209_e0": 1000,  # NO3-
    ###########
    # Vitamins
    ###########
    # Thiamine HCl
    # Cl- already included in NaCl
    # Assuming I do not need to add H
    "EX_cpd00305_e0": 1000,  # Thiamine
    # Biotin
    "EX_cpd00104_e0": 1000,  # Biotin
    # B12 (cyanocobalamin)
    "EX_cpd01826_e0": 1000,  # Cyanocobalamin
    # Folic acid
    "EX_cpd00393_e0": 1000,  # Folate
    # PABA
    "EX_cpd00443_e0": 1000,  # p-Aminobenzoic acid
    # Nicotinic acid (niacin)
    "EX_cpd00133_e0": 1000,  # Nicotinic acid
    # Inositol
    "EX_cpd00121_e0": 1000,  # Inositol
    # Ca Pantothanate
    # Ca2+ already included in CaCl2
    "EX_cpd00644_e0": 1000,  # Pantothenic acid
    # Pyridoxine HCl
    # Cl- already included in NaCl
    "EX_cpd00263_e0": 1000,  # Pyridoxine
    #################
    # Trace Elements
    #################
    # ZnSO4*7H2O
    # Sulfate already included in MgSO4
    "EX_cpd00034_e0": 1000,  # Zn2+
    # CoCl2*6H2O
    # Cl- already included in NaCl
    "EX_cpd00149_e0": 1000,  # Co2+
    # MnCl2*4H2O
    # Cl- already included in NaCl
    "EX_cpd00030_e0": 1000,  # Mn2+
    # Na2MoO4*2H2O
    # Na+ already included in NaCl
    "EX_cpd11574_e0": 1000,  # Molybdate (MoO4)
    # Na2SeO3
    # Na+ already included in NaCl
    "EX_cpd03387_e0": 1000,  # Selenite (O3Se)
    # NiCl2*6H2O
    # Cl- already included in NaCl
    "EX_cpd00244_e0": 1000,  # Ni2+
}

# Remove all nitrogen sources from Franzi'e medium for the experiments
# with nitrogen sources
marine_broth_no_n = marine_broth_wo_yeast_and_peptone.copy()
marine_broth_no_n.pop("EX_cpd00013_e0", None)  # Ammonia
marine_broth_no_n.pop("EX_cpd00209_e0", None)  # NO3-

# ProMM
# Media definition from Osnat's KBase narrative
# TODO: Check the media composition, by comparing to the media protocol
promm = {
    "EX_cpd00034_e0": 1000,  # Zn2+
    "EX_cpd03387_e0": 1000,  # Selenite
    "EX_cpd00020_e0": 1000,  # Pyruvate
    "EX_cpd00029_e0": 1000,  # Acetate
    "EX_cpd00007_e0": 20,  # O2
    "EX_cpd00244_e0": 1000,  # Ni2+
    "EX_cpd00971_e0": 1000,  # Na+
    "EX_cpd11574_e0": 1000,  # Molybdate
    "EX_cpd00030_e0": 1000,  # Mn2+
    "EX_cpd00254_e0": 1000,  # Mg
    "EX_cpd00205_e0": 1000,  # K+
    "EX_cpd00100_e0": 1000,  # Glycerol
    "EX_cpd00104_e0": 1000,  # BIOT
    "EX_cpd00048_e0": 1000,  # Sulfate
    "EX_cpd00009_e0": 1000,  # Phosphate
    "EX_cpd00001_e0": 1000,  # H2O
    "EX_cpd00067_e0": 1000,  # H+
    "EX_cpd00121_e0": 1000,  # L-Inositol
    "EX_cpd00133_e0": 1000,  # Nicotinamide
    "EX_cpd10516_e0": 1000,  # Fe+3
    "EX_cpd10515_e0": 1000,  # Fe+2
    "EX_cpd00159_e0": 1000,  # L-Lactate
    "EX_cpd00263_e0": 1000,  # Pyridoxol
    "EX_cpd00149_e0": 1000,  # Co2+
    "EX_cpd00058_e0": 1000,  # Cu2+
    "EX_cpd00099_e0": 1000,  # Cl-
    "EX_cpd00305_e0": 1000,  # Thiamin
    "EX_cpd00063_e0": 1000,  # Ca2+
    "EX_cpd00013_e0": 1000,  # NH3
    "EX_cpd00393_e0": 1000,  # Folate
    "EX_cpd00242_e0": 1000,  # H2CO3  # Is this a problem?
    "EX_cpd00644_e0": 1000,  # PAN
    "EX_cpd03424_e0": 1000,  # Vitamin B12
}

# Make a minimal version of ProMM without the carbon sources
promm_no_c = promm.copy()
promm_no_c.pop("EX_cpd00020_e0")  # Pyruvate
promm_no_c.pop("EX_cpd00029_e0")  # Acetate
promm_no_c.pop("EX_cpd00100_e0")  # Glycerol
promm_no_c.pop("EX_cpd00159_e0")  # L-Lactate

# Seawater minimal medium (SWM)- for Koch 2020 phenotypes
swm = {
    "EX_cpd00007_e0": 20,  # O2_e0
    "EX_cpd00067_e0": 1000,  # H+_e0
    "EX_cpd00001_e0": 1000,  # H2Omp (membrane purified)
    # 4.0 g NaSO4
    "EX_cpd00971_e0": 1000,  # Na+_e0
    "EX_cpd00048_e0": 1000,  # Sulfate (SO4)
    # 0.2 g KH2PO4
    "EX_cpd00205_e0": 1000,  # K+
    "EX_cpd00009_e0": 1000,  # Phosphate
    # 0.25 g NH4Cl
    "EX_cpd00013_e0": 1000,  # Ammonia
    "EX_cpd00099_e0": 1000,  # Cl-_e0
    # 20.0 g NaCl
    # Already added Na as part of NaSO4
    # Already added Cl as part of NH4Cl
    # 3.0 g MgCl2×6H2O
    "EX_cpd00254_e0": 1000,  # Mg_e0
    # Already added Cl as part of NH4Cl
    # 0.5 g KCl
    # Already added K+ as part of KH2PO4
    # Already added Cl as part of NH4Cl
    # 0.15 g CaCl2×2H2O
    "EX_cpd00063_e0": 1000,  # Ca2+_e0
    # Already added Cl as part of NH4Cl
    # 0.19 g NaHCO3
    # Already added Na as part of NaSO4
    # TODO: Should I add CO2?
    # 2.1 g FeSO4×7H2O (assuming this is iron(II))
    # FIXME: It required Fe+3 to grow, so I am incluiding Fe+3 in addition to Fe+2, but I am not sure if this is correct
    "EX_cpd10515_e0": 1000,  # Fe+2
    "EX_cpd10516_e0": 1000,  # Fe+3
    # Already added SO4 as part of NaSO4
    # 13.0 mL 25% HCl
    # Already added Cl as part of NH4Cl
    # 5.2 g Titriplex-(III) (Na2-EDTA)
    # Already added Na as part of NaSO4
    "EX_cpd00240_e0": 1000,  # EDTA (in FeEDTA)
    # 30.0 mg H3BO3
    # TODO: Should I add anything?
    # 100.0 mg MnCl2×4H2O
    "EX_cpd00030_e0": 1000,  # Mn2+_e0
    # 190.0 mg CoCl2×6H2O
    "EX_cpd00149_e0": 1000,  # Co2+_e0
    # 24.0 mg NiCl2×6H2O
    "EX_cpd00244_e0": 100,  # Ni2+
    # 2.0 mg CuCl2×2H2O
    "EX_cpd00058_e0": 1000,  # Cu2+_e0
    # 144.0 mg ZnSO4×7H2O
    "EX_cpd00034_e0": 100,  # Zn2+
    # 36.0 mg Na2MoO4×2H2O
    "EX_cpd11574_e0": 100,  # Molybdate
}

# HMB
# Media definition from Osnat's KBase narrative
# TODO: Check the media composition, by comparing to the media protocol
hmb = {
    "EX_cpd00034_e0": 100,  # Zn2+
    "EX_cpd03387_e0": 100,  # Selenite
    "EX_cpd00007_e0": 100,  # O2
    "EX_cpd00244_e0": 100,  # Ni2+
    "EX_cpd00971_e0": 100,  # Na+
    "EX_cpd11574_e0": 100,  # Molybdate
    "EX_cpd00030_e0": 100,  # Mn2+
    "EX_cpd00254_e0": 100,  # Mg
    "EX_cpd00205_e0": 100,  # K+
    "EX_cpd00104_e0": 0.1,  # BIOT
    "EX_cpd00048_e0": 100,  # Sulfate
    "EX_cpd00009_e0": 100,  # Phosphate
    "EX_cpd00001_e0": 100,  # H2O
    "EX_cpd00067_e0": 100,  # H+
    "EX_cpd00121_e0": 0.1,  # L-Inositol
    "EX_cpd00133_e0": 0.1,  # Nicotinamide
    "EX_cpd10516_e0": 100,  # Fe+3
    "EX_cpd10515_e0": 100,  # Fe+2
    "EX_cpd00209_e0": 10,  # Nitrate
    "EX_cpd00263_e0": 0.1,  # Pyridoxol
    "EX_cpd00149_e0": 100,  # Co2+
    "EX_cpd00058_e0": 100,  # Cu2+
    "EX_cpd00099_e0": 100,  # Cl-
    "EX_cpd00305_e0": 0.1,  # Thiamin
    "EX_cpd00063_e0": 100,  # Ca2+
    "EX_cpd00013_e0": 10,  # NH3
    "EX_cpd00393_e0": 0.1,  # Folate
    "EX_cpd00242_e0": 100,  # H2CO3
    "EX_cpd00644_e0": 0.1,  # PAN
    "EX_cpd03424_e0": 0.1,  # Vitamin B12
}

# MMB
# Media definition from Osnat's KBase narrative
# TODO: Check the media composition, by comparing to the media protocol
mmb = {
    "EX_cpd00034_e0": 100,  # Zn2+
    "EX_cpd03387_e0": 100,  # Selenite
    "EX_cpd00007_e0": 100,  # O2
    "EX_cpd00244_e0": 100,  # Ni2+
    "EX_cpd00971_e0": 100,  # Na+
    "EX_cpd11574_e0": 100,  # Molybdate
    "EX_cpd00030_e0": 100,  # Mn2+
    "EX_cpd00254_e0": 100,  # Mg
    "EX_cpd00205_e0": 100,  # K+
    "EX_cpd00048_e0": 100,  # Sulfate
    "EX_cpd00009_e0": 100,  # Phosphate
    "EX_cpd00001_e0": 100,  # H2O
    "EX_cpd00067_e0": 100,  # H+
    "EX_cpd10516_e0": 100,  # Fe+3
    "EX_cpd10515_e0": 100,  # Fe+2
    "EX_cpd00209_e0": 10,  # Nitrate
    "EX_cpd00149_e0": 100,  # Co2+
    "EX_cpd00058_e0": 100,  # Cu2+
    "EX_cpd00099_e0": 100,  # Cl-
    "EX_cpd00063_e0": 100,  # Ca2+
    "EX_cpd00013_e0": 10,  # NH3
    "EX_cpd00242_e0": 100,  # H2CO3
}

# PRO99
# Media definition from Osnat's KBase narrative
# TODO: Check the media composition, by comparing to the media protocol
pro99 = {
    "EX_cpd00034_e0": 1000,  # Zn2+
    "EX_cpd03387_e0": 1000,  # Selenite
    "EX_cpd00007_e0": 1000,  # O2
    "EX_cpd00244_e0": 1000,  # Ni2+
    "EX_cpd00971_e0": 1000,  # Na+
    "EX_cpd11574_e0": 1000,  # Molybdate
    "EX_cpd00030_e0": 1000,  # Mn2+
    "EX_cpd00254_e0": 1000,  # Mg
    "EX_cpd00205_e0": 1000,  # K+
    "EX_cpd00048_e0": 1000,  # Sulfate
    "EX_cpd00009_e0": 1000,  # Phosphate
    "EX_cpd00001_e0": 1000,  # H2O
    "EX_cpd00067_e0": 1000,  # H+
    "EX_cpd10516_e0": 1000,  # Fe+3
    "EX_cpd10515_e0": 1000,  # Fe+2
    "EX_cpd00149_e0": 1000,  # Co2+
    "EX_cpd00058_e0": 1000,  # Cu2+
    "EX_cpd00099_e0": 1000,  # Cl-
    "EX_cpd00063_e0": 1000,  # Ca2+
    "EX_cpd00013_e0": 1000,  # NH3
    "EX_cpd00242_e0": 1000,  # H2CO3
}

lb = {
    "EX_cpd00001_e0": 100,  # H2O
    "EX_cpd00007_e0": 100,  # O2
    "EX_cpd00009_e0": 100,  # Phosphate
    "EX_cpd00018_e0": 100,  # AMP
    "EX_cpd00023_e0": 100,  # L-Glutamate
    "EX_cpd00027_e0": 100,  # D-Glucose
    "EX_cpd00028_e0": 100,  # Heme
    "EX_cpd00030_e0": 100,  # Mn2+
    "EX_cpd00033_e0": 100,  # Glycine
    "EX_cpd00034_e0": 100,  # Zn2+
    "EX_cpd00035_e0": 100,  # L-Alanine
    "EX_cpd00039_e0": 100,  # L-Lysine
    "EX_cpd00041_e0": 100,  # L-Aspartate
    "EX_cpd00046_e0": 100,  # CMP
    "EX_cpd00048_e0": 100,  # Sulfate
    "EX_cpd00051_e0": 100,  # L-Arginine
    "EX_cpd00054_e0": 100,  # L-Serine
    "EX_cpd00058_e0": 100,  # Cu2+
    "EX_cpd00060_e0": 100,  # L-Methionine
    "EX_cpd00063_e0": 100,  # Ca2+
    "EX_cpd00065_e0": 100,  # L-Tryptophan
    "EX_cpd00066_e0": 100,  # L-Phenylalanine
    "EX_cpd00067_e0": 100,  # H+
    "EX_cpd00069_e0": 100,  # L-Tyrosine
    "EX_cpd00084_e0": 100,  # L-Cysteine
    "EX_cpd00091_e0": 100,  # UMP
    "EX_cpd00092_e0": 100,  # Uracil
    "EX_cpd00099_e0": 100,  # Cl-
    "EX_cpd00107_e0": 100,  # L-Leucine
    "EX_cpd00119_e0": 100,  # L-Histidine
    "EX_cpd00126_e0": 100,  # GMP
    "EX_cpd00129_e0": 100,  # L-Proline
    "EX_cpd00149_e0": 100,  # Co2+
    "EX_cpd00156_e0": 100,  # L-Valine
    "EX_cpd00161_e0": 100,  # L-Threonine
    "EX_cpd00182_e0": 100,  # Adenosine
    "EX_cpd00184_e0": 100,  # Thymidine
    "EX_cpd00205_e0": 100,  # K+
    "EX_cpd00215_e0": 100,  # Pyridoxal
    "EX_cpd00218_e0": 100,  # Niacin
    "EX_cpd00219_e0": 100,  # Prephenate
    "EX_cpd00220_e0": 100,  # Riboflavin
    "EX_cpd00226_e0": 100,  # HYXN
    "EX_cpd00239_e0": 100,  # H2S
    "EX_cpd00246_e0": 100,  # Inosine
    "EX_cpd00249_e0": 100,  # Uridine
    "EX_cpd00254_e0": 100,  # Mg
    "EX_cpd00311_e0": 100,  # Guanosine
    "EX_cpd00322_e0": 100,  # L-Isoleucine
    "EX_cpd00381_e0": 100,  # L-Cystine
    "EX_cpd00383_e0": 100,  # Shikimate
    "EX_cpd00393_e0": 100,  # Folate
    "EX_cpd00438_e0": 100,  # Deoxyadenosine
    "EX_cpd00531_e0": 100,  # Hg2+
    "EX_cpd00541_e0": 100,  # Lipoate
    "EX_cpd00644_e0": 100,  # PAN
    "EX_cpd00654_e0": 100,  # Deoxycytidine
    "EX_cpd00793_e0": 100,  # Thiamine phosphate
    "EX_cpd00971_e0": 100,  # Na+
    "EX_cpd01012_e0": 100,  # Cd2+
    "EX_cpd01048_e0": 100,  # Arssenate
    "EX_cpd03424_e0": 100,  # Vitamin B12
    "EX_cpd10515_e0": 100,  # Fe2+
    "EX_cpd10516_e0": 100,  # Fe3+
    "EX_cpd11595_e0": 100,  # Chromate
}

#: All media, keyed by the name used in known_growth_phenotypes.tsv and by
#: every consumer that used to read media_definitions.pkl.
MEDIA = {
    "minimal": minimal_media,
    "minimal_glucose": minimal_glucose,
    "minimal_acetate": minimal_acetate,
    "minimal_vitamins": minimal_vitamins,
    "mbm": mbm_media,
    "l1": l1_media,
    "bashir_c_free": bashir_c_free,
    "marine_broth_wo_yeast_and_peptone": marine_broth_wo_yeast_and_peptone,
    "marine_broth_wo_yeast_and_peptone_no_n": marine_broth_no_n,
    "promm": promm,
    "promm_no_c": promm_no_c,
    "hmb": hmb,
    "mmb": mmb,
    "pro99": pro99,
    "lb": lb,
    "swm": swm,
}

# --------------------------------------------------------------------------
# ModelSEED compound database
# --------------------------------------------------------------------------
#
# The scripts that derive the KBase and CarveMe tables from the media above need
# compound names, formulae and BiGG aliases, which come from a local clone of the
# ModelSEED database. That path used to be hardcoded to one developer's home
# directory in two separate files, so those scripts only ran on one machine.

import json  # noqa: E402  (kept below the definitions, which need no imports)
import os  # noqa: E402

#: Environment variable giving the path to ModelSEED's ``compounds.json``.
MODELSEED_COMPOUNDS_ENV = "MODELSEED_COMPOUNDS"

#: Where the ModelSEED database is looked for if the variable is not set.
MODELSEED_COMPOUNDS_DEFAULT = os.path.join(
    os.path.expanduser("~"),
    "Documents", "PhD", "Segre-lab", "ModelSEEDDatabase",
    "Biochemistry", "compounds.json",
)


def load_modelseed_compounds(path: str | None = None) -> dict:
    """Load ModelSEED's ``compounds.json``, keyed by compound ID.

    Looks at ``path``, then ``$MODELSEED_COMPOUNDS``, then
    :data:`MODELSEED_COMPOUNDS_DEFAULT`. Raises with an actionable message if
    none of those exist, rather than a bare ``FileNotFoundError`` on somebody
    else's home directory.
    """
    candidate = path or os.environ.get(MODELSEED_COMPOUNDS_ENV) or MODELSEED_COMPOUNDS_DEFAULT
    if not os.path.exists(candidate):
        raise FileNotFoundError(
            f"ModelSEED compounds.json not found at {candidate!r}. Clone "
            f"https://github.com/ModelSEED/ModelSEEDDatabase and point "
            f"${MODELSEED_COMPOUNDS_ENV} at its "
            f"Biochemistry/compounds.json, e.g.\n"
            f"  export {MODELSEED_COMPOUNDS_ENV}=/path/to/ModelSEEDDatabase/Biochemistry/compounds.json"
        )
    with open(candidate) as handle:
        return {met["id"]: met for met in json.load(handle)}


def convert_aliases_to_dict(alias_string) -> dict:
    """Parse a ModelSEED ``aliases`` list into ``{database: [ids]}``."""
    return {
        alias.split(":")[0]: [ak.strip() for ak in alias.split(":")[1].split(";")]
        for alias in alias_string
        if alias
    }


def modelseed_id_from_exchange(exchange_reaction: str) -> str:
    """``EX_cpd00027_e0`` -> ``cpd00027``."""
    return exchange_reaction.replace("EX_", "").replace("_e0", "")
