# Genera data/admet.rds con las funciones de admetshiny.
# Ejecutar UNA vez en local (requiere Java JDK para rcdk e internet para PubChem).
# El dashboard luego solo lee el .rds, así no necesita Java al renderizar.

library(admetshiny)

# ---- OPCIÓN A: desde nombres de compuestos (PubChem -> CDK) ----
nombres <- c("aspirin", "ibuprofen", "paracetamol", "caffeine", "metformin",
             "atorvastatin", "diazepam", "morphine", "warfarin", "omeprazole",
             "naproxen", "diclofenac", "propranolol", "loratadine", "fluoxetine",
             "sertraline", "amoxicillin", "ciprofloxacin", "digoxin", "verapamil")

ids   <- getSmilesFromIdentifiers(nombres, from = "name")
ids   <- ids[!is.na(ids$CanonicalSMILES), ]
desc  <- calcCDKDescriptors(ids$CanonicalSMILES)
df    <- mapCDKDescriptors(desc)   # añade violaciones + GI/BBB/P-gp

df$Name <- ids$query[match(df$SMILES, ids$CanonicalSMILES)]

# ---- OPCIÓN B: desde tu CSV/Excel propio ----
# raw <- read.csv("mis_datos.csv", check.names = FALSE)
# mapping <- c("CanonicalSMILES" = "SMILES", "MW" = "MW",
#              "iLOGP" = "LogP", "TPSA" = "TPSA")
# df <- mapADMETColumns(raw, mapping, calculate_cdk = TRUE)

dir.create("data", showWarnings = FALSE)
saveRDS(df, "data/admet.rds")
