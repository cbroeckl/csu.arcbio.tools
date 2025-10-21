
## bioactivity report demonstration

cids <- c(1983, 2519, 445154, 5793)
pc <- csu.pmf.tools::rc.cmpd.get.pubchem(cmpd.cid = cids, get.bioassays = TRUE)
bioactive <- pc$bioassays[which(pc$bioassays$Activity.Outcome == "Active"),]
bioactive <- bioactive[,c("cid", "Target.Accession", "Target.GeneID", "Activity.Value..uM.", "Activity.Name", "Assay.Name", "PubMed.ID")]
head(bioactive)
props <- pc$properties[,c("cid","Title", "MolecularFormula", "SMILES", "XLogP", "MonoisotopicMass", "TPSA")]
out <- merge(props, bioactive, by = 'cid', all.x = TRUE, all.y = TRUE)
dim(out)
head(out)

writexl::write_xlsx(out, "bioactivity.example.xlsx")
