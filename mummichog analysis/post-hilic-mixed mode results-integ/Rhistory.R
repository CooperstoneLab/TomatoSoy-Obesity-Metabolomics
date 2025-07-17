# PID of current job: 2588507
mSet<-InitDataObjects("mass_all", "mummichog", FALSE)
mSet<-SetPeakFormat(mSet, "rmp")
mSet<-UpdateInstrumentParameters(mSet, 10.0, "mixed", "yes", 0.02);
mSet<-Read.PeakListData(mSet, "Replacing_with_your_file_path");
mSet<-SanityCheckMummichogData(mSet)
add.vec <- c("M [1+]","M+H [1+]","M+2H [2+]","M+3H [3+]","M+Na [1+]","M+H+Na [2+]","M+H2O+H [1+]","M-H2O+H [1+]","M-H4O2+H [1+]","M-NH3+H [1+]","M-CO+H [1+]","M-CO2+H [1+]","M-HCOOH+H [1+]","M+HCOONa [1+]","M-HCOONa+H [1+]","M+NaCl [1+]","M-C3H4O2+H [1+]","M-H [1-]","M-2H [2-]","M-H2O-H [1-]","M-H+O [1-]","M+Na-2H [1- ]","M+Cl [1-]","M+ACN-H [1-]","M+HCOO [1-]","M+CH3COO [1-]")
mSet<-Setup.AdductData(mSet, add.vec);
mSet<-PerformAdductMapping(mSet, "mixed")
mSet<-SetPeakEnrichMethod(mSet, "integ", "v2")
mSet<-SetMummichogPval(mSet, 0.05)
mSet<-PerformPSEA(mSet, "hsa_mfn", "current", 5 , 100)
mSet<-PlotPSEAIntegPaths(mSet, "integ_peaks_0_", "png", 72, width=NA)
mSet<-SaveTransformedData(mSet)
