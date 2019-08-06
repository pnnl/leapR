#'clinical_data_table_1
#'
#'
#'

clinical_data_table_1 = function(patientlist, clinical_data) {
  # we presuppose certain column names
  a = mean(clinical_data[patientlist,"age_at_initial_pathologic_diagnosis"])
  b = sd(clinical_data[patientlist,"age_at_initial_pathologic_diagnosis"])
  c = min(clinical_data[patientlist,"age_at_initial_pathologic_diagnosis"])
  d = max(clinical_data[patientlist,"age_at_initial_pathologic_diagnosis"])
  e = sum(clinical_data[patientlist,"PlatinumStatus"]=="Sensitive")
  f = sum(clinical_data[patientlist,"PlatinumStatus"]=="Resistant")

  g = sum(clinical_data[patientlist,"tumor_stage"] %in% c("IIA", "IIB", "IIC", "IID", "II"))
  h = sum(clinical_data[patientlist,"tumor_stage"] %in% c("IIIA", "IIIB", "IIIC", "IIID", "III"))
  i = sum(clinical_data[patientlist,"tumor_stage"] %in% c("IVA", "IVB", "IVC", "IVD", "IV"))

  j = sum(clinical_data[patientlist,"tumor_grade"] == "G2")
  k = sum(clinical_data[patientlist,"tumor_grade"] == "G3")

  l = sum(clinical_data[patientlist,"tumor_residual_disease"] %in% c("1-10 mm", "No Macroscopic disease"))
  m = length(patientlist)-l

  n = sum(clinical_data[patientlist,"vital_status"]=="DECEASED")
  o = length(patientlist)-n

  p = median(clinical_data[patientlist,"days.to.death.or.last_followup"])
  q = sd(clinical_data[patientlist,"days.to.death.or.last_followup"])
  r = median(as.numeric(clinical_data[patientlist,"days_to_tumor_prog_or_recur"]), na.rm=T)
  s = sd(as.numeric(clinical_data[patientlist,"days_to_tumor_prog_or_recur"]), na.rm=T)

  return(c(paste(format(a, digits=1), " (", format(b, digits=3),")", sep=""),paste(c,d,sep="-"),
           "",
           paste(g, " (", format(g*100/length(patientlist), digits=2), "%)", sep=""),
           paste(h, " (", format(h*100/length(patientlist), digits=3), "%)", sep=""),
           paste(i, " (", format(i*100/length(patientlist), digits=3), "%)", sep=""),
           "",
           paste(j, " (", format(j*100/length(patientlist), digits=3), "%)", sep=""),
           paste(k, " (", format(k*100/length(patientlist), digits=3), "%)", sep=""),
           "",
           paste(l, " (", format(l*100/length(patientlist), digits=3), "%)", sep=""),
           paste(m, " (", format(m*100/length(patientlist), digits=3), "%)", sep=""),
           "",
           paste(e, " (", format(e*100/length(patientlist), digits=3), "%)", sep=""),
           paste(f, " (", format(f*100/length(patientlist), digits=3), "%)", sep=""),
           "",
           paste(o, " (", format(o*100/length(patientlist), digits=3), "%)", sep=""),
           paste(n, " (", format(n*100/length(patientlist), digits=3), "%)", sep=""),
           paste(r, " (", format(s, digits=4), ")", sep=""),
           paste(p, " (", format(q, digits=4), ")", sep="")
  ))
}
