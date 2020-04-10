
PCR_Neg_Limit <- 35 #From paper
PCR_Pos_Proxy <- 25 #Arbitrary

load_gautret_data <- function() {
  gautret <- readxl::read_excel(here("Gautret_et_al.xlsx"), 
                                sheet = "Patients", na = c("ND","Unknown", "NF")) %>% 
    mutate(Group = if_else(Hydroxychloroquine == "No", "Control", 
                           if_else(Azithromycin == "No", "Hydroxychloroquine", "Hydroxychloroquine + Azithromycin")
                           ),
           GroupShort = if_else(Hydroxychloroquine == "No", "Control", 
                           if_else(Azithromycin == "No", "HCQ", "HCQ + AZ")
           ),
           Days_From_Onset_Imputed = 
             as.integer(if_else(is.na(Days_From_Onset) | Days_From_Onset == "-", "0", Days_From_Onset)))
  
  gautret$Clinical_Status[gautret$Clinical_Status == "LTRI"] <- "LRTI"
  
  gautret %>% mutate(
    patient_label = paste0(ID, " - ", GroupShort),
    Treated = Group != "Control"
  )
  
}