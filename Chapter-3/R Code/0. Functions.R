
# ------------------------------------------------------------------------------
# Life expectancy computation functions
# ------------------------------------------------------------------------------

# Computes life expectancy values from an abridged life table built from age-specific mortality rates (mx).

LE_fun <- function (mx) {
  # 5 year age classes from 65 to 85+
  age = seq(0,90,by=5)
  # n and ax
  n = c(diff(age), 1)
  ax = 0.5 * n
  # probability of death
  qx = (n * mx)/(1 + (n - ax) * mx)
  qx = c(qx[-(length(qx))], 1)
  qx[qx > 1] = 1
  # survival probability
  px = 1 - qx
  # survivors 
  lx = c(100000,rep(0,(length(mx)-1)))
  for (i in 1:(length(mx) -1)){
    lx[i+1] <- lx[i]*px[i] }
  # deaths 
  dx = lx * qx
  # person-years 
  Lx = rep(0,length(mx))
  for (i in 1:length(mx) -1){
    Lx[i] = lx[i+1]*n[i] + ax[i]*dx[i] }
  Lx[length(mx)] = lx[length(mx)]/mx[length(mx)]  
  # life expectancy
  LE = rev(cumsum(rev(Lx)))/lx
  # life table
  LT <- data.frame(ax=ax, mx=mx, qx=qx, dx=dx, lx=lx, Lx=Lx, LE=LE)
  # return life expectancy
  return(LE)
}

ex_demotools <- function(nMx){
  require(DemoTools)
  
  # ex0 <- round((lt_abridged(nMx = nMx, 
  #                           Age = seq(0, 90, by = 5), 
  #                           axmethod = "un",
  #                           mod = FALSE) %>% 
  #                 as.data.frame() %>% 
  #                 pull(ex))[1],5)
  
  
  ex0 <- round((LE_fun(nMx)),5)
}   

ex0_demotools <- function(nMx){
  require(DemoTools)
  
  # ex0 <- round((lt_abridged(nMx = nMx, 
  #                           Age = seq(0, 90, by = 5), 
  #                           axmethod = "un",
  #                           mod = FALSE) %>% 
  #                 as.data.frame() %>% 
  #                 pull(ex))[1],5)
  
  
  ex0 <- round((LE_fun(nMx)),5)[1]
}   

ex65_demotools <- function(nMx){
  require(DemoTools)
  
  # ex0 <- round((lt_abridged(nMx = nMx, 
  #                           Age = seq(0, 90, by = 5), 
  #                           axmethod = "un",
  #                           mod = FALSE) %>% 
  #                 as.data.frame() %>% 
  #                 pull(ex))[1],5)
  
  
  ex0 <- round((LE_fun(nMx)),5)[14]
}   

ex25_demotools <- function(nMx){
  require(DemoTools)
  
  # ex0 <- round((lt_abridged(nMx = nMx, 
  #                           Age = seq(0, 90, by = 5), 
  #                           axmethod = "un",
  #                           mod = FALSE) %>% 
  #                 as.data.frame() %>% 
  #                 pull(ex))[1],5)
  
  
  ex0 <- round((LE_fun(nMx)),5)[6]
}   

# ------------------------------------------------------------------------------
# Utility functions for data manipulation and harmonization
# ------------------------------------------------------------------------------

sex_factor <- function(dt){
  dt <- dt %>% 
    mutate(Sex = case_when(
      Sex == 1 ~ "Male",
      Sex == 2 ~ "Female",
      .default = NA), Sex = factor(Sex, levels=c("Male", "Female")))
}

removedcause_factor <- function(dt){
  dt <-  dt %>%   
    mutate(RemovedCause = factor(RemovedCause,
                                 levels=c("Drugs", 
                                          "Alcohol",
                                          "Suicide",
                                          "Dod")))
}

from_code_to_region <- function(dt){
  dt <- dt %>% 
    mutate(Region = case_when( 
      Place == "01" ~ "Piemonte",
      Place == "02" ~ "Valle d'Aosta",
      Place == "03" ~ "Lombardia", 
      Place == "04" ~ "Trentino Alto-Adige",
      Place == "05" ~ "Veneto",
      Place == "06" ~ "Friuli-Venezia Giulia",
      Place == "07" ~ "Liguria",
      Place == "08" ~ "Emilia-Romagna",
      Place == "09" ~ "Toscana",
      Place == "10" ~ "Umbria",
      Place == "11" ~ "Marche",
      Place == "12" ~ "Lazio",
      Place == "13" ~ "Abruzzo",
      Place == "14" ~ "Molise",
      Place == "15" ~ "Campania",
      Place == "16" ~ "Puglia",
      Place == "17" ~ "Basilicata",
      Place == "18" ~ "Calabria",
      Place == "19" ~ "Sicilia",
      Place == "20" ~ "Sardegna" ,
      .default = NA)) %>% 
    relocate(Region, .before=Place) %>% 
    mutate(Region = factor(Region, 
                           levels=c("Piemonte",
                                    "Valle d'Aosta", "Lombardia", 
                                    "Trentino Alto-Adige", "Veneto", 
                                    "Friuli-Venezia Giulia", "Liguria",
                                    "Emilia-Romagna", "Toscana", "Umbria",
                                    "Marche", "Lazio", "Abruzzo", 
                                    "Molise", "Campania","Puglia",
                                    "Basilicata", "Calabria", "Sicilia",
                                    "Sardegna")))
}

make_nuts1 <- function(dt){
  dt <- dt %>% 
    mutate(Nuts1 = case_when(
      Region %in% c("Piemonte",
                    "Valle d'Aosta",
                    "Lombardia", 
                    "Liguria") ~ "Nord-Ovest",
      Region %in% c("Trentino Alto-Adige",
                    "Veneto",
                    "Friuli-Venezia Giulia",
                    "Emilia-Romagna") ~ "Nord-Est",
      Region %in% c("Toscana",
                    "Umbria",
                    "Marche",
                    "Lazio") ~ "Centro",
      Region %in% c("Abruzzo",
                    "Molise",
                    "Campania",
                    "Puglia",
                    "Basilicata",
                    "Calabria") ~ "Sud",
      Region %in% c("Sicilia", "Sardegna") ~ "Isole",
      .default = NA),
      Nuts1 = factor(Nuts1, levels=c("Isole", "Sud", "Centro","Nord-Est", "Nord-Ovest")))
}

filter_ages <- function(dt){
  dt <- dt %>% 
    filter(Age %in% c(#"25-29", "30-34", 
      "35-39", "40-44", 
      "45-49", "50-54", 
      "55-59", "60-64")) %>% 
    mutate(Age = case_when(
      #Age %in% c("25-29", "30-34") ~ "25-34",
      Age %in% c("35-39", "40-44") ~ "35-44",
      Age %in% c("45-49", "50-54") ~ "45-54",
      Age %in% c("55-59", "60-64") ~ "55-64",
      .default = NA))
} 

from_code_to_province <- function(dt){
  dt <- dt %>% 
    mutate(Province = case_when( 
      Place == "001" ~ "Torino",
      Place == "002" ~ "Vercelli",
      Place == "003" ~ "Novara",
      Place == "004" ~ "Cuneo",
      Place == "005" ~ "Asti",
      Place == "006" ~ "Alessandria",
      Place == "007" ~ "Valle d'Aosta / Vallée d'Aoste",
      Place == "008" ~ "Imperia",
      Place == "009" ~ "Savona",
      Place == "010" ~ "Genova",
      Place == "011" ~ "La Spezia",
      Place == "012" ~ "Varese",
      Place == "013" ~ "Como",
      Place == "014" ~ "Sondrio",
      Place == "015" ~ "Milano",
      Place == "016" ~ "Bergamo",
      Place == "017" ~ "Brescia",
      Place == "018" ~ "Pavia",
      Place == "019" ~ "Cremona",
      Place == "020" ~ "Mantova",
      Place == "021" ~ "Bolzano",
      Place == "022" ~ "Trento",
      Place == "023" ~ "Verona",
      Place == "024" ~ "Vicenza",
      Place == "025" ~ "Belluno",
      Place == "026" ~ "Treviso",
      Place == "027" ~ "Venezia",
      Place == "028" ~ "Padova",
      Place == "029" ~ "Rovigo",
      Place == "030" ~ "Udine",
      Place == "031" ~ "Gorizia",
      Place == "032" ~ "Trieste",
      Place == "033" ~ "Piacenza",
      Place == "034" ~ "Parma",
      Place == "035" ~ "Reggio nell'Emilia",
      Place == "036" ~ "Modena",
      Place == "037" ~ "Bologna",
      Place == "038" ~ "Ferrara",
      Place == "039" ~ "Ravenna",
      Place == "040" ~ "Forlì-Cesena",
      Place == "041" ~ "Pesaro e Urbino",
      Place == "042" ~ "Ancona",
      Place == "043" ~ "Macerata",
      Place == "044" ~ "Ascoli Piceno",
      Place == "045" ~ "Massa-Carrara",
      Place == "046" ~ "Lucca",
      Place == "047" ~ "Pistoia",
      Place == "048" ~ "Firenze",
      Place == "049" ~ "Livorno",
      Place == "050" ~ "Pisa",
      Place == "051" ~ "Arezzo",
      Place == "052" ~ "Siena",
      Place == "053" ~ "Grosseto",
      Place == "054" ~ "Perugia",
      Place == "055" ~ "Terni",
      Place == "056" ~ "Viterbo",
      Place == "057" ~ "Rieti",
      Place == "058" ~ "Roma",
      Place == "059" ~ "Latina",
      Place == "060" ~ "Frosinone",
      Place == "061" ~ "Caserta",
      Place == "062" ~ "Benevento",
      Place == "063" ~ "Napoli",
      Place == "064" ~ "Avellino",
      Place == "065" ~ "Salerno",
      Place == "066" ~ "L'Aquila",
      Place == "067" ~ "Teramo",
      Place == "068" ~ "Pescara",
      Place == "069" ~ "Chieti",
      Place == "070" ~ "Campobasso",
      Place == "071" ~ "Foggia",
      Place == "072" ~ "Bari",
      Place == "073" ~ "Taranto",
      Place == "074" ~ "Brindisi",
      Place == "075" ~ "Lecce",
      Place == "076" ~ "Potenza",
      Place == "077" ~ "Matera",
      Place == "078" ~ "Cosenza",
      Place == "079" ~ "Catanzaro",
      Place == "080" ~ "Reggio di Calabria",
      Place == "081" ~ "Trapani",
      Place == "082" ~ "Palermo",
      Place == "083" ~ "Messina",
      Place == "084" ~ "Agrigento",
      Place == "085" ~ "Caltanissetta",
      Place == "086" ~ "Enna",
      Place == "087" ~ "Catania",
      Place == "088" ~ "Ragusa",
      Place == "089" ~ "Siracusa",
      Place == "090" ~ "Sassari",
      Place == "091" ~ "Nuoro",
      Place == "092" ~ "Cagliari",
      Place == "093" ~ "Pordenone",
      Place == "094" ~ "Isernia",
      Place == "095" ~ "Oristano",
      Place == "096" ~ "Biella",
      Place == "097" ~ "Lecco",
      Place == "098" ~ "Lodi",
      Place == "099" ~ "Rimini",
      Place == "100" ~ "Prato",
      Place == "101" ~ "Crotone",
      Place == "102" ~ "Vibo Valenzia",
      Place == "103" ~ "Verbano-Cusio-Ossola",
      #Place == "104" ~ "Olbia Tempio",
      #Place == "105" ~ "Ogliastra",
      #Place == "106" ~ "Medio Campidano",
      #Place == "107" ~ "Carbonia-Iglesias",
      Place == "104" ~ "Sassari",
      Place == "105" ~ "Nuoro",
      # Place == "106" ~ "Sud Sardegna",
      # Place == "107" ~ "Sud Sardegna",
      Place == "106" ~ "Cagliari",
      Place == "107" ~ "Cagliari",
      Place == "108" ~ "Monza e della Brianza",
      Place == "109" ~ "Fermo",
      Place == "110" ~ "Barletta-Andria-Trani",
      #Place == "111" ~ "Sud Sardegna",
      Place == "111" ~ "Cagliari",
      .default = NA)) %>% 
    relocate(Province, .before=Place)
}

correct_name_province <- function(dt){
  dt <- dt %>% 
    mutate(Province = case_when( 
      Province == "Valle d'Aosta / Vallée d'Aoste" ~ "Aosta",
      Province == "Massa-Carrara" ~ "Massa Carrara",
      Province == "Vibo Valenzia" ~ "Vibo Valentia",
      .default = Province))
}





