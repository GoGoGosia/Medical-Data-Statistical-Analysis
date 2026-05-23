suppressPackageStartupMessages({
  library(Hmisc)      
  library(dplyr)      
  library(ggpubr)    
  library(car)        
  library(FSA)        
  library(dunn.test)  
  library(corrplot)  
})
options(width = 1000)
raport <- list()
argumenty <- commandArgs(trailingOnly = TRUE)

if(length(argumenty) == 0) {
  plik <- "przykladoweDane-Projekt.csv"
  cat("Nie podano pliku wejściowego\n")
  cat("Wczytano domyślny plik: ", plik, "\n\n")
} else {
  plik <- argumenty[1]
  cat("Wczytano plik:", plik, "\n\n")
}

dane <- read.csv2(plik, stringsAsFactors = FALSE, dec = ",")
raport$nazwaPliku <- plik

nazwaKolGrupa <- colnames(dane)[1]
dane[[nazwaKolGrupa]] <- as.factor(dane[[nazwaKolGrupa]])
poziomyGrup <- levels(dane[[nazwaKolGrupa]])

indeksyNumeryczne <- which(sapply(dane, is.numeric))
kolumnyNumeryczne <- names(dane)[indeksyNumeryczne]

tabelaZmian <- data.frame(Kolumna = character(), BrakiDanych = numeric(), WstawionaSrednia = numeric())
listaWartosciOdstajacych <- list()

for(i in indeksyNumeryczne) {
  nazwaKol <- colnames(dane)[i]
  
  # Imputacja
  liczbaNa <- sum(is.na(dane[, i]))
  if(liczbaNa > 0) {
    sredniaVal <- mean(dane[, i], na.rm = TRUE)
    tabelaZmian <- rbind(tabelaZmian, data.frame(
      Kolumna = nazwaKol, 
      BrakiDanych = liczbaNa, 
      WstawionaSrednia = round(sredniaVal, 2)
    ))
    dane[, i] <- as.numeric(impute(dane[, i], sredniaVal))
  }
  
  # Outliery
  outliers <- boxplot.stats(dane[[i]])$out
  if(length(outliers) > 0) {
    listaWartosciOdstajacych[[nazwaKol]] <- outliers
  }
}
raport$tabelaZmian <- tabelaZmian
raport$outliers <- listaWartosciOdstajacych

#Charakterystyka Grup
charSzeroka <- dane %>%
  group_by(.data[[nazwaKolGrupa]]) %>%
  summarise(across(all_of(kolumnyNumeryczne), list(
    Sr = ~round(mean(., na.rm = TRUE), 2),
    Md = ~round(median(., na.rm = TRUE), 2),
    SD = ~round(sd(., na.rm = TRUE), 2)
  ), .names = "{.col}_{.fn}"))

charLong <- as.data.frame(t(charSzeroka[,-1]))
colnames(charLong) <- charSzeroka[[1]]

charLong <- cbind(Parametr = rownames(charLong), charLong)
raport$charakterystyka <- charLong

#Analiza Porównawcza
listaWynikowTestow <- list()
for (param in kolumnyNumeryczne) {
  formulaTest <- as.formula(paste(param, "~", nazwaKolGrupa))
  
  shapiroP <- tapply(dane[[param]], dane[[nazwaKolGrupa]], function(x) shapiro.test(x)$p.value)
  normalny <- all(shapiroP > 0.05)
  wariancjaP <- leveneTest(formulaTest, data = dane)$`Pr(>F)`[1]
  jednorodna <- wariancjaP > 0.05
  
  opisPostHoc <- "Brak istotnych różnic"
  
  if (normalny && jednorodna) {
    testType <- "Parametryczny (ANOVA)"
    modelAov <- aov(formulaTest, data = dane)
    glownyP <- summary(modelAov)[[1]][["Pr(>F)"]][1]
    if (glownyP < 0.05) {
      ph <- TukeyHSD(modelAov)
      pary <- ph[[nazwaKolGrupa]][ph[[nazwaKolGrupa]][, "p adj"] < 0.05, , drop = FALSE]
      if (nrow(pary) > 0) opisPostHoc <- paste(rownames(pary), collapse = "; ")
    }
  } else {
    testType <- "Nieparametryczny (Kruskal-Wallis)"
    testKW <- kruskal.test(formulaTest, data = dane)
    glownyP <- testKW$p.value
    if (glownyP < 0.05) {
      ph <- dunnTest(formulaTest, data = dane, method = "bonferroni")
      pary <- ph$res[ph$res$P.adj < 0.05, ]
      if (nrow(pary) > 0) opisPostHoc <- paste(pary$Comparison, collapse = "; ")
    }
  }
  
  listaWynikowTestow[[param]] <- data.frame(
    Parametr = param, Normalnosc = ifelse(normalny, "OK", "Brak"), 
    Wariancja = ifelse(jednorodna, "OK", "Brak"), Test = testType, 
    Pvalue = round(glownyP, 5), Istotne_roznice = ifelse(glownyP < 0.05, "TAK", "NIE"),
    Miedzy_grupami = opisPostHoc
  )
}
raport$tabelaPorownawcza <- do.call(rbind, listaWynikowTestow)

#Analiza korelacji
listaKorelacji <- list()
for(grp in poziomyGrup) {
  podgrupa <- dane[dane[[nazwaKolGrupa]] == grp, ]
  pary <- combn(kolumnyNumeryczne, 2)
  for(i in 1:ncol(pary)) {
    p1 <- pary[1,i]; p2 <- pary[2,i]
    testKor <- cor.test(podgrupa[[p1]], podgrupa[[p2]], method = "spearman", exact = FALSE)
    if(testKor$p.value < 0.05) {
      listaKorelacji[[length(listaKorelacji) + 1]] <- data.frame(
        Grupa = grp, Parametr_1 = p1, Parametr_2 = p2, 
        Wspolczynnik_R = round(testKor$estimate, 3), Pvalue = round(testKor$p.value, 5)
      )
    }
  }
}
raport$tabelaKorelacji <- if(length(listaKorelacji) > 0) do.call(rbind, listaKorelacji) else "Brak istotnych statystycznie korelacji."

#RAPORT
sink("Raport.txt")
cat("RAPORT STATYSTYCZNEJ ANALIZY DANYCH\n")
cat("Plik wejściowy:", raport$nazwaPliku, "\n")

cat("\n1. PRZYGOTOWANIE DANYCH\n\n")
if(nrow(raport$tabelaZmian) > 0) {
  print(raport$tabelaZmian)
}else {
  cat("Brak braków danych.\n")
}

cat("\nWARTOŚCI ODSTAJĄCE:\n\n")
for(n in names(raport$outliers)) {
  cat(paste0(n, ": ", paste(raport$outliers[[n]], collapse=", "), "\n"))
}

cat("\n2. CHARAKTERYSTYKA GRUP\n\n")
print.data.frame(as.data.frame(raport$charakterystyka), row.names = FALSE)

cat("\n3. ANALIZA PORÓWNAWCZA\n\n")
print.data.frame(raport$tabelaPorownawcza, row.names = FALSE)

cat("\n4. ANALIZA KORELACJI W GRUPACH\n\n")
if(is.data.frame(raport$tabelaKorelacji)) {
  print.data.frame(raport$tabelaKorelacji, row.names = FALSE)
} else {
  cat(raport$tabelaKorelacji, "\n")
}
sink()

#WYKRESY
pdf("Wykresy.pdf", width = 10, height = 7)
liczbaGrup <- length(poziomyGrup)
kolory <- colorRampPalette(c("#ffc8dd", "#ffafcc", "#fb6f92"))(liczbaGrup)
koloryMacierzy <- colorRampPalette(c("#a2d2ff", "white", "#ffafcc"))(200)

# Boxplot
for(param in kolumnyNumeryczne) {
  Pval <- raport$tabelaPorownawcza$Pvalue[raport$tabelaPorownawcza$Parametr == param]
  boxplot(as.formula(paste(param, "~", nazwaKolGrupa)), 
          data = dane,
          main = paste("Rozkład:", param, "\np-value =", Pval),
          col = kolory,
          border = "#5D5D5D",
          xlab = "Grupa", ylab = "Wartość")
}

# Heatmapy
for(grp in poziomyGrup) {
  danePodgrupa <- dane[dane[[nazwaKolGrupa]] == grp, kolumnyNumeryczne]
  macierzKor <- cor(danePodgrupa, method = "spearman")
  
  corrplot(macierzKor, 
           method = "color", 
           type = "upper", 
           col = koloryMacierzy,   
           addCoef.col = "#4D4D4D", 
           tl.col = "black", tl.srt = 45, 
           title = paste("\nMacierz korelacji - Grupa:", grp),
           mar = c(0,0,3,0))
}

# Scatterplots
if(is.data.frame(raport$tabelaKorelacji)) {
  for(i in 1:nrow(raport$tabelaKorelacji)) {
    g <- raport$tabelaKorelacji$Grupa[i]
    p1 <- raport$tabelaKorelacji$Parametr_1[i]
    p2 <- raport$tabelaKorelacji$Parametr_2[i]
    
    kolorId <- which(poziomyGrup == g)
    dSub <- dane[dane[[nazwaKolGrupa]] == g, ]
    
    plot(dSub[[p1]], dSub[[p2]], 
         main = paste("Korelacja:", g, "|", p1, "vs", p2),
         xlab = p1, ylab = p2, 
         pch = 21, bg = kolory[kolorId], 
         col = "#5D5D5D", cex = 1.6)
    
    abline(lm(dSub[[p2]] ~ dSub[[p1]]), col = "#a2d2ff", lwd = 2.5) 
  }
}

dev.off()

cat("\nPROJEKT ZAKOŃCZONY!\nSprawdź pliki Raport.txt oraz Wykresy.pdf\n")