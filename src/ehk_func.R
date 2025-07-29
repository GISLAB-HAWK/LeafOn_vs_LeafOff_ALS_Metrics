#--------------------------------------------------------------------------
# Name:         ehk_func.R
# Author:       Northwest German Forest Research Institute
#               Department Forest Growth
#--------------------------------------------------------------------------

# Funktionen für Einzelbaum-Höhenimputation über Einheitshöhenkurve

# 1. Erzeugen der Eingangstabelle im notwendingen Format

# 2. Eingangsgrößen für Höhenimputation
# 	dm: dg über alle Bäume
# 	ds: dg der Bäume mit Höhenmessung
# 	hs: hg der Bäume mit Höhenmessung
# 	hm: hg aller Bäume (über invertierte Funktion, vgl. Formel 5.2.4.9 BWI-Methodenband 2012)
# 	Berechnungsfolge für versch Ebenen mit Höhenmessung Plot/Schicht/LN/Bagr 
#   --> wenn für Ebene Plot/Schicht/LN/Bagr keine Höhenmessung verfügbar, wird Plot/Schicht/LN genommen usw
# 	globales Modell wenn am Plot keine Höhenmessung vorliegt
#_________________________________________________________________________________
# erstellt Eingangsdaten im benötigten Format
# id = eindeutige Plot-Id, bei BWI-Daten aus Tnr,Enr erstellen
# bagr = Baumartengruppe mit Koeffizienten für EHK (Stand BWI3)
# bhd in cm, hoe in m

input_ehk=function(id, bnr, bs, bhd, hoe, nha, bagr){
	
	ein = data.frame(id=id, bnr=bnr, bs=bs, bhd=bhd, hoe=hoe, nha=nha, bagr=bagr)

# Baumartengruppen und Koefizienten EHK (BWI 2012 Methodenband)
# bagr_ehk: ALN und ALH --> BU bei Höhenmessung, dh. ggf keine separaten Messungen
# ba_ln: Laubbaum/Nadelbaum. Wichtig bei Unterstand oder wenn nicht für alle bagr_ehk Höhenmessungen vorliegen 	
	z = data.frame(bagr=c("FI", "TA", "DGL", "KI", "LAE", "BU", "EI", "ALH", "ALN"),
							   ba_ln=c("N", "N", "N", "N", "N", "L", "L", "L", "L"),
								 bagr_ehk=c("FI", "TA", "DGL", "KI", "LAE", "BU", "EI", "BU", "BU"),
								 k0=c(0.183,0.079,0.24,0.29,0.074,0.032,0.102,0.122,0.032),
								 k1=c(5.688,3.992,6.033,1.607,3.692,6.04,3.387,5.04,4.42))

	ein=merge(ein,z, by="bagr")

# Ebenen festlegen 
# wird später benötigt um zu prüfen für welche Ebenen Messhöhen vorliegen.
# bei fehlenden Messhöhen innerhalb einer Ebene wird bei der Ergänzung die nächst höhere Ebene gewählt
# 1.Ebene: Baumartengruppe/LN/Schicht/Plot. 
# 2.Ebene: LN/Schicht/Plot. 
# 3.Ebene: Schicht/Plot. 
# 4.Ebene: Plot
	ein$id4 <-paste(ein$id, sep="")        					 # Plot
	ein$id3 <-paste(ein$id4, ein$bs, sep="_")        # Plot und Schicht
	ein$id2 <-paste(ein$id3, ein$ba_ln, sep="_")     # Plot, Schicht und Baumartengruppe (LN)                                                       
	ein$id1 <-paste(ein$id2, ein$bagr_ehk, sep="_")  # Plot, Schicht, Baumartengruppe (LN) und Baumart

	ein=ein[order(ein$id, ein$bnr), 
		c("id", "id4", "id3", "id2", "id1", "bnr", "bs", "bhd", "hoe", "nha", "bagr", "ba_ln", "bagr_ehk", "k0", "k1")]
	return(ein)
}


#_________________________________________________________________________________
# Berechnet dm, hm, ds, hs je Ebene (s.o) mit Gewichtung der Bäume (nha)
# dann Einzelbaumhöhen über EHK 
# für Plots ohne Messhöhe globales Modell h=f(d,bagr,bs) anwenden

ehk=function(ein){

  y=ein
	y$g=((pi/4)*y$bhd^2)*y$nha
  
	x1=aggregate(g ~ id1, data=y, FUN=sum)
	x1$n=aggregate(nha ~ id1, data=y, FUN=sum)$nha
	x1$dm1=sqrt((4*x1$g/x1$n)/pi)
	x1$g=NULL
	x1$n=NULL
	y=merge(y, x1, by="id1", all.x=T)
	
	x1=aggregate(g ~ id2, data=y, FUN=sum)
	x1$n=aggregate(nha ~ id2, data=y, FUN=sum)$nha
	x1$dm2=sqrt((4*x1$g/x1$n)/pi)
	x1$g=NULL
	x1$n=NULL
	y=merge(y, x1, by="id2", all.x=T)
	
	x1=aggregate(g ~ id3, data=y, FUN=sum)
	x1$n=aggregate(nha ~ id3, data=y, FUN=sum)$nha
	x1$dm3=sqrt((4*x1$g/x1$n)/pi)
	x1$g=NULL
	x1$n=NULL
	y=merge(y, x1, by="id3", all.x=T)

	x1=aggregate(g ~ id4, data=y, FUN=sum)
	x1$n=aggregate(nha ~ id4, data=y, FUN=sum)$nha
	x1$dm4=sqrt((4*x1$g/x1$n)/pi)
	x1$g=NULL
	x1$n=NULL
	y=merge(y, x1, by="id4", all.x=T)

# für Ausgabe merken	
	ein=y

# dg und hg Höhenmessbäume
# nur Höhenmessbäume behalten, am Ende aber mit Gesamtdaten verschneide
  y=y[!is.na(y$hoe) & y$hoe>0,] 
	y$gh=y$g*y$hoe
  
# dg
	x1=aggregate(g ~ id1, data=y, FUN=sum)
	x1$n=aggregate(nha ~ id1, data=y, FUN=sum)$nha
	x1$ds1=sqrt((4*x1$g/x1$n)/pi)
	x1$g=NULL
	x1$n=NULL
	ein=merge(ein, x1, by="id1", all.x=T)
	
	x1=aggregate(g ~ id2, data=y, FUN=sum)
	x1$n=aggregate(nha ~ id2, data=y, FUN=sum)$nha
	x1$ds2=sqrt((4*x1$g/x1$n)/pi)
	x1$g=NULL
	x1$n=NULL
	ein=merge(ein, x1, by="id2", all.x=T)
	
	x1=aggregate(g ~ id3, data=y, FUN=sum)
	x1$n=aggregate(nha ~ id3, data=y, FUN=sum)$nha
	x1$ds3=sqrt((4*x1$g/x1$n)/pi)
	x1$g=NULL
	x1$n=NULL
	ein=merge(ein, x1, by="id3", all.x=T)

	x1=aggregate(g ~ id4, data=y, FUN=sum)
	x1$n=aggregate(nha ~ id4, data=y, FUN=sum)$nha
	x1$ds4=sqrt((4*x1$g/x1$n)/pi)
	x1$g=NULL
	x1$n=NULL
	ein=merge(ein, x1, by="id4", all.x=T)

# hg
  x1=aggregate(gh ~ id1, data=y, FUN=sum)
  x1$gsum=aggregate(g ~ id1, data=y, FUN=sum)$g
  x1$hs1=x1$gh/x1$gsum
	x1$gh=NULL
	x1$gsum=NULL
	ein=merge(ein, x1, by="id1", all.x=T)
	
  x1=aggregate(gh ~ id2, data=y, FUN=sum)
  x1$gsum=aggregate(g ~ id2, data=y, FUN=sum)$g
  x1$hs2=x1$gh/x1$gsum
	x1$gh=NULL
	x1$gsum=NULL
	ein=merge(ein, x1, by="id2", all.x=T)
	
  x1=aggregate(gh ~ id3, data=y, FUN=sum)
  x1$gsum=aggregate(g ~ id3, data=y, FUN=sum)$g
  x1$hs3=x1$gh/x1$gsum
	x1$gh=NULL
	x1$gsum=NULL
	ein=merge(ein, x1, by="id3", all.x=T)
	
  x1=aggregate(gh ~ id4, data=y, FUN=sum)
  x1$gsum=aggregate(g ~ id4, data=y, FUN=sum)$g
  x1$hs4=x1$gh/x1$gsum
	x1$gh=NULL
	x1$gsum=NULL
	ein=merge(ein, x1, by="id4", all.x=T)

# Ebene mit Höhenmessung Messung behalten
	ein$dm=ein$dm4
	ein$dm=ifelse(!is.na(ein$ds3),ein$dm3,ein$dm)
	ein$dm=ifelse(!is.na(ein$ds2),ein$dm2,ein$dm)
	ein$dm=ifelse(!is.na(ein$ds1),ein$dm1,ein$dm)	
	
	ein$ds=ein$ds4
	ein$ds=ifelse(!is.na(ein$ds3),ein$ds3,ein$ds)
	ein$ds=ifelse(!is.na(ein$ds2),ein$ds2,ein$ds)
	ein$ds=ifelse(!is.na(ein$ds1),ein$ds1,ein$ds)

	ein$hs=ein$hs4
	ein$hs=ifelse(!is.na(ein$hs3),ein$hs3,ein$hs)
	ein$hs=ifelse(!is.na(ein$hs2),ein$hs2,ein$hs)
	ein$hs=ifelse(!is.na(ein$hs1),ein$hs1,ein$hs)
	
# Spalten löschen
	del=c("id1", "id2", "id3", "id4", "g", "dm1", "dm2", "dm3", "dm4", "ds1", "ds2", "ds3", "ds4", "hs1", "hs2", "hs3", "hs4")
		for(i in del)
		{
			ein[,i]=NULL
		}

# hm berechnen über invertierte EHK berechnen
	ein$hm=((ein$hs-1.3)/(exp(ein$k0*(1-(ein$dm/ein$ds)))*exp(ein$k0*((1/ein$dm)-(1/ein$ds)))))+1.3

# Höhenimputation für Einzelbäume über EHK
	ein$hoe_mod = 1.3 + (ein$hm-1.3) * exp(ein$k0*(1-(ein$dm/ein$bhd))) * exp(ein$k1*((1/ein$dm)-(1/ein$bhd)))


# globales Modell für Höhenimputation, wenn keine Messungen am Plot vorliegen
	l1=lm(hoe ~ log(bhd) + factor(bagr) + factor(bs), data=y) 
  ein$hoe_mod2=predict(l1, newdata=ein) 
  
  if(any(is.na(ein$hoe_mod))){
  	ein[is.na(ein$hoe_mod),]$hoe_mod=ein[is.na(ein$hoe_mod),]$hoe_mod2
  }
  ein$hoe_mod2=NULL
	
  return(ein)
}



#_________________________________________________________________________________
# globales Modell für Höhenimputation




