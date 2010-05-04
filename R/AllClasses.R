

###SeqLogo###

setClass("pwm",representation(
pwm="matrix",
consensus="character",
ic="numeric",
width="numeric",
alphabet="character"))

###motiv###

setClass("transcriptionFactor", representation(
name="character",
pwm="matrix"
))

setClass("alignments", representation(
TF="transcriptionFactor",
evalue="numeric",
sequence="character",
match="character",
strand="character"
))	

setClass("matches", representation(
name="character",
aligns="list",
similarity="character",
valid="numeric"
))

setClass("motiv", representation(
input="list",
bestMatch="list",
argv="character"
))

###filter###

setClass("filter",
representation(
name="list",
tfname="list",
top="list",
evalueMax="list",
lengthMax="list",
valid="list"
),
prototype=list(
name=list(""),
tfname=list(""),
top=list(10),
evalueMax=list(1),	
lengthMax=list(100),
valid=list(1)
))

setClass("filters",
representation(
filters="list"
))

setClass("position", representation(
motifName="character",
positionVector="RangedData",
pwm="list",
similarity="character"
))

