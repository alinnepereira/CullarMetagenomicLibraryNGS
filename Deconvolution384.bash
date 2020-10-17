#!/bin/bash

###BUILDING THE DECONVOLUTION INDEX

echo "Time Start BUILDING THE DECONVOLUTION INDEX"
date

##Finding Plate IDs sucessfully deconvoluted

echo "Time Start Finding Plate IDs successfully deconvoluted"
date

echo "Finding Plate IDs successfully deconvoluted"

awk '{print $1}' PlateAllxCOLUMN.RESULTS.txt | sort | uniq > PlateAllxCOLUMN.deconvolutedIDs.temp1.txt
awk '{print $1}' PlateAllxROW.RESULTS.txt | sort | uniq > PlateAllxROW.deconvolutedIDs.temp1.txt

cat *deconvolutedIDs.temp1.txt > PlateAllconcat.deconvolutedIDs.temp1.txt

sort PlateAllconcat.deconvolutedIDs.temp1.txt | uniq -d > PlateAll.deconvolutedIDs.temp.txt

echo "Time End Finding Plate IDs successfully deconvoluted"
date

###BULDING INDEX

echo "Time Start BULDING INDEX"
date

echo "Extracting best hits from Row and Column to build index."

##Adding BLAST headers to the index
# echo "Query Subject %id AlignmentLength Mismatches GapOpenings QueryStart QueryEnd SubjectStart SubjectEnd Evalue Bitscore" > Index.txt


awk 'FNR==NR { array[$0]++; next } $1 in array { print; delete array[$1] }' PlateAll.deconvolutedIDs.temp.txt PlateAllxCOLUMN.RESULTS.txt > IndexColumn.temp.txt

awk 'FNR==NR { array[$0]++; next } $1 in array { print; delete array[$1] }' PlateAll.deconvolutedIDs.temp.txt PlateAllxROW.RESULTS.txt > IndexRow.temp.txt

# ##Fixing subject ID to match names

perl -pi -e 's/:/|/g' IndexColumn.temp.txt
perl -pi -e 's/:/|/g' IndexRow.temp.txt

## Sorting Index

sort -k1 IndexColumn.temp.txt > IndexColumn.sorted.temp.txt
sort -k1 IndexRow.temp.txt > IndexRow.sorted.temp.txt

echo "Done!"

# Build translation table
	 
# echo "Building Name Index"
	 
# grep ">" PlateAll.fasta > ContigsNames.txt
# grep ">" ColumnAll.fasta >> ContigsNames.txt
# grep ">" RowAll.fasta >> ContigsNames.txt
	 
grep "|" ContigsNames.txt | awk '{print $1, $2}' | sed 's/>gnl|//g' > IDtranslationTable.txt
echo "Done!"
	
echo "Time End BULDING INDEX"
date

###BUILDING DECONVOLUITON INDEX

echo "Time Start BUILDING DECONVOLUITON INDEX"
date

echo "Building Deconvolution index"

##Printing the columns side by side

paste IndexColumn.temp.txt IndexRow.temp.txt | cut -f 1,2,14 > Deconvolution.temp.txt

##Translating file

head Deconvolution.temp.txt
P01A1.236	S1C04.3 	S1RH.2
P01A1.691	S1C04.42 	S1RC.20
P01B1.6		S1C05.4 	S1RD.1
P01B1.3440 	S1C04.42 	S1RC.20

###translating Column IDs
awk 'FNR==NR { a[$1]=$2; next } $2 in a { $2=a[$2] }1' IDtranslationTable.txt Deconvolution.temp.txt > ColumnTranslation.temp.txt

###Translating Row IDs

awk 'FNR==NR { a[$1]=$2; next } $3 in a { $3=a[$3] }1' IDtranslationTable.txt ColumnTranslation.temp.txt > Deconvolution.temp.txt

##Generating triplets and coordinates:

# Generating 96 coordintes
awk '{print substr($1, 1, 3), substr($1, 4, 2) substr($2, 3, 3) substr($3, 3, 2)}' Deconvolution.temp.txt > 96positions.temp.txt

# Generating 384 coordintes
awk 'NR == FNR { a[$1] = $2; next } $1 in a { $1 = a[$1] } $2 in a { $2 = a[$2] } 1' translationtable384.txt 96positions.temp.txt > 384positions.temp.txt

#Generating final Index
paste 96positions.temp.txt 384positions.temp.txt Deconvolution.temp.txt | awk '{print $1 substr($2, 1, 2), substr($2, 3, 3), substr($2, 6, 2), $3, substr($4, 1, 3), substr($4, 4, 2), $5, $6, $7}' | sed -e 's/ /-/1' -e 's/ /-/1' | sed -e 's/ /-/2' -e 's/ /-/2' | sed '1 i\96wellCloneCoordinates 384wellCloneCoordinates PlateContig ColumnContig RowContig Function Flag' > DeconvolutionIndex.txt

######BUILDING PARTIAL DECONVOLUTION INDEX

##Extracting partially deconvoluted IDs

echo "Building Partial Deconvolution index"

#Generating list of contigs
awk '{print $1}' PlateAllxCOLUMN.RESULTS.txt | sort | uniq > PlateAllxCOLUMN.deconvolutedIDs.temp1.txt
awk '{print $1}' PlateAllxROW.RESULTS.txt | sort | uniq > PlateAllxROW.deconvolutedIDs.temp1.txt

cat *deconvolutedIDs.temp1.txt > PlateAllconcat.deconvIDs.temp1.txt

sort PlateAllconcat.deconvIDs.temp1.txt | uniq -u > PlateAll.partiallydeconvolutedIDs.temp.txt

#Extracting
awk 'FNR==NR { array[$0]++; next } $1 in array { print; delete array[$1] }' PlateAll.partiallydeconvolutedIDs.temp.txt PlateAllxCOLUMN.RESULTS.txt > IndexColumnPartial.temp.txt

awk 'FNR==NR { array[$0]++; next } $1 in array { print; delete array[$1] }' PlateAll.partiallydeconvolutedIDs.temp.txt PlateAllxROW.RESULTS.txt > IndexRowPartial.temp.txt

##Fixing subject ID to match names

perl -pi -e 's/:/|/g' IndexColumnPartial.temp.txt
perl -pi -e 's/:/|/g' IndexRowPartial.temp.txt

## Sorting Index

sort -k1 IndexColumnPartial.temp.txt > IndexColumnPartial.sorted.temp.txt
sort -k1 IndexRowPartial.temp.txt > IndexRowPartial.sorted.temp.txt

##Translating file

###translating Column IDs
awk 'FNR==NR { a[$1]=$2; next } $2 in a { $2=a[$2] }1' IDtranslationTable.txt IndexColumnPartial.sorted.temp.txt > ColumnTranslation.Partial.temp.txt

###Translating Row IDs

awk 'FNR==NR { a[$1]=$2; next } $2 in a { $2=a[$2] }1' IDtranslationTable.txt IndexRowPartial.sorted.temp.txt > RowTranslation.Partial.temp.txt

##Printing the Column Partial Deconvolution

awk '{print $1, $2, "NA"}' ColumnTranslation.Partial.temp.txt | awk '{print substr($1, 1, 5), substr($2, 1, 5), $3, $0}' | sed -e 's/ /-/1' -e 's/ /-/1' > Deconvolution.partial.column.temp.txt

##Printing the Row Partial Deconvolution

awk '{print $1, "NA", $2}' RowTranslation.Partial.temp.txt | awk '{print substr($1, 1, 5), $2, substr($3, 1, 4), $0}' | sed -e 's/ /-/1' -e 's/ /-/1' > Deconvolution.partial.row.temp.txt

##Concatenate Partial Deconvolution

cat Deconvolution.partial.column.temp.txt Deconvolution.partial.row.temp.txt | sed '1 i\CloneCoordinates PlateContig ColumnContig RowContig Flags Function' > DeconvolutionIndex.partial.txt

#######EXTRACTING NOT-DECONVOLUTED CONTIGS LIST

echo "Building list of not-Deconvoluted contigs"

#Create an output file with three columns one for Plate contigs, Column Contigs, Row Contigs that was not deconvoluted! In order to do that I will use the Deconvolution Index:

##generating a list of all deconvoluted contigs(fully and partially)

awk 'FNR > 1 {print $2, $3, $4 }' DeconvolutionIndex.txt | awk 1 RS=" |\n" > DeconvolutedContigs.temp.txt

awk 'FNR > 1 {print $2, $3, $4 }' DeconvolutionIndex.partial.txt | awk 1 RS=" |\n" | sed '/NA/d' >> DeconvolutedContigs.temp.txt

sort DeconvolutedContigs.temp.txt | uniq > DeconvolutedContigs.sorted.temp.txt

awk '{print $2}' ContigsNames.txt > AllIds.temp.txt

cat DeconvolutedContigs.sorted.temp.txt AllIds.temp.txt | sort | uniq -u > NotDeconvolutedContigs.temp.txt

awk 'FNR==NR { array[$0]++; next } $2 in array { print; delete array[$2] }' NotDeconvolutedContigs.temp.txt ContigsNames.txt > DeconvolutionIndex.singletons.txt

echo "Time End BULDING INDEX"
date

echo "Removing Intermidiate files"

rm *temp*

