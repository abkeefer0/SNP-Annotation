#SNP annotation pipeline.
#Will need to update lines 8, 33, 36, 39, and 45 to match your setup.

#Stop script if an error occurs.
set -uxeo pipefail

#Prepare bcf_bgzip files downloaded from SNVPhyl.
for d in ~/Andrea/SNVPhyl/Copenhagen1/Files/* ;
do
	( cd $d && 
	#Convert bcf_bgzip files to vcf files.
	bcftools view *.bcf_bgzip -O v -o ${d}.vcf

	#Zip newly created vcf files.
	bgzip -c ${d}.vcf > ${d}.vcf.gz

	#Tabix newly created vcf files.
	tabix -fp vcf ${d}.vcf.gz

	#Pull out only variant positions.
	vcf-subset -t SNPs ${d}.vcf.gz > ${d}_snps.vcf

	#Zip new snps vcf file.
	bgzip -c ${d}_snps.vcf > ${d}_snps.vcf.gz

	#Tabix new snps vcf file.
	tabix -fp vcf ${d}_snps.vcf.gz

	);
done

#Combine all SNP only vcfs into one multisample vcf.
bcftools merge ~/Andrea/SNVPhyl/Copenhagen1/Files/*_snps.vcf.gz --force-samples -o ~/Andrea/SNVPhyl/Copenhagen1/combined_vcf.vcf

#Edit combined vcf so that reference name matches the SnpEff database name.
perl -pi -w -e 's/SC09_Chromosome/Chromosome/g;' ~/Andrea/SNVPhyl/Copenhagen1/combined_vcf.vcf

#Remove lines from combined vcf that start with '##'.
sed '/^##/ d' ~/Andrea/SNVPhyl/Copenhagen1/combined_vcf.vcf > ~/Andrea/SNVPhyl/Copenhagen1/final.vcf

#Set SnpEff database
Database=$1

#Run SnpEff to annotate variants.
java -Xmx4g -jar ~/Andrea/snpEFF/snpEff/snpEff.jar -v -strict ${Database} final.vcf > SnpEff.vcf

#Cut out important annotation information.
cat SnpEff.vcf | cut -f 8 | cut -d '|' -f 2,3,4,5,6,7,8,9,10,11,12,13,14 > annotation1.txt

#Cut out other pertinent information from original vcf file.
cat SnpEff.vcf | cut -f 1,2,3,4,5,6,7 > annotation2.txt

#Combine the two new annotation files.
paste annotation2.txt annotation1.txt | column -s $'\t' -t > annotated_all_snps.txt

#File with list of positions of "valid" SNPs derived from SNV table from SNVPhyl.
Valid=$2

#Grep for valid SNP positions in combined annotation file.
grep -f ${Valid} annotated_all_snps.txt > annotated_valid_snps.txt

#Count the number of annotation effects.
cat annotated_valid_snps.txt | tr -s ' ' | cut -d ' ' -f 8 | cut -d '|' -f 1 | sort | uniq -c > Effect_counts.txt

#Extract peg column from annotated_valid_snps.txt file.
cat annotated_valid_snps.txt | tr -s ' ' | cut -d ' ' -f 8 | cut -d '|' -f 4 > Peg_column.txt

#Check to see if some Pegs are duplicated, which can throw off the 3 below commands.
cat annotated_valid_snps.txt | tr -s ' ' | cut -d ' ' -f 8 | cut -d '|' -f 4 | sort | uniq -c | sort -nr > Multiple_pegs.txt

#Perform the 4 following commands to extract gene information. Best to run separately in terminal. Pegs_valid.txt file can be generated using the Peg_column.txt file in excel. It may require some editing due to GO values and not peg values being present. GFF file is from RAST.
#sed 's/$/;/' Pegs_valid.txt > Pegs_valid_2.txt
#grep -f Pegs_valid_2.txt SC09_PacBio_RAST.gff > Gene_names_1.txt 
#cat Gene_names_1.txt | cut -f 9 | cut -d ';' -f 2 > Gene_names_2.txt
#sed 's/Name=//g' Gene_names_2.txt > Gene_names_3.txt

#End of script.