#! /bin/bash

gunzip sREGION.vcf.gz
vcfwave sREGION.vcf | gzip >  sREGION-vcfw.vcf.gz
gzip -cd sREGION-vcfw.vcf.gz | grep "^#" | grep -v "^#CHROM" > header1.tmp
gzip -cd sREGION-vcfw.vcf.gz | grep "^#CHROM" | awk '{print $0"\tsREF"}' > header2.tmp
gzip -cd sREGION-vcfw.vcf.gz | grep -v "^#" | awk -v len=sMIN_L '{if (length($4) >= len || length($5) >= len) print $0"\t"0}' > vars.tmp
cat header1.tmp header2.tmp vars.tmp | bgzip > sREGION-vcfw-sMIN_KBkb.vcf.gz
rm header*.tmp
rm vars.tmp

SET="sREGION-vcfw-sMIN_KBkb"

bcftools norm \
	-m- \
	-o sSET.norm.vcf \
	sSET.vcf.gz

grep "^#" sSET.norm.vcf > header.tmp
grep -v "^#" sSET.norm.vcf | awk '{if (length($4) > length($5)) print $0}' > dels.tmp
grep -v "^#" sSET.norm.vcf | awk '{if (length($4) < length($5)) print $0}' > ins.tmp

awk -v len=sMIN_L 'BEGIN { OFS = "\t" } { if (length($5) < len) $5=substr($4,1,1) } 1' dels.tmp > dels1.tmp
awk -v len=sMIN_L 'BEGIN { OFS = "\t" } { if (length($5) >= len) print $0 }' dels.tmp > dels2.tmp

awk -v len=sMIN_L 'BEGIN { OFS = "\t" } { if (length($4) < len) $4=substr($4,1,1) } 1' ins.tmp > ins1.tmp
awk -v len=sMIN_L 'BEGIN { OFS = "\t" } { if (length($4) >= len) print $0 }' ins.tmp > ins2.tmp

cat header.tmp dels1.tmp dels2.tmp ins1.tmp ins2.tmp | bcftools sort > sSET.norm.simplified.vcf
rm *.tmp

bgzip sSET.norm.vcf
bgzip sSET.norm.simplified.vcf

rsync -avuP sSET.norm*.gz aharder@mando.hagsc.org:sOUTDIR
