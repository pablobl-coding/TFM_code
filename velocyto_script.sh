#!/bin/bash
BAM=$1
barcodesfile=$2
gtf=$3

SUB='/'
if [[ "$BAM" == *"$SUB"* ]]; then # Sees if there are "/", to see if it is a route to the .BAM file, or just the file
	echo "option 1"
	BAMF=$(echo $BAM | awk -F "/" '{print $(NF)}') # If it's a route, it uses "/" as a delimiter and keeps just the last element, which is the file
	samtools view -h $BAM | grep -v '|stAb' | grep -v '|pAb' | samtools view -bS -o noAb_$BAMF # To eliminate all possible alignments with antibodies
	samtools index -b noAb_$BAMF # In order to perform the following steps, is necessary to index the file
	cat <(samtools view -HS noAb_$BAMF) <(samtools view noAb_$BAMF | grep "MA:Z:*"  | sed  "s/MA:Z:/UB:Z:/" ) | samtools view -b > final_$BAMF # It changes the "MA" tag to an "UB" tag and keeps only the lines that contain UMIs
	samtools index -b final_$BAMF # We, once again, index the file
	if [[ "$barcodesfile" == *"$SUB"* ]]; then # Checks if we have a route to the barcode file, or just the file
		echo "option 1.1"
		BF=$(echo $barcodesfile | awk -F "/" '{print $(NF)}') # Same thing as in the previous case
		cut -d ',' -f 1 $barcodesfile | grep -v '#' |  sed '1d' > new_$BF # To modify the barcodes file to an appropriate format
		samtools sort -t CB -O BAM -o cellsorted_final_$BAMF final_$BAMF # Ordering based on mapping
		velocyto run -b new_$BF final_$BAMF $gtf # We perform the velocyto run
	else
		echo "option 1.2" # Since we just have the file, it doesnt modify the input (the rest is exactly the same)
		cut -d ',' -f 1 $barcodesfile | grep -v '#' |  sed '1d' > new_$barcodesfile 
		samtools sort -t CB -O BAM -o cellsorted_final_$BAMF final_$BAMF
		velocyto run -b new_$barcodesfile final_$BAMF $gtf
	fi
# To remove the noAb files, which are unnecessary: rm noAb_$BAMF
else
	echo "option2" # Since we just have the file, it doesnt modify the input (the rest is exactly the same as in the other case)
	samtools view -h $BAM | grep -v '|stAb' | grep -v '|pAb' | samtools view -bS -o noAb_$BAM
	samtools index -b noAb_$BAM
	cat <(samtools view -HS noAb_$BAM) <(samtools view noAb_$BAM | grep "MA:Z:*"  | sed  "s/MA:Z:/UB:Z:/" ) | samtools view -b > final_$BAM
	samtools index -b final_$BAM
	if [[ "$barcodesfile" == *"$SUB"* ]]; then
		echo "option 2.1"
		BF=$(echo $barcodesfile | awk -F "/" '{print $(NF)}')
		cut -d ',' -f 1 $barcodesfile | grep -v '#' |  sed '1d' > new_$BF
		samtools sort -t CB -O BAM -o cellsorted_final_$BAM final_$BAM
		velocyto run -b new_$BF final_$BAM $gtf
	else
		echo "option 2.2"
		cut -d ',' -f 1 $barcodesfile | grep -v '#' |  sed '1d' > new_$barcodesfile
		samtools sort -t CB -O BAM -o cellsorted_final_$BAM final_$BAM
		velocyto run -b new_$barcodesfile final_$BAM $gtf
	fi
# To remove the noAb files, which are unnecessary: rm noAb_$BAM
fi
