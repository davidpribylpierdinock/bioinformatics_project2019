#creates list of proteomes that have the McrA gene and the number of copies of the HSP70 gene in that proteome
#usage: bash methanogens.sh 

#STEP 1: clear working directory of conflicting files and compile all reference sequences into one fasta file

rm hsp70* mcrA*
cat ref_sequences/mcrA* >> mcrA.ref
cat ref_sequences/hsp70* >> hsp70.ref

#STEP 2: run muscle and hmmbuild on each compiled reference fasta

for gene in *.ref
do
	executables/muscle -in $gene -out $gene.fasta
	executables/hmmbuild $gene.hmm $gene.fasta
done

#STEP 3: use the hmm to find methane producing proteomes

for proteome in proteomes/*.fasta
do 
        executables/hmmsearch mcrA.ref.hmm $proteome >> mcrA_searchOut.txt
done

grep "Domain search space" mcrA_searchOut.txt | egrep -n [1-9][0-9]* | cut -d : -f 1 > mcrA_positive.txt

#STEP 4: use the hmm to search for pH resistance and refine to only methane producing proteomes

for proteome in proteomes/*.fasta
do
	executables/hmmsearch hsp70.ref.hmm $proteome >> hsp70_searchOut.txt
done

grep "Domain search space" hsp70_searchOut.txt | grep -n "Domain search space" > hsp70_counts.txt

for i in $(cat mcrA_positive.txt)
do 
	egrep "^$i:" hsp70_counts.txt | cut -d : -f 1,3 | sed -e 's/ \+/ /g' | sed 's/\[number of targets reported over threshold\]//g' >> output.txt
done

#STEP 5: add header to output file

echo "Methanogen proteome number: number of pH resistance genes" > finalOutput.txt
cat output.txt >> finalOutput.txt
cat finalOutput.txt
