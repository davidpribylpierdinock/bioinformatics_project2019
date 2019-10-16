#creates list of proteomes that have the McrA gene and the number of copies of the HSP70 gene in that proteome
#usage: bash methanogens.sh 

#STEP 1: compile all reference sequences into one fasta file

cat ref_sequences/mcrA* >> mcrA.ref
cat ref_sequences/hsp70* >> hsp70.ref

#STEP 2: run muscle and hmmbuild on each compiled reference fasta

for gene in *.ref
do
	executables/muscle -in $gene -out $gene_aligned.fasta
	executables/hmmbuild $gene.hmm $gene_aligned.fasta
done

#STEP 3: use the hmm to find methane producing proteomes

for proteome in proteome_*.fasta
do 
        executables/hmmsearch mcrA.ref.hmm $proteome >> mcrAsearchOut.txt
done

#STEP 4: 

for proteome in proteome_*.fasta
do
        executables/hmmsearch hsp70.ref.hmm $proteome >> hsp70searchOut.txt
done
