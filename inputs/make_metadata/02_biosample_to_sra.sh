conda activate sandbox
conda install entrez-direct

cat refseq_prokaryotes_biosamples.txt | \
join-into-groups-of 500 | \
xargs -n 1 sh -c 'epost -db biosample -id "$0" -format acc | \
elink -target sra |  \
efetch -db sra -format runinfo -mode xml | \
xtract -pattern Row -def "NA" -element Run spots bases spots_with_mates avgLength \
size_MB download_path Experiment LibraryStrategy LibrarySelection LibrarySource \
LibraryLayout InsertSize InsertDev Platform Model SRAStudy BioProject ProjectID \
Sample BioSample SampleType TaxID ScientificName SampleName CenterName \
Submission Consent >> refseq_prokaryotes_metadata_sra.txt'
