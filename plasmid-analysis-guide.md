# Plasmid analysis

## Downloading plasmids:

I downloaded plasmids from the following three sources:

PLSDB: version 2023_11_03 fasta and tab separated metadata file
Lerminiaux KPC paper: 202 KPC plasmids
Lerminiaux OXA paper: 70 OXA plasmids

PLSDB plasmids were downloaded directly from their website.

Lerminiaux plasmids were downloaded using BatchEntrez on NCBI by submitting text file lists of accessions and downloading results as a fasta file.

## Creating input files for creating mob-suite db:

### Create taxonomy files:

For oxa plasmids:
`bioawk -t -c fastx '{split($comment, a, " "); print $name" "$comment, a[1]"_"a[2]}' oxa-plasmids.fasta >> oxa-plasmids-taxonomy.txt`

For kpc plasmids:
`bioawk -t -c fastx '{split($comment, a, " "); print $name" "$comment, a[1]"_"a[2]}' kpc-plasmids.fasta >> kpc-plasmids-taxonomy.txt`

For plsdb plasmids:
`bioawk -t -c fastx '{split($comment, a, " "); print $name" "$comment, a[1]"_"a[2]}' plsdb.fna >> plsdb-plasmids-taxonomy.txt`

Combine all tax files:
`cat oxa-plasmids-taxonomy.txt kpc-plasmids-taxonomy.txt plsdb-plasmids-taxonomy.txt > all-plasmids-taxonomy.txt`

### Run mob-suite typer on all downloaded plasmids:

Create combined fasta file of plasmids:
`cat oxa-plasmids.fasta kpc-plasmids.fasta plsdb.fna > all-plasmids.fasta`

Run mobsuite typer on all plasmids:

It crashed when I tried to run all at once, so I split the fasta of all plasmids into files of 1000 plasmid seq each and ran mobtyper on each, then I combined them into one results file:

```
(mobsuite-env) user@user-ThinkStation-P3-Tower:~/cpe/plasmid-analysis/all-plasmid-seqs-1000$ bioawk -c fastx '{
    if (NR % 1000 == 1) {
        file = sprintf("chunk_%03d.fasta", int((NR-1)/1000) + 1)
    }
    print ">"$name" "$comment"\n"$seq > file
}' all-plasmids.fasta
```

`(mobsuite-env) user@user-ThinkStation-P3-Tower:~/cpe/plasmid-analysis/all-plasmid-seqs-1000$ for i in *.fasta; do mob_typer --multi --infile ${i} --out_file ${i}-mobtyper-results.txt --num_threads 24 -a ${i}-temp/; done`

Combine results into one file:
`(base) user@user-ThinkStation-P3-Tower:~/cpe/plasmid-analysis/all-plasmid-seqs-1000$ awk 'FNR==1 && NR!=1 {next} {print}' chunk*.txt > all-plasmids-mobtyper-results.txt`

### Set up run of mob cluster --mode build (build database):

After messing around with it, I had to rename plasmids only by their accessions in the fasta, taxonomy, and mobtyper results files. This was done like so:

Taxonomy:
`(mobsuite-env) user@user-ThinkStation-P3-Tower:~/cpe/plasmid-analysis$ awk 'NR==1 {print $0; next} {sub(" .*", "", $1); print $0}' FS='\t' OFS='\t' all-plasmids-taxonomy.txt > edited-all-plasmids-taxonomy.txt`

mobtyper results:
`(mobsuite-env) user@user-ThinkStation-P3-Tower:~/cpe/plasmid-analysis$ awk 'NR==1 {print $0; next} {sub(" .*", "", $1); print $0}' FS='\t' OFS='\t' all-plasmid-seqs-1000/all-plasmids-mobtyper-results.txt > edited-all-plasmids-mobtyper-results.txt`

fasta:
`(mobsuite-env) user@user-ThinkStation-P3-Tower:~/cpe/plasmid-analysis$ bioawk -c fastx '{split($name, a, " "); print ">"a[1]"\n"$seq}' all-plasmids.fasta > edited-all-plasmids.fasta`

Then I ran mobcluster:
`(mobsuite-env) user@user-ThinkStation-P3-Tower:~/cpe/plasmid-analysis$ mob_cluster --mode build -f edited-all-plasmids.fasta -t edited-all-plasmids-taxonomy.txt -p edited-all-plasmid`

Which seems to have built the database properly...

## Run mobsuite recon on my genomes

`(mobsuite-env) user@user-ThinkStation-P3-Tower:~/cpe/plasmid-analysis$ for i in $(cat ~/cpe/npa-list.txt); do mob_recon --infile ~/cpe/de-novo-assemblies/polished-assemblies/${i}-polypolish-assembly.fasta --outdir mobrecon-results/${i}/ --plasmid_meta all-plasmids-db/clusters.txt --plasmid_db all-plasmids-db/references_updated.fasta --plasmid_mash_db all-plasmids-db/references_updated.fasta.msh --num_threads 28 -u -s ${i}; done`

## Get contigs with carbapenemases for each isolate:

### Start with unicycler assemblies + card rgi output genomes:
I wrote a script to parse this:

```
(base) user@user-ThinkStation-P3-Tower:~/cpe/plasmid-analysis$ for i in $(cat ~/cpe/npa-list.txt); do python ~/postdoc-python-scripts/cpe11_get_carba_contigs.py --sample ${i} --card_input ~/cpe/genome-annotations/cpe9-card-rgi/${i}-CARD-RGI.txt --gff ~/cpe/genome-annotations/bakta-annotations/${i}/${i}.gff3 --outfile ~/cpe/plasmid-analysis/cpe11-contigs-with-carbapenemases.txt; done
Traceback (most recent call last):
  File "/home/user/postdoc-python-scripts/cpe11_get_carba_contigs.py", line 91, in <module>
    main()
  File "/home/user/postdoc-python-scripts/cpe11_get_carba_contigs.py", line 85, in main
    parse_rgi(args.card_input, args.gff, args.sample, args.outfile)
  File "/home/user/postdoc-python-scripts/cpe11_get_carba_contigs.py", line 59, in parse_rgi
    with open(card_handle, "r") as infile1:
         ^^^^^^^^^^^^^^^^^^^^^^
FileNotFoundError: [Errno 2] No such file or directory: '/home/user/cpe/genome-annotations/cpe9-card-rgi/KSTARO065-KS-CARD-RGI.txt'
```

KSTARO065 is known to have failed assembly.

I then filtered the output file to remove blank lines and oxa-1, -2, -9, and -10 (not carbapenemases):
`(mobsuite-env) user@user-ThinkStation-P3-Tower:~/cpe/plasmid-analysis$ awk -F'\t' 'NF && !($4 == "OXA-1" || $4 == "OXA-9" || $4 == "OXA-2" || $4 == "OXA-10")' cpe11-contigs-with-carbapenemases.txt > EDITED-cpe11-contigs-with-carbapenemases.txt`

### Now do the same for SKESA assemblies (but annotate with Bakta first!):

Run polypolish on 6 skesa assemblies:
```
# bwa index
(polypolish-polca-env) user@user-ThinkStation-P3-Tower:~/cpe/de-novo-assemblies/cpe9-skesa-missing-carbapenemases-reassemblies$ for i in $(cat ~/cpe/genome-annotations/skesa-isolates.txt); do bwa index ${i}-skesa-assembly.fasta; done

# bwa mem
(polypolish-polca-env) user@user-ThinkStation-P3-Tower:~/cpe/de-novo-assemblies/cpe11-polypolish-skesa$ for i in $(cat ~/cpe/genome-annotations/skesa-isolates.txt); do bwa mem -t 24 -a ../cpe9-skesa-missing-carbapenemases-reassemblies/${i}-skesa-assembly.fasta /media/user/DATA/cpe/trimmed-fastq-files/${i}_1.fastq.gz > ${i}-alignments_1.sam; bwa mem -t 24 -a ../cpe9-skesa-missing-carbapenemases-reassemblies/${i}-skesa-assembly.fasta /media/user/DATA/cpe/trimmed-fastq-files/${i}_2.fastq.gz > ${i}-alignments_2.sam; done

# run polypolish steps:
(polypolish-polca-env) user@user-ThinkStation-P3-Tower:~/cpe/de-novo-assemblies/cpe11-polypolish-skesa$ for i in $(cat ~/cpe/genome-annotations/skesa-isolates.txt); do polypolish filter --in1 ${i}-alignments_1.sam --in2 ${i}-alignments_2.sam --out1 ${i}-filtered_1.sam --out2 ${i}-filtered_2.sam; polypolish polish ../cpe9-skesa-missing-carbapenemases-reassemblies/${i}-skesa-assembly.fasta ${i}-filtered_1.sam ${i}-filtered_2.sam > ${i}-skesa-polypolish-assembly.fasta; done

# manually deleted intermediate files
```

Annotate with bakta:
`(bakta-1.9.3) user@user-ThinkStation-P3-Tower:~/cpe/genome-annotations/cpe11-skesa-bakta$ for i in $(cat ../skesa-isolates.txt); do bakta --db /home/user/bakta-dbs/db --min-contig-length 200 --prefix ${i} --output ${i}/ --strain ${i} --compliant --threads 24 --verbose ~/cpe/de-novo-assemblies/cpe11-polypolish-skesa/${i}-skesa-polypolish-assembly.fasta; done`

Re-run CARD on bakta annotations (ffn) (we'll see if bakta annotated the carbapenemases even...):
`(base) user@user-ThinkStation-P3-Tower:~/cpe/genome-annotations/cpe11-skesa-bakta$ for i in $(cat ../skesa-isolates.txt); do sudo docker run -v $PWD:/data quay.io/biocontainers/rgi:6.0.3--pyha8f3691_0 rgi main -i /data/${i}/${i}.ffn -o /data/${i}-CARD-RGI -n24 --clean; done`

Then parse:
`(base) user@user-ThinkStation-P3-Tower:~/cpe/plasmid-analysis$ for i in $(cat ~/cpe/genome-annotations/skesa-isolates.txt); do python ~/postdoc-python-scripts/cpe11_get_carba_contigs.py --sample ${i} --card_input ~/cpe/genome-annotations/cpe11-skesa-bakta/${i}-CARD-RGI.txt --gff ~/cpe/genome-annotations/cpe11-skesa-bakta/${i}/${i}.gff3 --outfile ~/cpe/plasmid-analysis/cpe11-skesa-with-carbapenemases.txt; done`

Manually deleted blank lines/wrong oxa enzymes...

Yes! It seems to have recovered the correct carbapenemases from the bakta annotations.

### Run mob-recon for skesa assemblies
(mobsuite-env) user@user-ThinkStation-P3-Tower:~/cpe/plasmid-analysis$ for i in $(cat ~/cpe/genome-annotations/skesa-isolates.txt); do mob_recon --infile ~/cpe/de-novo-assemblies/cpe11-polypolish-skesa/${i}-skesa-polypolish-assembly.fasta --outdir mobrecon-results-skesa/${i}/ --plasmid_meta all-plasmids-db/clusters.txt --plasmid_db all-plasmids-db/references_updated.fasta --plasmid_mash_db all-plasmids-db/references_updated.fasta.msh --num_threads 24 -s ${i}; done


## Parse mob-recon results based on contigs with carbapenemases -- NOTE that isolates that only have a carbapenemase by KMA are not included here --

I will parse this separately for unicycler/skesa-based mobrecon results.

First for unicycler:
`(base) user@user-ThinkStation-P3-Tower:~/postdoc-python-scripts$ for i in $(cat ~/cpe/plasmid-analysis/carbapenemase-isolate-list-unicycler.txt); do python cpe11_parse_mobrecon_results.py --sample ${i} --mobrecon_contigs ~/cpe/plasmid-analysis/mobrecon-results/${i}/contig_report.txt --mobtyper_results ~/cpe/plasmid-analysis/mobrecon-results/${i}/mobtyper_results.txt --input_file ~/cpe/plasmid-analysis/EDITED-cpe11-contigs-with-carbapenemases.txt --assembler unicycler --output_contigs ~/cpe/plasmid-analysis/cpe11-mobrecon-contigs-RESULTS.txt --output_plasmids ~/cpe/plasmid-analysis/cpe11-mobrecon-plasmids-RESULTS.txt; done`

Then skesa:
`(base) user@user-ThinkStation-P3-Tower:~/postdoc-python-scripts$ for i in $(cat ~/cpe/genome-annotations/skesa-isolates.txt); do python cpe11_parse_mobrecon_results.py --sample ${i} --mobrecon_contigs ~/cpe/plasmid-analysis/mobrecon-results-skesa/${i}/contig_report.txt --mobtyper_results ~/cpe/plasmid-analysis/mobrecon-results-skesa/${i}/mobtyper_results.txt --input_file ~/cpe/plasmid-analysis/cpe11-skesa-with-carbapenemases.txt --assembler skesa --output_contigs ~/cpe/plasmid-analysis/cpe11-mobrecon-contigs-RESULTS2.txt --output_plasmids ~/cpe/plasmid-analysis/cpe11-mobrecon-plasmids-RESULTS2.txt; done`

Note there are 2 sets of results files, one for unicycler and one for skesa. I will polish them up/look at them later.






## Files created during these analyses:

all-plasmids.fasta: fasta file of all plasmids used to create the mobsuite db; this includes ~55K PLSDB plasmids, 202 Lerminiaux KPC plasmids, and 70 Lerminiaux OXA plasmids

all-plasmids-taxonomy.txt: taxonomy file for mobrecon for all plasmids used to create mobsuite db; caveat is that this file has plasmids named by full names, which did not work when running mob_cluster; see edited-all-plasmids-taxonomy.txt

carbapenemase-isolate-list-unicycler.txt: list of isolates with carbapenemases detected from unicycler assemblies by card rgi; used to parse mobrecon results to avoid errors for isolates that do not have any carbapenemases but may have other plasmids

cpe11-contigs-with-carbapenemases.txt: output of running the cpe11_get_carba_contigs.py python script; unedited output that contains non-carbapenemase OXA genes and blank lines where isolates were missing carbapenemases; generated using the mentioned python script and CARD RGI output

EDITED-cpe11-contigs-with-carbapenemases.txt: edited version of this file that removed oxa-1, -2, -9, and -10, as well as blank lines; this is the good version used downstream (for running cpe11_parse_mobrecon_results.py)

cpe11-skesa-with-carbapenemases.txt: similar to cpe11-contigs-with-carbapenemases.txt but based on SKESA assemblies for isolates that had no carbapenemases recovered by unicycler but did by skesa; only includes isolates that DID have a carbapenemase recovered by skesa (N=6); was manually edited to remove blanks; this is the good version used downstream (for running cpe11_parse_mobrecon_results.py)

edited-all-plasmids.fasta: edited version of this file that renamed all contigs to just their accession; used to generate mobsuite db

edited-all-plasmids.fasta.msh: not actually sure about this file...do not remove! I think generated by running mobtyper prior to building db using mob_cluster --mode build...

edited-all-plasmids-mobtyper-results.txt: results of mobtyper run on all plasmids, where plasmids have been renamed only by their accessions

edited-all-plasmids-taxonomy.txt: taxonomy file for building mobsuite db, where plasmids have been renamed only by their accessions

kpc-plasmids.fasta: Lerminiaux fasta of KPC plasmids
kpc-plasmids.txt: Lerminiaux accession list (one per line) of KPC plasmids
kpc-plasmids.xlsx: Lerminiaux KPC plasmids metadata
kpc-plasmids-taxonomy.txt: Lerminiaux KPC plasmids taxonomy file; contatenated with others to make all-plasmids-taxonomy.txt

OXA and PLSDB similarly named files to set of 4 above are equivalents.




