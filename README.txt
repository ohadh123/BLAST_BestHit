The following documentation describes the instructions and logic behind the BLAST best-hit algorithm. The purpose of this program is to allow the user an alternative to sorting through hundreds of sequences manually; rather, the user can achieve the same (and better) results by allowing a fast, accurate, and consistent computer to perform such tasks. 

Note: For further assistance, the demo folder attached displays the command line and the files created as a result of running the program with blastp

#Prerequisites

Installation of BLAST+: The BLAST best-hit algorithm runs on a tabular file with comment lines, containing a series of matches given as output from performing BLAST through the command-line. The BLAST+ user manual can be found at https://www.ncbi.nlm.nih.gov/books/NBK279690/ and can be installed to your system by following the Installation tab.

Database Creation: Read this step if you wish to create your own custom database; otherwise, you may download existing databases from the web and skip to the next step. Begin by retrieving your files that contains your relevant information; these may be downloaded from sites as UniProt. Next, mask your file by using a masking program that came with your BLAST+ package; we recommend segmasker. While in the directory of your data files, type: 

segmasker -in (your_database) -infmt (your_database_format) -parse_seqids -outfmt (your_output_database_format, recommended maskinfo_asn1_bin) -out (your_database_name.asnb, name you wish to provide) 

Then, create your database by typing: 

makeblastdb -in (your_database) -input_type (your_data_file_type, remove this arg for fasta database) -dbtype (molecule_type_of_db, either ‘prot’ or ‘nucl’) -parse_seqids -mask_data (your_masked_data_file, the one created by segmasker) -out (your_output_database, name you wish to provide) -title (your_database_title, name you wish to provide) 

Your custom database is now successfully made. If ever stuck, please refer to the -h and -help commands for segmasker and makeblastdb, as well as to the user manual for further information.

Running BLAST+: The algorithm has compatibility for five types of BLAST, including blastp, blastn, blastx, tblastn, tblastx. Therefore, you may execute any of programs listed and obtain a result that could be run through the best-hit algorithm. To run BLAST simply, type your BLAST command followed by the four required arguments. For example, if you wished to test a protein sequence against a protein database, you would write:

blastp -db (your_database_name, made at previous section) -query (your_query_sequences_to_test.fasta) -out (your_blast_output_matches.txt, name you wish to provide) -outfmt “7 qseqid sseqid stitle sacc qstart qend sstart send qseq sseq evalue bitscore score length pident nident mismatch positive gapopen gaps ppos qlen slen”

The keywords “qseqid, sseqid, stitle, sacc, evalue, bitscore, length, pident, gaps, qlen” are mandatory for the best-hit algorithm to work at full capacity; you may find other keywords using the -help command and add them to -outfmt to receive additional information about your query. The number 7 indicates that the output will be a tabular text file with comment lines, which is the necessary file format for the BLAST best-hit algorithm. Other arguments may be added to the BLAST run to specialize your search. For example, you may add -matrix pam30 -threshhold 10 to a blastp search to specify a scoring matrix and a custom threshhold. Additional arguments can be found using the -h or -help commands, or at https://www.ncbi.nlm.nih.gov/books/NBK279675/. 

For a blastp search, we recommend you add “-matrix pam30 -comp_based_stats F” to the above to get the most optimized results.

*MAKE SURE TO SELECT CORRECT MATRIX, THRESHOLD, AND OTHER SETTINGS SO THAT BLAST+ RETRIEVES YOUR DESIRED RESULTS* 

Running the Algorithm 

You have now generated a tabular file that is ready to be processed for best-hits. To run the algorithm on your file, type: 

python BLASTCustomOutputScript.py -i (your_blast_output_matches.txt, file generated from running BLAST+) -s (species preference, type in your desired species in order (all caps) separated by commas; i.e. HUMAN,RAT,MOUSE)

The resulting .csv file should display the best match(es) for each species specified of every query sequence. You will notice that protein ID’s are listed in their respective species columns, but only one top match is displayed in full extent. Which top match is displayed depends primarily on which species the user typed first under the -s command, giving that species the highest priority of the group. The queries are also grouped into clusters; the queries in each cluster are highly similar and should often contain similar top picks. 

As of the moment, the algorithm processes every query independently and displays top picks independently of other queries in the cluster. This leaves a significant room for error if the user is not careful and only looks at the subject id column to generalize the best singular match for the whole cluster. To avoid this, please refer to the entire list of matches under the top pick column to determine the match you best see fit for each one of your queries. *Such uniform clustering analysis is currently in the process of being added to the program.*

The program’s optional arguments allow for broad customization and can be found under the -h command when running the program. They provide the user with options to specify an output file, to set values and cutoffs for certain match data bits, and to choose between the order of the filtering methods.

Best-Hit Logic

The following section describes the logic behind the best-hit algorithm. The program is split into a sequence of 10 filters. Each filter has the task of removing matches that have less desired databit values than those of other matches. Thus, one by one, the filters shrink the pool of matches for each query entry until the top picks are finally revealed. The general order of the filters is as follows: 

1) User Preferences: The program examines the user’s specified thresholds for percent identity, alignment length, and evalue, and throws out any matches that do no meet the user’s criteria. 

2) Species: The matches are split by species, and each group goes through the rest of the algorithm independently so as to produce a best match for each species indicated.

3) Identity: The matches are filtered by a percent identity and percent coverage filter. Percent coverage is a measure defined by the alignment length (minus gaps) over the query length. The two metrics are combined to provide a wholistic score. This filter is “soft” in a sense that matches up to 20% below the match with the best wholistic score are not filtered out. 

4) Bitscore: The matches go through a bitscore filter. Likewise, it is a “soft” filter in that matches up to 20% below the match with the best bitscore are not removed from the group.

5) E-Value: The matches are filtered by e-value. It is also a “soft” filter in that it allows matches with worse e-values than the match with the best e-value to remain; this threshold increases with decreasing value, up to 1.0e10 for very small value comparisons.

The purpose of these soft filters is to give more leverage between matches and reduce the probability of prematurely removing significant matches. The same would not be accomplished with “strict” or “hard” filters.

6) Review Status: The matches are sorted by their reviewed/unreviewed status, so that either all unreviewed matches are removed or all reviewed matches are removed. The latter case only happens if there exist unreviewed matches with higher identity scores than the identity score of the top reviewed match. This can be countered by specifying a mismatch tolerance (-m) that allows reviewed matches to trump unreviewed matches even if they contain a number of mismatches (which is reflected in identity score).

7) Ranking: The matches go through a strict identity filter to remove matches earlier approved by the soft identity filter.

8) Annotation: The matches are filtered by their titles and alphabetic IDs. Titles without “uncharacterized”,”putative”,”hypothetical” as well as IDs that begin with letters “P”,”O”,”Q” are given preference, and the other matches are filtered accordingly.

9) Isoform: (For Proteins) The matches go through an isoform filter that is meant to give preference to canonical matches over isoforms of those proteins. Furthermore, isoforms are reported as isoforms only when one isoform of a protein exists (in the match group) with an identity score of 100%; otherwise, the mentions of the match being an isoform are removed. 

10) Top Pick: The matches are processed through a hard evalue filter that helps make the final decision between near-identical matches. Note: multiple matches may still remain after this point, yet it will be very hard to distinguish among them.

This concludes the description of the filters and the general logic behind the algorithm.

Additional notes:
Two queries that share equivalent titles (in their .fasta files) will be considered as unique sequences and will be presented at different lines in the resulting file. Likewise, duplicate query sequences will not be merged as they will also be treated as different queries. 

Remember to update your databases every month or so by re-downloading them from online and creating them as discussed above. 

For further inquiries, feel free to contact me at ohad.koronyo@gmail.org
