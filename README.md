# CKMeansTreeClustering
A new fast method for building multiple consensus trees using k-means

# About
	=> =============================================================================================================
	=> Program   : KMeansTreeClustering - 2017
	=> Authors   : Nadia Tahiri and Vladimir Makarenkov (Universite du Quebec a Montreal)
	=> This program computes a clustering of phylogenetic trees based on the k-means partitioning algorithm.
	=> The set of trees in the Newick format should be provided as input.
	=> The optimal partitioning in K classes is returned as output. The number of classes can be determined by the 
	=> Silhouette and Calinski-Harabasz cluster validity indices adapted for tree clustering. The non-squared and 
	=> squared Robinson and Foulds topological distance can be used. 
	=> The recommended option: Silhouette + non-squared Robinson and Foulds distance.
	=> =============================================================================================================

# Installation
	$ git clone https://github.com/TahiriNadia/CKMeansTreeClustering.git
	$ make 
	or
	$ make install
	    
	clean project
	$ make clean
	 
# Help
	$ make help

# Examples
	Please execute the following command line:
	=> For trees: ./KMTC -tree input_file criterion
	
	Input_file - the input file for the program 
	criterion - the criterion for the k-means algorithm (1, 2, 3 or 4, see below)
	
	Command line execution:
	./KMTC -tree ../input/input.tre 3
	
# Input
	=> See the folder "data"
	Phylogenetic trees in the Netwick format (see the example in: input/input.tre)
	
# Output
	=> See the folder "output"
	The output is in the following files:
	1) stat.csv - for clustering statistics;
	2) output.txt - for cluster content.
