// K-means clustering using the SAS Fastclust algorithm,
// described on p. 352 of Numerical Ecology (1998).
// 
// // hoice of transformations for species abundance data,
// from Legendre and Gallagher (submitted).
// 
// Loop with different random starting configurations.
// 
// References:
// 
// Legendre, P. and E. Gallagher. 2001. Ecologically meaningful transformations 
// for ordination of species data. Oecologia (in press).
// 
//  Legendre, P. and L. Legendre, 1998. Numerical Ecology. 2nd English edition.
//  Elsevier Science BV, Amsterdam.
// 
//  Milligan, G. W. and M. C. Cooper. 1988. A study of standardization of 
//  variables in cluster analysis. Journal of Classification 5: 181-204.
// 
//                                           Pierre Legendre, August 1999
// 
//  345678901234567890123456789012345678901234567890123456789012345678901234567890
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <map>
#include <sstream>
#include <limits>


FILE *Output4;

//--Read the data
void ReadData1(int &n,int &nmax,int &p,int &pmax,double** mat/* ,double* coord */,int* ishort,double* weight,double* colsum,int &ntran,char* nameb,int N);
//--Calculate the kmeans
void CheckCentr(int &n,int &nmax,int &p,int &pmax,int &k1,int &kmax,double** mat,double** xbar,int* ishort,double** sx,int &idebug);
void Assign(int &iran,int &n,int &nmax,int &k1,int &kmax,int* list,int* howmany,int* no,int &idebug,int &iassign,int &iseed, int random_number);
//--Squared distances to group centroids. Assign objects to nearest one
void Centroids(int &n,int &nmax,int &p,int &pmax,int &kmax,double** mat,double** sx,double** xbar,double* mean,int* list,int* howmany,int &kk,int* ishort,int &idebug);
void CompSST(int &n,int &nmax,int &p,int &pmax,double** mat,double* weight,int* ishort,double &SST);
void Distances(int &n,int &nmax,int &p,int &pmax,int &kmax,double** mat,double** xbar,double* Dvec,int* list,int* howmany,double &SSE,int &kk,double* weight,int* ishort,int &idebug, double &W, double ** n_identique);
double DistanceBH(int &n,int &nmax,int &p,int &pmax,int &kmax,double** mat,double** xbar, int* howmany, int &kk,int* list,double* weight,int* ishort,map<int,string>  mapIndicesTreesFinal);
void Standard(int &n,int &nmax,int &kmax,int &p,int &pmax,double** mat,double** sx,double** sx2,double* vect,double** var,double* weight,int* ishort,int &istand);		
void Transform(int &n,int &nmax,int &p,int &pmax,double** mat,double* weight,double* colsum,int &ntran);
void Permute(int &iseed,int &n,int &nmax,int *iordre);
void Euclidean(int &i,int &k,int &nmax,int &kmax,int &p,int &pmax,double** mat,double** xbar,double* weight,int* ishort,double &D1);
void dist(int &i,int &k,double** mat,double &D1,int &n,int* list);
double Euclidean_point(int &i, int&j, int &p,int &pmax,double** mat,double* weight,int* ishort);
int Split(vector<string>& vecteur, string chaine, char separateur);
void listePartition(int Strouve [],map<int,string>  mapIndicesTreesFinal,vector <string> indicesTrees);
int fact (int valeur);
double f_RI(int Strouve[],int Sref[],int N);
double f_ARI(int Strouve[],int Sref[],const char *K_real,int group,int N);
void outStat(int Strouve[],int Sref[],char *criteria,int N,char *N_especes,char *percent,const char *K_real,int group,double score,int **listr);
void convert_list(int *&list, int n, vector<string> &centroid_k,vector<int> &centroid_k_pos, int &kk);
void RF_tree2tree(double &D1,string tree_i_k,string centro_k);
double NouvelleCallRF(double nb_esp_identique);
void delete_space(string &str);
double DistanceCH_old(int &n,int &kmax,double** mat,int* list,double** Ww);

//fonctions sans passer par le centroid
double DistanceCH(int &n,int &kmax,double** mat,int* list,double** Ww,double FO_new, double facteur);

//Borne INFERIEUR
double FO_borne_inf(int &n,int &kmax,double** mat,double* Dvec,int* list,int* howmany,double &SSE,int &kk,vector <string> monTableau, double ** n_identique);

//Borne MOYENNE
double FO_borne_avg(int &n,int &kmax,double** mat,double* Dvec,int* list,int* howmany,double &SSE,int &kk,vector <string> monTableau, double ** n_identique);

//Borne SUPERIEUR
double FO_borne_supp(int &n,int &kmax,double** mat,double* Dvec,int* list,int* howmany,double &SSE,int &kk,vector <string> monTableau);

//SUPERTREES
double FO_super_tree(int &n,int &kmax,double** mat,double* Dvec,int* list,int* howmany,double &SSE,int &kk,vector <string> monTableau, double ** n_identique, double ** tree_cluster_leaves);
//fonctions sans passer par le centroid pour supertree
double DistanceCH_supertree(int &n,int &kmax,double** mat,int* list,double** Ww,double FO_new, double facteur, double ** tree_cluster_leaves, double ** n_identique);

//Borne EUCLIDEAN
double FO_euclidean(int &n,int &kmax,double** mat,double* Dvec,int* list,int* howmany,double &SSE,int &kk,vector <string> monTableau, double ** n_identique);

//fonctions en passant par le pseudo-centroid
//BORNE SUPERIEUR (similaire a K-medoid)
void Centroids_approx(int &n,double** mat,double** xbar,int* list,int &kk,vector<int> &centroid_k_pos);
double FO_approx_Kmedoid(int &n,int &kmax,double** mat,double* Dvec,int* list,int* howmany,double &SSE,int &kk,vector<int> &centroid_k_pos, double ** n_identique);
double DistanceCH_approx_Kmedoid(int &n,int &kk,double** mat,vector<int> &centroid_k_pos,double FO_new,int indice_arbre_cons,int* list);

//SILHOUETTE INDEX
double DistanceSilhouette(int &n,int &nmax,int &p,int &pmax,int &kmax,double** mat,double** xbar,int* howmany, int &kk,int* list,double* weight,int* ishort);
double FO_silhouette(int &n,int &kmax,double** mat,double* Dvec,int* list,int* howmany,double &SSE,int &kk,vector<int> &centroid_k_pos);

//fonctions en passant par le centroid (consense)
void Centroids_consensus(int &n,double** mat,double** xbar,int* list,int &kk, vector<string> &centroid_k, vector<string> monTableau);
double FO_consense(int &n,int &kk,vector<string> &centroid_k, vector<string> monTableau,double** mat,double* Dvec,int* list,int* howmany,double &SSE,vector<int> &centroid_k_pos,char *N_especes);
double DistanceCH_consense(int &n,int &kk,vector<string> &centroid_k,string centroid_C_min,double FO_new);

void Distances_approx(int &n,int &nmax,int &p,int &pmax,int &kmax,double** mat,double** xbar,double* Dvec,int* list,int* howmany,double &SSE,int &kk,double* weight,int* ishort,int &idebug, double &W,vector <string> monTableau,vector<string> &centroid_k,vector<int> &centroid_k_pos);
void DFS(int i, int *visited, int **mat_adjacence, int n, int &nb_visited);
void buildMatAdjecence(int k, int n, int *list, double **n_identique, int **mat_adjacence);
void buildMatAdjacenceFusion(int k1, int k2, int n, int *list, double **n_identique, int **mat_adjacence);
double arrondir(double num,int digits);
void conv2sameRef(int *Strouve,int *Sref, int N);

//fonction W
double FO_W(int &n,int &kmax,double** mat,double* Dvec,int* list,int* howmany,double &SSE,int &kk,vector <string> monTableau);
double DistanceW(int &n,int &kmax,double** mat,int* list,double** Ww,double FO_new,double facteur, double ** n_identique);
double DistanceW2(int &n,int &kmax,double** mat, int* list, double** Ww);

//fonctions sans passer par le centroid -- GAP STATISTIQUE
double DistanceGAP(int &n,int &kmax,double** mat,int* list,double** Ww,double FO_new, double facteur, double ** n_identique);

//GAP STATISTIQUE
double FO_gap_statistique(int &n,int &kmax,double** mat,double* Dvec,int* list,int* howmany,double &SSE,int &kk,vector <string> monTableau);


//  Parameters of the program:
//  n    = number of observations (ex. number of trees)   
//  nmax = maximum number of observations (default: 10000)
//  p    = number of variables     (ex. number of variables describing each trees)
//  pmax = maximum number of variables  (default: 10000)  
//  k    = number of groups (centroids)
//  kmax = maximum number of groups
//  kk   = ???
//  niter = maximum iteration for convergeance of centroid (fixed=100)
//  Parameter (nmax=100000,pmax=10,kmax=10)
//  critera = (0,1,2)
//            0: C-H  
//            1: logSS  
//            2: Silhouette   
//			  3: W
//   mat       =  double[nmax+1][pmax+1] (data matrix))
//   weight         =  double vector[p] of weigth vector for variable p
//   xbar      =  vector [k][i] = mean distance of i to centroid k
//   list      =  vector[i]=k contains group assignments for the objects.
//   howmany[k]=  contains the number of objects in each group centroid (k). 
//   ishort[p] = vector containing the position(p) of valide variables (containing non zero value for weight -- see ReadData). 
//   nobest
//   var       = double [kmax+1]
//   vect      = double [pmax+1]
//   mean[p]      = vector of mean value not weigthed for the variable [p] double [pmax+1] 
//   Dvec      = double [kmax+1]
//   CHr       = vector[k] of the best CH this k
//   BHr       = vector[k] of the best BH this k
//   Silr      = vector[k] of the best Silr this k
//   LogSSr    = vector[k] of the best LogSS this k
//   Wr        = vector[k] of the best W this k
//   SSEr      = double [kmax+1]  sum of squared error statistic (SSE) = within-group sum of squares
//   listr     = int[kmax+1][nmax+1];
//   howmanyr  = int[kmax+1][nmax+1];	
// sx = new double*[kmax+1];
//  sx2 = new double*[kmax+1];

//   no = new int [nmax+1]; 
//	iordre = new int [nmax+1];
//   nobest = new int [kmax+1];
//   nnitr = new int [kmax+1];
//
// Output file for summary results	
//
// Output2   = file.statistics.txt
// Output3   = file.sum.txt
// Output4   = stat.txt
// Output7    = file.groups.txt
//
// Note, we expect the input file to be a matrix in the form (below))
// were 5 is the n, and 4 is the p and the following line are the value (double))
//5	4
//0.0	0.0	0.0	1.0	0.0	
//0.0	1.0	0.0	0.0	1.0	
//1.0	1.0	1.0	0.0	1.0	
//0.0	0.0	0.0	0.0	1.0	
//
// IMPORTANT: See also the variables below

//Variables globales
/* map<int,int> nbNi; */


int main_kmeans(char **argv,vector <string> monTableau, double ** mat, double ** n_identique,double ** Ww, vector<int> tabIndices, int intParam)
{
	//*****************Define variables******************************************//
	 // Variables
    map<int,string>  mapIndicesTreesFinal;
	vector <string> indicesTrees;
	time_t tbegin2,tend2;
    double texec2=0.;
	 
	double W = 0.0;
	double CH = -10000000000.0;
	double silhouette=0.0;
	double LogSS=0.0;
	double BH=0.0;
	double GapStat=0.0;

	bool use_weight=false;

	double CHr_max=-10000000000.0;
	int CHr_group=0;

	double BHr_min=1000000.0;
	int BHr_group=0;

	double Silr_max=-10000000000.0;
	int Silr_group=0;

	double LogSS_min=1000000.0;
	int LogSS_group=0;

	double W_min=1000000.0;
	double W_max=-10000000000;
	int W_group=0;
	double FO_new = 10000000000.0;
	
	double GAP = -10000000000.0;
	double Gapr_max=-10000000000.0;
	int Gapr_group=0;

    // Start timer
    tbegin2=time(NULL);				// get the current calendar time
	
	int N = monTableau.size(); //quantity of initial tree
	int i=0, j=0;		//Counters
	int n=N, p=N;		//,pmax,kmax; //      Integer p,pmax,kmax
	int ntran=0, itemp=0, iseed=0, niter=0, kk=0, nit=0;		//added declarations for variables
	int nnit=0, k=0, i1ref=0, i2ref=0, igr1=0, igr2=0;		//added declarations for variables
	int idebug=0 ; // 0, no debug, 1 debug
	int k1=0, k2=0;  //added declarations for variables
	int hard_max_k=0; //--Setting the max k1
        
    int random_number=5; //--Fixed random number
	int istand=0;   //--0 No standardization
	int iassign=2;  // 1 equal, 2 random
    int iran=5;   //--Number of random position
    int nran=100;  //--Number of Random start  
		
	int nmax=N;    //--Maximum number of object -Parameter (nmax=10000,pmax=250,kmax=100)
	int pmax=N;      //--Maximum data point (variable))
 	int kmax=N;	  // Maximum number of groups
	
	char *criteria = argv[0];
	char *N_especes = argv[1];
	const char *K_real = argv[2];
	char *percent = argv[3];
	
	int Strouve[N];
	
	for(int linej=0;linej<N;linej++){
		Strouve[linej]= 0;
	}
	
	double	**sx,**sx2,**xbar,**var;	//sx(kmax,pmax),sx2(kmax,pmax),xbar(kmax,pmax),var(kmax,pmax)
	double **tree_cluster_leaves = new double *[n];
	for(int i=0;i<n;i++){
		tree_cluster_leaves[i]=new double [4];
	}
	
	sx = new double*[kmax+1];
	sx2 = new double*[kmax+1];
	xbar = new double*[kmax+1];
	var = new double*[kmax+1];

	for (i=0;i<=kmax;i++)
	{
		sx[i] = new double[pmax+1];
		sx2[i] = new double[pmax+1];
		xbar[i] = new double[pmax+1];
		var[i] = new double[pmax+1];
	}

	//variables centroids (trees and indices)
	vector <string> centroid_k; 
	vector <int> centroid_k_pos; 
	string centroid_C_min = "";
	
	double *distances_RF_norm = new double[4];
	
	for(int linej=0;linej<4;linej++){
		distances_RF_norm[linej]= 0.0;
	}
	
	for (i=0; i<=kmax; i++)
	{
		for (j=0; j<=pmax; j++)
		{
			sx[i][j] = 0.0;
			sx2[i][j] = 0.0;
			xbar[i][j] = 0.0;
			var[i][j] = 0.0;
		}
	}
	

	double *Dvec,*CHr, *BHr,*SSEr, *Silr, *LogSSr, *Wr, *diff_W, *V_W, *Wr_ln, *Gapr;
	Dvec = new double [kmax+1];
	SSEr = new double [kmax+1];
	
	CHr = new double [kmax+1];
	Gapr = new double [kmax+1];
	BHr = new double [kmax+1];
	Silr = new double [kmax+1];
	LogSSr = new double [kmax+1];
	Wr = new double [kmax+1];
	Wr_ln = new double [kmax+1];
	diff_W = new double [kmax+1];
	V_W = new double [kmax+1];
	

	for (i=0; i<=kmax; i++)
	{
		Dvec[i] = 0.0;
		SSEr[i] = 0.0;		
	}

	double *vect,*mean,*weight;		//vect(pmax),mean(pmax),weight(pmax),
	vect = new double [pmax+1];
	mean = new double [pmax+1];
	weight = new double [pmax+1];
	for (i=0; i<=pmax; i++)
	{
		vect[i] = 0.0;
		mean[i] = 0.0;
		weight[i] = 0.0;
	}


	double D1=0,Dref=0,SSE=0,SSEref=0,SST=0;		
    double temp=0;					

	int **listr;					//listr(kmax,nmax),
	listr = new int*[kmax+1];
	for (i=0;i<=kmax;i++)
	{
		listr[i] = new int [nmax+1];
	}

	for (i=0; i<=kmax; i++)
	{
		for (j=0; j<=nmax; j++)
		{
		  listr[i][j] = 1;
		}
	}


	int	**howmanyr;		//howmanyr(kmax,kmax)
	howmanyr = new int*[kmax+1];
	for (i=0;i<=kmax;i++)
	{
		howmanyr[i] = new int [kmax+1];
	}

	for (i=0; i<=kmax; i++)
	{
		for (j=0; j<=kmax; j++)
		{
		  howmanyr[i][j] = 0;
		}
	}

	
	int *list,*no, *iordre;		//list(nmax),no(nmax), iordre(nmax)
	list = new int [nmax+1];
	no = new int [nmax+1];
	iordre = new int [nmax+1];

	for (i=0; i<=nmax; i++)
	{
		list[i] = 0;
		no[i] = 0;
		iordre[i] = 0;
	}

	int *howmany,*nobest/*, *nobestSilhouette, *nobestLogSS ,*nobestCH, *nobestBH, *nobestW */, *nnitr;		//howmany(kmax),,nobest(kmax), nnitr(kmax);	
	howmany = new int [kmax+1];
	nobest = new int [kmax+1];
    nnitr = new int [kmax+1];
	
	for (i=0; i<=kmax; i++)
	{
		howmany[i] = 0;
		nobest[i] = 0;
		nnitr[i] = 0;
	}

	int *ishort;			//ishort(pmax);										
	ishort = new int [pmax+1];
	
	for (i=0; i<=pmax; i++)
	{
		ishort[i] = 0;
	}


//  Modification Centroids: add ",nameb" to next line
	char *nameb; 
	nameb = new char [300];
        

//***********************  Read data file  **********************************
         
   	int max_k1 = 5;
	
	if (N<=5) max_k1=N-1;
	
	k1=max_k1;
	double facteur = 1.0;
	k2=2;
	if(intParam==1 || intParam==6 || intParam==8 || intParam==9 || intParam==10 || intParam==11 || intParam==12 || intParam==14 || intParam==20){
		if(intParam==11){
			facteur = (2.0/(N*1.0));
		}else if(intParam==9){
			facteur = (1.0/(N-1.0));
		}else if(intParam==12){
			facteur = (1.0/(N*1.0));
		}else if(intParam==10){
			facteur = ((3.0*N)-2.0)/(2.0*N*(N-1.0));
		}
		
		for (i=0; i<=kmax; i++){
			CHr[i] = -10000000000.0;
		}
		
		if (k1<=2) {
			printf("*** Warning, not enough trees (k1:%d) k1 set to 3\n",k1);
			k1=max_k1;      
		}
	}else if(intParam==2){
		for (i=0; i<=kmax; i++){
			BHr[i] = 0.0;
		}
		
		if (k1<=2) {
			printf("*** Warning, not enough trees (k1:%d) k1 set to 3\n",k1);
			k1=max_k1;      
		}
	}else if(intParam==3 || intParam==16 || intParam==17 || intParam==18 || intParam==19){
		for (i=0; i<=kmax; i++){
			Silr[i] = -10000000000.0;
		}
		
		if (k1<=2) {
			printf("*** Warning, not enough trees (k1:%d) k1 set to 3\n",k1);
			k1=max_k1;      
		}
	}else if(intParam==4){
		for (i=0; i<=kmax; i++){
			LogSSr[i] = 0;
		}
		
		if (k1<=2) {
			printf("*** Warning, not enough trees (k1:%d) k1 set to 3\n",k1);
			k1=max_k1;      
		}
	}else if(intParam==5){
		for (i=0; i<=kmax; i++){
			Wr[i] = 1000000000.0;
			Wr_ln[i] = -10000000000.0;
			diff_W[i] = 0.0;
			V_W[i] = 0.0;
		}
		k2=1;
		//k2=1;
		if (k1<=1) {
			printf("*** Warning, not enough trees (k1:%d) k1 set to 3\n",k1);
			k1=max_k1;      
		}
	}else if(intParam==13 || intParam==15){
		for(int dist1=0; dist1<n;dist1++){
			for(int dist2=0; dist2<n;dist2++){
				mat[dist1][dist2]=mat[dist1][dist2]*mat[dist1][dist2];
			}
		}
		for (i=0; i<=kmax; i++)
		{
			CHr[i] = -10000000000.0;
		}

		if (k1<=2) {
			printf("*** Warning, not enough trees (k1:%d) k1 set to 3\n",k1);
			k1=max_k1;      
		}
	}else if(intParam==7){
		for(int dist1=0; dist1<n;dist1++){
			for(int dist2=0; dist2<n;dist2++){
				mat[dist1][dist2]=mat[dist1][dist2]*mat[dist1][dist2];
			}
		}
		
		for (i=0; i<=kmax; i++)
		{
			Wr[i] = 1000000000.0;
			Wr_ln[i] = -10000000000.0;
			diff_W[i] = 0.0;
			V_W[i] = 0.0;
		}
		
		//k2=1;
		if (k1<=1) {
			printf("*** Warning, not enough trees (k1:%d) k1 set to 3\n",k1);
			k1=max_k1;      
		}
	}else if (intParam==21){
		for (i=0; i<=kmax; i++){
			Gapr[i] = -10000000000.0;
		}
		
		k2=1;
		if (k1<=2) {
			printf("*** Warning, not enough trees (k1:%d) k1 set to 3\n",k1);
			k1=max_k1;      
		}
	}
		
	if (k1>kmax) {
	  printf("*** Warning, limiting groups to %d \n",kmax);
	  k1=max_k1-1;          
	}

	if (hard_max_k!=0) {
	  k1=max_k1;          
	}

	 //--Read the data from files
    ReadData1(n,nmax,p,pmax,mat/* ,coord */,ishort,weight,mean,ntran,nameb,N);	//Call ReadData11(n,nmax,p,pmax,mat,coord,ishort,w,mean,ntran,namea)
	     
    CompSST(n,nmax,p,pmax,mat,weight,ishort,SST);
	
	double dist_min = 10000000000.0;
	double dist_i = 0.0;
	int indice_arbre_cons = 0;
	
	if(intParam==1 || intParam==3 || intParam==14){
		for (int i=0;i<N;i++){
			dist_i = 0.0;
			for (int j=0;j<N;j++){
				dist_i += mat[i][j];	
			}	
			
			if(dist_i<dist_min){
				centroid_C_min = monTableau[i];
				indice_arbre_cons = i+1;
				dist_min = dist_i;		
			}	
		}	
	}else if(intParam==8 || intParam==15){
		// Ouvertute du fichier contenant les parametres du programme Consense
		
		ofstream myfileTree;
		myfileTree.open ("trees.txt");
		
		for(int n_tree=0; n_tree<monTableau.size(); n_tree++){
			myfileTree << monTableau[n_tree];
			myfileTree << endl;
		}
		myfileTree.close();
		
		ofstream myfile;
		myfile.open ("parametresConsense");
		myfile << "trees.txt";
		myfile << "\n";
		myfile << "C\n";
		myfile << "C\n";
		myfile << "Y\n";
		myfile.close();
	 
		
		//appel du logiciel consense
		system("./consense <parametresConsense >outTmp");
		
		//Format consensus tree C for only one line
		system("cat outtree|tr \"\n\" \" \">tmp.txt");
		system("cp tmp.txt outtree");

		//Recuperer le string du consensus des arbres de départ par la variable centroid_C_min
		ifstream fileC;
		fileC.open("outtree");
		
		getline(fileC,centroid_C_min);
		fileC.close(); 
		
		system("rm outtree outfile outTmp trees.txt parametresConsense");
	}
	
	for(int i1=0; i1<n; i1++){
		for(int i2=0; i2<n; i2++){
			mat[i1][i2] = arrondir(mat[i1][i2],3);
		}
	}
	
	// Compute vector 'mean' of overall means
	for (j=1;j<=p;j++)  { mean[j]=0;  }

	for (i=1;i<=n;i++) 
	{
	  for (j=1;j<=p;j++) 
	  {
		mean[j]=mean[j]+mat[i-1][ishort[j]-1];
	  }
	}

	for (j=1;j<=p;j++)
	{
		mean[j]=mean[j]/(n*1.0);
	}
       

	iseed=0;
	double CH_new = -10000000000.0;
	double Gap_new = -10000000000.0;
	double W_new = 1000000000.0;
	double SH_new = -10000000000.0;
	
	//initialize Sref
	int Sref [N];
	for (int j=0; j<N; j++){
		Sref[j]=0;
	}
	
	int number_cluster = 0;
	int k_cluster = 0;
	 
	int nbInit =0;
	int nbFin =0;
	map <int, int> CH_conversion;
	map <int, int> Gap_conversion;
	map <int, int> BH_conversion;
	map <int, int> Silhouette_conversion;
	map <int, int> LogSS_conversion;	
	map <int, int> W_conversion;

	int realk = 0;
	int unique = 0;
	int CHk = 0;
	int BHk = 0;
	int Gapk = 0;
	int Silhouettek = 0;
	int LogSSk = 0;
	int wk = 0;
	for (i=0; i<=kmax; i++){
		howmany[i] = 0;
		nobest[i] = 0;
		nnitr[i] = 0;
	}
	
	int **mat_adjacence = new int *[n];
	
	for(int i=0;i<n;i++){
		mat_adjacence[i]=new int [n];
	}
	//Boucle interne de la recherche du pivot de matrice d'adjecence (celui ne contient pas 2)
	int aret_de_recherche = 0;
	int depart_graphe_i = 0;
	int depart_graphe_j = 0;
	
	int * visited = new int[n];//indique si les elements ont ete visites
	int nb_visited = 0; //compte le nombre d'element visite 
	int *nk = new int [kmax+1];
	int bon_partition_initiale = 0;
	int nb_partition_OK = 0;
	int max_possibility_parti_init = 0;
	
	//initialisation des variables
	for(int k=1;k<=kmax; k++){
		nk[k]=0;
	}
	
    for (iran=1;iran<=nran; iran++) {
		CH_new = -10000000000.0;
		Gap_new = -10000000000.0;
		SH_new = -10000000000.0;
		silhouette = -10000000000.0;
		BHk = 0;
		Silhouettek = 0;
		CHk = 0;
		Gapk = 0;
		use_weight = false;
		LogSSk = 0;
		wk = 0;
		realk = 0;
		unique = 0;
		number_cluster = 0;

		if(intParam==20){
			if(iassign!=4){		
				bon_partition_initiale=0;
				while(bon_partition_initiale==0 && max_possibility_parti_init<10){
					nb_partition_OK = 0;	
					//initialisation des variables
					for(int k=1;k<=kmax; k++){
						nk[k]=0;
					}	
					Assign(iran,n,nmax,k1,kmax,list,howmany,no,idebug,iassign,iseed, random_number);
					//Compter le nombre d'éléments par cluster
					for(int k=1;k<=kmax; k++){
						nk[list[k]]++;
					}

					for(int k_init=1;k_init<=k1;k_init++){
						nb_visited=0;	
						//Inferer la matrice d'adjacence du cluster k_init
						buildMatAdjecence(k_init,n,list,n_identique,mat_adjacence);						
						//visited is initialized to zero
						for(int i=0;i<n;i++)
							visited[i]=0;
									
						//connaitre la premiere ligne et la premiere (la premiere cellule)
						//qui ne contient pas un 2
						aret_de_recherche = 0;
						depart_graphe_i = 0;
						while(depart_graphe_i<n && aret_de_recherche==0){
							depart_graphe_j = 0;
							while(depart_graphe_j<n && aret_de_recherche==0){
								if(mat_adjacence[depart_graphe_i][depart_graphe_j]!=2){
									aret_de_recherche = 1;
								}
								depart_graphe_j++;
							}
							depart_graphe_i++;
						}
						DFS(depart_graphe_i-1,visited,mat_adjacence,n,nb_visited);
						if(nb_visited==nk[k_init]){
							nb_partition_OK++;
						}
					}
									
					if(nb_partition_OK==k1){
						bon_partition_initiale=1;
					}
					max_possibility_parti_init++;	
				}
				if(bon_partition_initiale!=1){
					k1=k2-1;
				}
			}
		}else{
			if(iassign!=4){
				Assign(iran,n,nmax,k1,kmax,list,howmany,no,idebug,iassign,iseed, random_number);
			}
		}
        // Big loop on number of groups, downwards from k1 to k2 (k1>=k2) - - - - - -
		niter=100;

		for (kk=k1;kk>=k2;kk--){
			SSEref=1.0e20;		
			W = 1000000000.0;
			CH = -10000000000.0;
			GAP = -10000000000.0;
			LogSS=0.0;
			BH=0.0;
			FO_new = 10000000000.0;
			W_new = 10000000000.0;

			for (nit=1;nit<=niter;nit++) 
			{    
				
				if(idebug==1){
					printf ("Iteration = %d",nit);	
					printf ("SSEref = %lf",SSEref);	
					for (i=1;i<=n;i++){
						printf ("%lf  ",list[i]);	
					}
				}					
				nnit=nit;			

				// If centroids are read, we want to use then during the first iteration
                // (nit-1) of the first level of the cascade (kk=k1).
				  if(!((iassign==4)&&(kk==k1)&&(nit==1))) 
				  {
                      // Compute centroids for all groups; write to matrix 'xbar'.
					  if(intParam==6 || intParam==7){
						Centroids(n,nmax,p,pmax,kmax,mat,sx,xbar,mean,list,howmany,kk,ishort,idebug);
					  }else if(intParam==8 || intParam==15){
						Centroids_consensus(n,mat,xbar,list,kk,centroid_k,monTableau);
					  }else if(intParam==1 || intParam==3){
						Centroids_approx(n,mat,xbar,list,kk,centroid_k_pos);
					  }
				  }	
				
			
                  // Compute distances to group centroids and assign objects to nearest one   
					if(intParam==6 || intParam==7){
						Distances(n,nmax,p,pmax,kmax,mat,xbar,Dvec,list,howmany,SSE,kk,weight,ishort,idebug, W, n_identique);//Call Distances(n,nmax,p,pmax,kmax,mat,xbar,Dvec,list,howmany,SSE,kk,w,ishort,idebug)					
					}else if(intParam==9 || intParam==16){
						FO_new = FO_borne_inf(n,kmax,mat,Dvec,list,howmany,SSE,kk,monTableau, n_identique);
					}else if(intParam==10 || intParam==17){
						FO_new = FO_borne_avg(n,kmax,mat,Dvec,list,howmany,SSE,kk,monTableau, n_identique);
					}else if(intParam==11 || intParam==18){
						FO_new = FO_borne_supp(n,kmax,mat,Dvec,list,howmany,SSE,kk,monTableau);
					}else if(intParam==1 || intParam==3){
						FO_new = FO_approx_Kmedoid(n,kmax,mat,Dvec,list,howmany,SSE,kk,centroid_k_pos, n_identique);
					}else if(intParam==8 || intParam==15){
						FO_new = FO_consense(n,kk,centroid_k,monTableau,mat,Dvec,list,howmany,SSE,centroid_k_pos,N_especes);
					}else if(intParam==12 || intParam==19){
						FO_new = FO_euclidean(n,kmax,mat,Dvec,list,howmany,SSE,kk,monTableau, n_identique);						
					}else if(intParam==20){
						FO_new = FO_super_tree(n,kmax,mat,Dvec,list,howmany,SSE,kk,monTableau, n_identique,tree_cluster_leaves);
					}else if(intParam==21){
						FO_new = FO_gap_statistique(n,kmax,mat,Dvec,list,howmany,SSE,kk,monTableau);						
					}else if(intParam==5){
						FO_new = FO_W(n,kmax,mat,Dvec,list,howmany,SSE,kk,monTableau);
					}  
				  
				  number_cluster = 0;
				  		  
					if(intParam==1){
						CH_new = DistanceCH_approx_Kmedoid(n,kk,mat,centroid_k_pos,FO_new,indice_arbre_cons,list);
						
						if(CH_new>=CHr[kk]){
							SSEr[kk]=SSE;
							nobest[kk]=iran;
						
							nnitr[kk]=nnit;	
							CH=CH_new;
							CHr[kk]=CH;
							for (int i=1;i<=n;i++)			
							{
								listr[kk][i]=list[i];
							}	
							
							for (int i=1;i<=kk;i++){
								howmanyr[kk][i]=howmany[i];
							}	
						}
					}else if(intParam==2){
						BH = DistanceBH( n,nmax,p,pmax,kmax,mat,xbar,howmany, kk,list,weight,ishort,mapIndicesTreesFinal);
					}else if(intParam==3 || intParam==16 || intParam==17 || intParam==18 || intParam==19){
						SH_new = DistanceSilhouette( n,nmax,p,pmax,kmax,mat,xbar,howmany, kk,list,weight,ishort);
						if(SH_new>=Silr[kk]){
							SSEr[kk]=SSE;		
							nobest[kk]=iran;
						
							nnitr[kk]=nnit;
							Silr[kk]=SH_new;
							for (int i=1;i<=n;i++){
								listr[kk][i]=list[i];
							}
							
							for (int i=1;i<=kk;i++){
								howmanyr[kk][i]=howmany[i];
							}
						}
					}else if(intParam==4){
						LogSS = log10((SST)/SSE);
						if(LogSS>=LogSSr[kk]){
							SSEr[kk]=SSE;		
							nobest[kk]=iran;	
						
							nnitr[kk]=nnit;	
							LogSSr[kk]=LogSS;
							for (int i=1;i<=n;i++)			
							{
								listr[kk][i]=list[i];
							}	
							
							for (int i=1;i<=kk;i++){
								howmanyr[kk][i]=howmany[i];
							}
						}
					}else if(intParam==5 || intParam==7){
						W_new = DistanceW(n,kmax,mat,list,Ww,FO_new,facteur,n_identique);
														
						if(W_new<=Wr[kk]){
							W=W_new;
							Wr[kk]=W;
							
							SSEr[kk]=SSE;		
							nobest[kk]=iran;	

							for (int i=1;i<=n;i++){
								listr[kk][i]=list[i];
							}
							
							for (int i=1;i<=kk;i++){
								howmanyr[kk][i]=howmany[i];
							}
						}
					}else if(intParam==6){
						CH_new = DistanceCH_old(n,kmax,mat,list,Ww);
													
						if(CH_new>=CHr[kk]){
							SSEr[kk]=SSE;		
							nobest[kk]=iran;
						
							nnitr[kk]=nnit;	
							CH=CH_new;
							CHr[kk]=CH;  
							for (int i=1;i<=n;i++){
								listr[kk][i]=list[i];
							}
															
							for (int i=1;i<=kk;i++){
								howmanyr[kk][i]=howmany[i];
							}
						}
					}else if(intParam==8 || intParam==15){
						CH_new = DistanceCH_consense(n,kk,centroid_k,centroid_C_min,FO_new);
						if(CH_new>=CHr[kk]){
							SSEr[kk]=SSE;	
							nobest[kk]=iran;	
						
							nnitr[kk]=nnit;
							CH=CH_new;
							CHr[kk]=CH;  
							for (int i=1;i<=n;i++)			
							{
								listr[kk][i]=list[i];
							}
							
							for (int i=1;i<=kk;i++){
								howmanyr[kk][i]=howmany[i];
							}	
						}
					}else if(intParam==9 || intParam==10 || intParam==11 || intParam==12){
						CH_new = DistanceCH(n,kmax,mat,list,Ww,FO_new,facteur);
						if(CH_new>=CHr[kk]){
							SSEr[kk]=SSE;	
							nobest[kk]=iran;
						
							nnitr[kk]=nnit;
							CH=CH_new;
							CHr[kk]=CH;  
							for (int i=1;i<=n;i++)
							{
								listr[kk][i]=list[i];
							}	
															
							for (int i=1;i<=kk;i++)
							  {howmanyr[kk][i]=howmany[i];}	
						}
					}else if(intParam==13){
						
						if(CH_new>=CHr[kk]){
							SSEr[kk]=SSE;
							nobest[kk]=iran;
						
							nnitr[kk]=nnit;	
							CH=CH_new;
							CHr[kk]=CH;  
							for (int i=1;i<=n;i++)			
							{
								listr[kk][i]=list[i];
							}
															
							for (int i=1;i<=kk;i++)			
							  {howmanyr[kk][i]=howmany[i];}	
						}
					}else if(intParam==20){
						CH_new = DistanceCH_supertree(n,kmax,mat,list,Ww,FO_new,facteur,tree_cluster_leaves,n_identique);
						if(CH_new>=CHr[kk]){
							SSEr[kk]=SSE;		
							nobest[kk]=iran;	
						
							nnitr[kk]=nnit;	
							CH=CH_new;
							CHr[kk]=CH; 
							for (int i=1;i<=n;i++)
							{
								listr[kk][i]=list[i];
							}
															
							for (int i=1;i<=kk;i++)	{
								howmanyr[kk][i]=howmany[i];
							}	
						}
					}else if(intParam==21){
						Gap_new = DistanceGAP(n,kmax,mat,list,Ww,FO_new,facteur,n_identique);
						
						if(Gap_new>=Gapr[kk]){
							SSEr[kk]=SSE;		
							nobest[kk]=iran;
						
							nnitr[kk]=nnit;
							GAP=Gap_new;
							Gapr[kk]=GAP;  
							for (int i=1;i<=n;i++)		
							{
								listr[kk][i]=list[i];
							}	
															
							for (int i=1;i<=kk;i++){
								howmanyr[kk][i]=howmany[i];
							}	
						}
					}
				  
				  // Compute sum of squared error statistic (SSE) = within-group sum of squares 
				 if(fabs(SSEref-SSE)>(SSE/1000.0))	
				  {
					  SSEref=SSE;
				  }
				  else			
				  {
					goto m60;  
				  }	
					
			  }	

              // Compute the Calinski-Harabasz (1974) index 'CH' and
m60:          
              // Concatenate the two closest groups before going to the next value of kk
			  Dref=1000000.0;	
			  D1=0.0;		
			  i1ref=1;		
			  i2ref=2;		
			  if(intParam==7){
				  for (igr1=1;igr1<kk;igr1++)
				  {
					  for (igr2=(igr1+1);igr2<=kk;igr2++)
					  {
						  D1=0.0;	
						  for (j=1;j<=p;j++)	
						  {
							  D1=D1+( (xbar[igr1][j]-xbar[igr2][j])*(xbar[igr1][j]-xbar[igr2][j]));	
						  } 
						  if(D1<Dref) 
						  {
							  i1ref=igr1;		
							  i2ref=igr2;		
							  Dref=D1;			
						  }	//endif
					  }
				  }
			  }else if(intParam==1 || intParam==3 || intParam==5 || intParam==6 || intParam==8 || intParam==9 ||intParam==10 || intParam==11 || intParam==12 || intParam==13 ||  intParam==15 || intParam==16 || intParam==17 || intParam==18 || intParam==19 || intParam==20 || intParam==21){
				if(intParam==20){ // Cas of supertree
					//Initialisation of variables
					int i_fusion = 1;
					int j_fusion = kk;
					int fusion_OK = 0;
				  
					//Loop until each partition is connected
					while (fusion_OK==0 && i_fusion<kk){
						j_fusion = kk;
						while (fusion_OK==0 && j_fusion>i_fusion){
							
							//Build adjacence matrix
							buildMatAdjacenceFusion(i_fusion,j_fusion,n,list,n_identique,mat_adjacence);
						
							//visited is initialized to zero
							nb_visited=0;	
							for(int i=0;i<n;i++)
								visited[i]=0;
									

							//connaitre la premiere ligne et la premiere (la premiere cellule)
							//qui ne contient pas un 2
							aret_de_recherche = 0;
							depart_graphe_i = 0;
							while(depart_graphe_i<n && aret_de_recherche==0){
								depart_graphe_j = 0;
								while(depart_graphe_j<n && aret_de_recherche==0){
									if(mat_adjacence[depart_graphe_i][depart_graphe_j]!=2){
										aret_de_recherche = 1;
									}
									depart_graphe_j++;
								}
								depart_graphe_i++;
							}
						
							//Depth-first search (Recherche en profondeur)
							DFS(depart_graphe_i-1,visited,mat_adjacence,n,nb_visited);
							if(nb_visited==(howmany[i_fusion]+howmany[j_fusion])){
								fusion_OK=1;
							}
							j_fusion--;
						}
						i_fusion++;
					}
				  
					if(fusion_OK==1){
						i1ref=i_fusion-1;		
						i2ref=j_fusion+1;		
					}else{
						kk=k2-1;
					} 
				}else{	
				  i1ref=1;		
				  i2ref=kk;		
				}
			  }else{
				for (igr1=1;igr1<kk;igr1++){
					for (igr2=(igr1+1);igr2<=kk;igr2++){
						if(centroid_k[igr1]!="" && centroid_k[igr2]!=""){
							D1=0.0;	
							RF_tree2tree(D1,centroid_k[igr1],centroid_k[igr2]);	
							if(D1<Dref) //if(D1.lt.Dref) then
							{
								i1ref=igr1;		
								i2ref=igr2;		
								Dref=D1;			
							}	//endif
						}
					  }
				  }
			  }
				
              //Group "i2ref" disappears
			  for (i=1;i<=n;i++)		
			  {
				  if(list[i]==i2ref){
					list[i]=i1ref;	
				  }
				  if(list[i]==i2ref) list[i]=list[i]-1;			
			  }		//70 continue
			  
			  howmany[i1ref]=howmany[i1ref]+howmany[i2ref];		
			  
			  for (k=(i2ref+1);k<=kk;k++)			
			  {
				  howmany[k-1]=howmany[k];		
			  } 
		  }	//end for each k
		//--------------------------------------------------------------------
        // Affichage des organisation des groupes pour chaque nran (random start)
        //--------------------------------------------------------------------
		
		nbInit =0;
		nbFin =0;

		double diff = 0.0;
		
		for(int i=0;i<tabIndices.size();i++){
			nbFin+=tabIndices.at(i);
			for (int j=nbInit; j<nbFin; j++){
				Sref[j]=i+1;
			}
			nbInit+=tabIndices.at(i);
		}
		
		if(intParam==1 || intParam==8 || intParam==9 || intParam==10 || intParam==11 || intParam==12 || intParam==13 || intParam==15 || intParam==20){
			for (k=k1;k>=k2;k--)	
			{
				if (CHr[k]>=CHr_max) {
					CHr_group=k;
					CHr_max=CHr[k];
											
					//Pour évider les clusters vides
					realk=0;
					unique=0;
					CHk = 0;
					
					for(int i=1; i<=k; i++){
						unique=0;
						for(int j=1; j<=n; j++){
							if(listr[CHr_group][j]==i && unique==0){
								CHk++;
								CH_conversion[i] = CHk;
								unique=1;
								realk++;
							}
						}
					}
					
					for(int iz=1; iz<=N; iz++){
						listr[CHr_group][iz]=CH_conversion[listr[CHr_group][iz]];
						Strouve[iz-1]=CH_conversion[listr[CHr_group][iz]];
					}
					CHr_group=realk;	
					
				}	
			}
		}else if(intParam==6){
			for (k=k1;k>=k2;k--)	
			{
				if (CHr[k]>=CHr_max) {
					CHr_group=k;
					CHr_max=CHr[k];
					
					//Pour évider les clusters vides
					realk=0;
					unique=0;
					CHk = 0;
					
					for(int i=1; i<=k; i++){
						unique=0;
						for(int j=1; j<=n; j++){
							if(listr[CHr_group][j]==i && unique==0){
								CHk++;
								CH_conversion[i] = CHk;
								unique=1;
								realk++;
							}
						}
					}
					
					for(int iz=1; iz<=N; iz++){
						listr[CHr_group][iz]=CH_conversion[listr[CHr_group][iz]];
						Strouve[iz-1]=CH_conversion[listr[CHr_group][iz]];
					}
					CHr_group=realk;					
				}	
			}
		}else if(intParam==2){
			for (k=k2;k<=k1;k++)
				{
					if (BHr[k]<BHr_min) {
						BHr_group=k;
						BHr_min=BHr[k];
						
						for(int iz=1; iz<=N; iz++){
							Strouve[iz-1]=listr[BHr_group][iz];
						}
					}	
				}
		}else if(intParam==3 || intParam==16 || intParam==17 || intParam==18 || intParam==19){
			for (k=k1;k>=k2;k--){
				if (Silr[k]>=Silr_max){
					Silr_group=k;
					Silr_max=Silr[k];
											
					//Pour évider les clusters vides
					realk=0;
					unique=0;
					Silhouettek = 0;
					
					for(int i=1; i<=k; i++){
						unique=0;
						for(int j=1; j<=n; j++){
							if(listr[Silr_group][j]==i && unique==0){
								Silhouettek++;
								CH_conversion[i] = Silhouettek;
								unique=1;
								realk++;
							}
						}
					}
					
					for(int iz=1; iz<=N; iz++){
						listr[Silr_group][iz]=CH_conversion[listr[Silr_group][iz]];
						Strouve[iz-1]=CH_conversion[listr[Silr_group][iz]];
					}
					Silr_group=realk;				
				}	
			}
		}else if(intParam==4){
			for (k=k2;k<=k1;k++){
				diff=fabs(LogSSr[k]-LogSSr[k+1]);
		
				if (diff<LogSS_min) {
					LogSS_group=k;
					LogSS_min=diff;
					
					for(int iz=1; iz<=N; iz++){
						Strouve[iz-1]=listr[LogSS_group][iz];
					}
				}	
			}
		}else if(intParam==5 || intParam==7){
			for (k=k1;k>=k2;k--){
				if (Wr[k]>=W_max){
					//pour connaitre le nombre de partition adéquate.
					W_group=k;
					W_max=Wr[k];
					
					realk=0;
					unique=0;
					wk = 0;
					
					for(int i=1; i<=k; i++){
						unique=0;
						for(int j=1; j<=n; j++){
							if(listr[W_group][j]==i && unique==0){
								wk++;
								W_conversion[i] = wk;
								unique=1;
								realk++;
							}
						}
					}
					
					for(int iz=1; iz<=N; iz++){
						listr[W_group][iz]=W_conversion[listr[W_group][iz]];
						Strouve[iz-1]=W_conversion[listr[W_group][iz]];
					}
					W_group=realk;
					
				}	
				
			}
		}else if(intParam==21){
			for (k=k1;k>=k2;k--)	
			{
				if (Gapr[k]>=Gapr_max) {
					Gapr_group=k;
					Gapr_max=Gapr[k];
			
					for(int iz=1; iz<=N; iz++){
						Strouve[iz-1]=listr[Gapr_group][iz];
					}
					
				}	
			}
		}	
		
				
	}  //fin random start
// Print results
	  
		
	switch (intParam){
		case 1:
		{		
			strcpy(criteria, "K-medoid -- CH");
			outStat(Strouve,Sref,criteria,N,N_especes,percent,K_real,CHr_group,CHr_max,listr);
		}break;

		case 3:
		{		
			strcpy(criteria, "K-medoid -- SH");
			outStat(Strouve,Sref,criteria,N,N_especes,percent,K_real,Silr_group,Silr_max,listr);
		}break;
		
		case 5:
		{
			strcpy(criteria, "W");
			int interval = 0;
			Wr_ln[1] = 0.0;
				
			for (k=k2;k<k1;k++){
				Wr_ln[1] += k*abs(Wr[k]-Wr[k+1]);
				interval++;
			}
			
			Wr_ln[1] = (1.8*Wr_ln[1])/interval;
			double terme = Wr_ln[1];
			
			W_max = -10000000000;
			Wr_ln[1] = Wr_ln[1]-abs(Wr[1]-Wr[2]);
			
				
			for (k=k2+1;k<k1;k++){
				Wr_ln[k] = ((k-1)*(abs(Wr[k-1]-Wr[k])-abs(Wr[k]-Wr[k+1])));		
			}
			
 			for (k=k1;k>=k2;k--){
				if (Wr_ln[k]>=W_max){
					W_group=k;
					W_max=Wr_ln[k];
					for(int iz=1; iz<=N; iz++){
						listr[W_group][iz]=W_conversion[listr[W_group][iz]];
						Strouve[iz-1]=W_conversion[listr[W_group][iz]];
					}
				}	
			} 
			outStat(Strouve,Sref,criteria,N,N_especes,percent,K_real,W_group,W_max,listr);
		}break;
		
		case 6:
		{
			strcpy(criteria, "CH -- OLD VERSION");
			outStat(Strouve,Sref,criteria,N,N_especes,percent,K_real,CHr_group,CHr_max,listr);			
		}break;
		
		case 7:
		{
			strcpy(criteria, "W(RF^2)");
			int interval = 0;
			Wr_ln[1] = 0.0;
				
			for (k=k2;k<k1;k++){
				Wr_ln[1] += k*abs(Wr[k]-Wr[k+1]);
				interval++;
			}
			
			Wr_ln[1] = (1.8*Wr_ln[1])/interval;
			double terme = Wr_ln[1];
			
			W_max = -10000000000;
			Wr_ln[1] = Wr_ln[1]-abs(Wr[1]-Wr[2]);
			
				
			for (k=k2+1;k<k1;k++){
				Wr_ln[k] = ((k-1)*(abs(Wr[k-1]-Wr[k])-abs(Wr[k]-Wr[k+1])));		
			}
			
			for (k=k1;k>=k2;k--){
				if (Wr_ln[k]>=W_max){
					W_group=k;
					W_max=Wr_ln[k];
					for(int iz=1; iz<=N; iz++){
						listr[W_group][iz]=W_conversion[listr[W_group][iz]];
						Strouve[iz-1]=W_conversion[listr[W_group][iz]];
					}
				}	
			}
			outStat(Strouve,Sref,criteria,N,N_especes,percent,K_real,W_group,W_max,listr);
			
		}break;
		
		case 8:
		{		
			strcpy(criteria, "Consense");
			outStat(Strouve,Sref,criteria,N,N_especes,percent,K_real,CHr_group,CHr_max,listr);
		}break;
		
		case 9:
		{
			strcpy(criteria, "Lower bound -- CH");
			outStat(Strouve,Sref,criteria,N,N_especes,percent,K_real,CHr_group,CHr_max,listr);			
		}break;
		
		case 10:
		{
			strcpy(criteria, "Average bound -- CH");
			outStat(Strouve,Sref,criteria,N,N_especes,percent,K_real,CHr_group,CHr_max,listr);			
		}break;
		
		case 11:
		{
			strcpy(criteria, "Upper bound -- CH");
			outStat(Strouve,Sref,criteria,N,N_especes,percent,K_real,CHr_group,CHr_max,listr);	
		}break;
		
		case 12:
		{
			strcpy(criteria, "Euclidean -- CH");
			outStat(Strouve,Sref,criteria,N,N_especes,percent,K_real,CHr_group,CHr_max,listr);			
		}break;
		
		case 13:
		{
			strcpy(criteria, "Euclidean(RF^2)");
			outStat(Strouve,Sref,criteria,N,N_especes,percent,K_real,CHr_group,CHr_max,listr);			
		}break;
		
		case 14:
		{
			strcpy(criteria, "CH (Borne avg) -- VERSION 2");
			outStat(Strouve,Sref,criteria,N,N_especes,percent,K_real,CHr_group,CHr_max,listr);			
		}break;
		
		case 15:
		{		
			strcpy(criteria, "Consense(RF^2)");
			outStat(Strouve,Sref,criteria,N,N_especes,percent,K_real,CHr_group,CHr_max,listr);
		}break;
		
		case 16:
		{
			strcpy(criteria, "Lower bound -- SH");
			outStat(Strouve,Sref,criteria,N,N_especes,percent,K_real,Silr_group,Silr_max,listr);			
		}break;
		
		case 17:
		{
			strcpy(criteria, "Average bound -- SH");
			outStat(Strouve,Sref,criteria,N,N_especes,percent,K_real,Silr_group,Silr_max,listr);			
		}break;
		
		case 18:
		{
			strcpy(criteria, "Upper bound -- SH");
			outStat(Strouve,Sref,criteria,N,N_especes,percent,K_real,Silr_group,Silr_max,listr);
		}break;
		
		case 19:
		{
			strcpy(criteria, "Euclidean -- SH");
			outStat(Strouve,Sref,criteria,N,N_especes,percent,K_real,Silr_group,Silr_max,listr);		
		}break;
		
		case 20:
		{
			strcpy(criteria, "EuCH");
			outStat(Strouve,Sref,criteria,N,N_especes,percent,K_real,CHr_group,CHr_max,listr);	
		}break;
		
		case 21:
		{
			strcpy(criteria, "Gap statistic");
			outStat(Strouve,Sref,criteria,N,N_especes,percent,K_real,Gapr_group,Gapr_max,listr);	
		}break;
		
	}
	
	// End timer
    tend2=time(NULL);				// get the current calendar time

	// Compute execution time
    texec2=difftime(tend2,tbegin2);	// tend-tbegin (result in second)
	fprintf (Output4,"%.3f;\n",texec2);
	
	//}// fin de random start
	
	// Print results
	
	
	// *********************Close output files ***************************
	fclose(Output4);
	//*********************** Remove matrix ******************************
	  
	for (i=0;i<=kmax;i++)
	{
		delete [] sx[i];
		delete [] sx2[i];
		delete [] xbar[i];
		delete [] var[i];
		delete [] listr[i];
		delete [] howmanyr[i];
	}
	
	delete [] sx;
	delete [] sx2;
	delete [] xbar;
	delete [] var;
	delete [] listr;
	delete [] howmanyr;
	
	delete [] Dvec;
	
	delete [] CHr;
	delete [] BHr;
	delete [] Gapr;
	delete [] LogSSr;
	delete [] Silr;
	delete [] Wr;
	delete [] Wr_ln;
	delete [] diff_W;
	delete [] V_W;
	
	delete [] SSEr;

	delete [] vect;
	delete [] mean;
	delete [] weight;

	delete [] list;
	delete [] no;
	delete [] iordre;

	delete [] howmany;
	delete [] nobest;
	delete [] nnitr;

	delete [] ishort;

	delete [] nameb;
	delete [] nk;
	delete [] distances_RF_norm;
	
	for (int i=0;i<n;i++)
	{
		delete [] mat_adjacence[i];
		delete [] tree_cluster_leaves[i];
	}
	delete [] mat_adjacence;
	delete [] tree_cluster_leaves;
	delete [] visited;

		
	return 0;
}


 //      end
//************************End of Main

//******************************************************************************
//**********************************FUNCTIONS***********************************
//******************************************************************************

//////////////////////
//***********CheckCentr
//////////////////////

// Modification Centroids: this whole subroutine

//stat output
void outStat(int Strouve[],int Sref[],char *criteria,int N,char *N_especes,char *percent,const char *K_real,int group/* ,double RI,double ARI */,double score,int **listr)
{

	//Compute Rand index between Strouve and Sref
	double RI = f_RI(Strouve,Sref,N);

	//Compute Rand index Adjusted between Strouve and Sref
	double ARI = f_ARI(Strouve,Sref,K_real,group,N);
		
	fprintf (Output4,"%s;",criteria);
	fprintf (Output4,"%i;",N);
	fprintf (Output4,"%s;",N_especes);
	fprintf (Output4,"%s;",percent);
	fprintf (Output4,"%s;",K_real);
	fprintf (Output4,"%d;",group);
	int diff = atoi(K_real)-group;
	diff = fabs(diff);
	const int const_K_real = atoi(K_real);
	const int const_group = group;
	
	double max_k = max(const_K_real,const_group);
	double diff_norm = (diff*1.0)/(max_k*1.0);
	fprintf (Output4,"%i;",diff);
	fprintf (Output4,"%.3f;",diff_norm);
	if(atoi(K_real) == group){
		fprintf (Output4,"%i;",1);
	}else{
		fprintf (Output4,"%i;",0);
	}
	fprintf (Output4,"%.3f;",RI);
	fprintf (Output4,"%.3f;",ARI);	
	fprintf (Output4,"%.3f;",score);
	
	fprintf (Output4,"part(");
	for (int p=1; p<=N; p++){
		if(p==N){
			fprintf (Output4,"%i%s",Strouve[p-1]," ");
		}else{
			fprintf (Output4,"%i%s",Strouve[p-1]," <> ");
		}
		
	}
	fprintf (Output4,");");
	
		fprintf (Output4,"part(");
	for (int p=1; p<=N; p++){
		if(p==N){
			fprintf (Output4,"%i%s",Sref[p-1]," ");
		}else{
			fprintf (Output4,"%i%s",Sref[p-1]," <> ");
		}
		
	}
	fprintf (Output4,");");
}

//compute rand index
double f_RI(int Strouve[],int Sref[],int N)
{
	double comb = 1.0;

	for (int i=N; i>=(N-2+1); i--)
	{
		comb*=i;
	}

	comb/=2;
	 
	double a=0.0;
	double b=0.0; 
							
	for (int i=0; i<N-1; i++){
		for (int j=i+1; j<N; j++){
			if(Sref[i]!=Sref[j]){
				if(Strouve[i]!=Strouve[j]){
					b++;
				}
			}else{
				if(Strouve[i]==Strouve[j]){
					a++;
				}
			}
		}
	}						
	double RI = (a+b)/comb;
	
	return RI;
}

//Convert the partition found by the same number cluster that the partition ref
void conv2sameRef(int *Strouve,int *Sref, int n)
{
	int k = 0;

	//To know the number of cluster
	for(int i=0; i<n; i++){
		if(Strouve[i]>k){
			k=Strouve[i];
		}
	}
	
	int *nk_trouve = new int [k];
	int *nk_ref = new int [k];
	
	for(int i=0;i<k; i++){
		nk_trouve[i]=0;
		nk_ref[i]=0;
	}
	
	//Compter le nombre d'éléments par cluster
	//For Strouve
	for(int i=0;i<k; i++){
		nk_trouve[Strouve[i]]++;
	}

	//For Sref
	for(int i=0;i<k; i++){
		nk_ref[Sref[i]]++;
	}
	
	delete []nk_trouve;
	delete []nk_ref;
}

//compute adjusted rand index
double f_ARI(int Strouve[],int Sref[],const char *K_real,int group,int N)
{
	int kReal = atoi(K_real);
	int tabCongruence [kReal+1][group+1];
	int sumLigne [kReal+1];
	int sumColonne[group+1];
	
	//initialisation du tableau sumLigne à 0
	for(int i=0;i<=kReal;i++){
		sumLigne[i]=0;
	}
	
	//initialisation du tableau sumColonne à 0
	for(int i=0;i<=group;i++){
		sumColonne[i]=0;
	}
	
	//initialisation du tableau des congruence à 0
	for(int i=0;i<=kReal;i++){
		for(int j=0;j<=group;j++){
			tabCongruence[i][j]=0;
		}
	}
	
	
	for(int i=0;i<N;i++){
		tabCongruence[Sref[i]][Strouve[i]]++;
	}
	
	for(int i=1;i<=kReal;i++){
		for(int j=1;j<=group;j++){
			sumLigne[i]+=tabCongruence[i][j];
			sumColonne[j]+=tabCongruence[i][j];
		}
	}
	
	double a=0.0;
	double b=0.0; 
	double c=0.0;
	double d=0.0; 
	
	double comb = 1.0;

	for (int i=N; i>=(N-2+1); i--)
	{
		comb*=i;
	}

	comb/=2;
							
	for (int i=0; i<N-1; i++){
		for (int j=i+1; j<N; j++){
			if(Sref[i]!=Sref[j]){
				if(Strouve[i]!=Strouve[j]){
					d++;
				}else{
					c++;
				}
			}else{
				if(Strouve[i]==Strouve[j]){
					a++;
				}else{
					b++;
				}
			}
		}
	}
	
	double ARI = 0.0;

	if(a*2.0==((b+a)+(c+a))){
		ARI = 1.0;
	}else{
		ARI = a - ((b+a)*(c+a))/comb;
		ARI=ARI/((((b+a)+(c+a))/2.0)-(((b+a)*(c+a))/comb));
	}

	return ARI;
}

void CheckCentr(int &n,int &nmax,int &p,int &pmax,int &k1,int &kmax,double** mat,double** xbar,int* ishort,double** sx,int &idebug)
{

	double temp=0;		//      Real*8 mat(nmax,pmax),sx(kmax,pmax),xbar(kmax,pmax),temp
	int j=0, i=0, kk=0;
//      Integer ishort(pmax)
// Check that the k1 centroids provided are within the range of data values.
// 'sx' is used to store the minimum and maximum values of the variables.
	for (j=1;j<=p;j++){
		sx[1][j]=mat[1-1][j-1];	
		sx[2][j]=mat[1-1][j-1];	
	}						//enddo


	for (i=2;i<=n;i++){
		for (j=1;j<=p;j++){
		    temp=mat[i-1][ishort[j]-1];				//temp=mat(i,ishort(j))
			if(temp<sx[1][j]) sx[1][j]=temp;		//if(temp.lt.sx(1,j)) sx(1,j)=temp
		    if(temp<sx[2][j]) sx[2][j]=temp;		//if(temp.gt.sx(2,j)) sx(2,j)=temp
		}										//enddo
	}											//enddo

	for (kk=1;kk<=k1;kk++){
		for (j=1;j<=p;j++){
			if((xbar[kk][j]<sx[1][j])||(xbar[kk][j]>sx[2][j]))		//if((xbar(kk,j).lt.sx(1,j)).or.(xbar(kk,j).gt.sx(2,j))) then
				{
					printf("Centroids not within the range of the (transformed?) data.");	//stop 'Centroids not within the range of the (transformed?) data.'
					exit(1);
				}				//endif
			}					//enddo
		}						//enddo
    return;
      
}//end
//************************End of CheckCentr

//////////////////////
//***********Assign
//////////////////////

void Assign(int &iran,int &n,int &nmax,int &k1,int &kmax,int* list,int* howmany,int* no,int &idebug,int &iassign,int &iseed, int random_number)
{
      int k=0, i=0, ii=0, itemp=0, kk=0, how=0, isum=0;	
      char namea[255];
      double turn=0;
	  
        // Assign objects to groups.
        // On output, vector 'list' gives group assignments for the objects
        // whereas vector 'howmany' gives the number of objects in each group

   if ((iassign==1) || (iassign==2))   		
   {
	   how=n/k1;
	   for (k=1;k<=(k1-1);k++) {howmany[k]=how;}				   
	   howmany[k1]=n-(k1-1)*how;	
	   ii=0;			

	   for (k=1;k<=k1;k++)		
	   {
		   for (kk=1;kk<=howmany[k];kk++)	
		   {
			   ii++;			
			   list[ii]=k;
		   }
	   }						


	   if(iassign==1) return;	
           // Assign objects at random to the groups
	   if(iran==1)			
	   {
		  for (i=1;i<=(random_number+100);i++)  turn=rand()/(rand() % 32767);
	   }							//end if
           Permute(iseed,n,nmax,list);	
           return;
   }
   else if (iassign==3)
   {
// Read file of group assignments.
// First line: how many objects in each group?
// Then, read members of each group on one line (list of object numbers).
   
	   printf ("Name of file of group assignments?");		//60 write(*,*) 'Name of file of group assignments?'
	   scanf ("%s",namea);		//read(*,*) namea

	   FILE *Input3;
	   if ((Input3 = fopen(namea,"r"))==0) { printf("\n %s :Open Failed....",namea); exit(1); }   	
          
           printf ("File of group assignments: %s\n",namea);
		  		 
	  for (k=1;k<=k1;k++)				
	  {
		  fscanf(Input3,"%d",&howmany[k]);
	  }

      isum=0;							
	  for (k=1;k<=k1;k++)				
	  {
		  isum=isum+howmany[k];			
	  }

      if(isum!=n)					
	  {
		  printf("Objects assigned to groups do not sum to n.");
		  exit(1);
	  }
      	  
	  for (i=1;i<=n;i++) {
		  list[i]=-1;				
	  }

	  for (k=1;k<=k1;k++)			
	  {
		  for (i=1;i<=howmany[k];i++)				
		  {
			  fscanf(Input3, "%d", no[i]);
		  }
     	  for (i=1;i<=howmany[k];i++)		
		  {									
			  list[no[i]]=k;					
		  }
	  }

	  for (i=1;i<=n;i++)	
	  {	
		  if(list[i]==-1)	
		  {
			  printf("Overlapping assignments to groups.");
			  exit(1);
		  }									
	  }
      fclose(Input3);			
      return;
   }
   else 
   {
          printf("Wrong perameter <iassign> in function <Assign>.");
          exit(1);
   }									

} // end
//************************End of Assign


//////////////////////
//***********Centroids
//////////////////////
void Centroids(int &n,int &nmax,int &p,int &pmax,int &kmax,double** mat,double** sx,double** xbar,double* mean,int* list,int* howmany,int &kk,int* ishort,int &idebug)
{  
	int k=0, j=0, i=0, itemp=0;

// Compute centroids for all groups; write to matrix 'xbar'.
// Vector 'list' contains group assignments for the objects.
// Vector 'howmany' contains the number of objects in each group.

	
	//initialisation de sx[][] par 0.0
	for (k=1;k<=kk;k++)		
	{
		for (j=1;j<=p;j++)		
		{
			sx[k][j]=0.0;		
		}
	}

	for (i=1;i<=n;i++)		
	{
		for (j=1;j<=p;j++)		
		{
			itemp=list[i];		
			sx[itemp][j]=sx[itemp][j]+mat[i-1][ishort[j]-1];		//5 sx(itemp,j)=sx(itemp,j)+mat(i,ishort(j))
		}
	}

	for (k=1;k<=kk;k++)		
	{
		itemp=howmany[k];	
		if(itemp==0)			//if(itemp.eq.0) then
		{
			for (j=1;j<=p;j++)		
			{
				xbar[k][j]=mean[j];		//6    xbar(k,j)=mean(j)
			}
		}
        else	
		{
			for (j=1;j<=p;j++)		
			{
				xbar[k][j]=sx[k][j]/static_cast<float>(itemp);		//7    xbar(k,j)=sx(k,j)/dfloat(itemp)
			}
		}			
	}
	

	return;
}			//  end
//************************End of Centroids


//fonctions pour la variante consense

double DistanceCH_consense(int &n,int &kk,vector<string> &centroid_k,string centroid_C_min,double FO_new){	
	double SSB=0.0;
	double SSW=FO_new;
	double distance_total = 0.0;
	int k_cluster = kk;
	double RF = 0.0;

	//compute SSB
	for(int k=1;k<=kk; k++){
		if(centroid_k[k]!="" && centroid_C_min!=""){	
			ofstream myTrees;
			myTrees.open ("myTrees");
			myTrees << centroid_k[k]<<"\n";
			myTrees << centroid_C_min<<"\n";
			myTrees.close();
				
			RF=NouvelleCallRF(16);
			//RF_tree2tree(RF,centroid_k[k],centroid_C_min);
			SSB+=RF;
		}
	}
	
 	if((fabs(SSW)>0.000001) && (k_cluster>1)){
		distance_total=(SSB/SSW)*((n-k_cluster)/(k_cluster-1.0));
	}
	else if(fabs(SSW)<=0.000001 && (k_cluster>1)){
		distance_total=10000000.0*SSB*((n-k_cluster)/(k_cluster-1.0));
	}

	return distance_total;	
    
}//  end ************************End of Distances DistanceCH_consense


void Centroids_consensus(int &n,double** mat,double** xbar,int* list,int &kk, vector<string> &centroid_k, vector<string> monTableau)
{  

	int cluster_k=0;
	char *filename = new char [100];
	string str_Ci = "";	
	//vider les vecteurs des infos des centroids 
	centroid_k.clear();
	
	//initialisation des vecteurs des centroids
	for(int k=0; k<=kk; k++){
		centroid_k.push_back("");
	}
	
	//Mettre tous les arbres de la meme repartition dans des mêmes fichier pour inférer leur consensus majoritaire
	for (int i=1;i<=n;i++){
		cluster_k=list[i];
		
		sprintf(filename,"%s%d","outtree",cluster_k);
		
		/* Ouvertute du fichier qui contiendra tous les arbres de chaque partition */
		ofstream myfile;
		myfile.open (filename, ios::out | ios::app);
		myfile <<monTableau[i-1];
		myfile << "\n";
		myfile.close();
	}

	//Inférer leur consensus majoritaire de chaque partition	
	for (int i=1;i<=kk;i++){
		sprintf(filename,"%s%d","outtree",i);
		ifstream fichierk(filename, ios::in);  // on ouvre en lecture

		if(fichierk){		
			// Ouvertute du fichier contenant les parametres du programme Consense
			ofstream myfile;
			myfile.open ("parametresConsense");
			myfile << filename;
			myfile << "\n";
			myfile << "C\n";
			myfile << "C\n";
			myfile << "Y\n";
			myfile.close();
		 
			//appel du logiciel consense
			system("./consense <parametresConsense >outTmp");
			
			//Format consensus tree Ci for only one line
			system("cat outtree|tr \"\n\" \" \">tmp.txt");
			system("cp tmp.txt outtree");

			//Recuperer le string du consensus des arbres du cluster i par la variable Ci
			ifstream fileCi;
			fileCi.open("outtree");
			
			getline(fileCi,str_Ci);
			centroid_k[i]=(str_Ci);
			fileCi.close(); 
			system("rm outfile outTmp outtree tmp.txt parametresConsense");
			fichierk.close();
		}
		
	}
	system("rm outtree*");
	return;
}			//  end
//************************End of Centroids_consensus

double FO_consense(int &n,int &kk,vector<string> &centroid_k, vector<string> monTableau,double** mat,double* Dvec,int* list,int* howmany,double &SSE,vector<int> &centroid_k_pos,char *N_especes)
{
	double *clusterK_same = new double [kk+1];
	int *nk_CH = new int [kk+1];
	double RF = 0.0;
	double Dref=0,D1=0;		//Real*8 Dref,D1,SSE,Dvec(kk),weight(pmax)
	int	kref=0;		//Integer list(nmax),howmany(kk),kref
	//Integer ishort(pmax)
	// Compute squared distances to group centroids. Assign objects to nearest one
	SSE=0;		//SSE=0.0

	int new_k = 0;
	int old_k = 0;
	int nb_cluster_dest = 0;   
	int nb_cluster_source = 0;  
	double FO_old = 0.0;
	double FO_new = 0.0;
	double tmp_calc = 0.0;
	double tmp_calc_dest = 0.0;
	double tmp_calc_source = 0.0;
	
	double min_dist = 1000000000;
	int k_source = 0;
	
	double nb_esp_identique = atoi(N_especes);
	char *filename = new char [100];
	string str_Ci = "";
	
	double distances[4];	
	for (int j=0; j<4; j++){
		distances[j]=0.0;
	}
	
	//initialisation des variables
	for(int k=0;k<=kk; k++){
		nk_CH[k]=0;
		clusterK_same[k]=0.0;
	}
	
	//Compter le nombre d'éléments par cluster
	for(int k=1;k<=n; k++){
		nk_CH[list[k]]++;
	}
	
	//compute for each cluster initially, SSW value (intra groupe distance)
	//compute SSW
	
	for (int i=1;i<=n;i++){	
		if(centroid_k[list[i]]!=""){	
			ofstream myTrees;
			myTrees.open ("myTrees");
			myTrees << centroid_k[list[i]]<<"\n";
			myTrees << monTableau[i-1]<<"\n";
			myTrees.close();
				
			RF=NouvelleCallRF(16);
			
			/*main_hgt(monTableau[i-1],centroid_k[list[i]],distances);
			RF = distances[0];*/	
			
			//RF_tree2tree(RF,monTableau[i-1],centroid_k[list[i]]);
			
			clusterK_same[list[i]]+=RF;
		}
	} 

	//compute SSW at the beginning (distance intra groupe)
	for (int k=1;k<=kk;k++){
		if(nk_CH[k]>1){
			FO_old += clusterK_same[k];
		}
	}
	
	for (int i=1;i<=n; i++)			//do 20 i=1,n
	{
		k_source = list[i];
		if(nk_CH[list[i]]>1){
			for (int k=1;k<=kk;k++)			//do 12 k=1,kk
			{				
				//Calcul de la distance RF de chaque point i
				// et assignation du point i au bon cluster
				// Compute a RF distance to the centroid k
				
				//test si le point i n'appartenait pas initiallement à k 
				//ET que le point i n'est pas le centroid de centroid_k_pos[i] ---  && i!=centroid_k_pos[k]
				if(k_source!=k && list[i]!=k){
					tmp_calc_source = 0.0;
					tmp_calc_dest = 0.0;
					nb_cluster_dest = nk_CH[k];
					nb_cluster_source = nk_CH[list[i]];
					
					FO_new = FO_old - clusterK_same[k] - clusterK_same[list[i]];	

					//Mettre tous les arbres de la meme repartition dans des mêmes fichier pour inférer leur consensus majoritaire
					
					//Cluster source (list[i])
					
					//Récupérer tous les arbres du cluster source (list[i]) sans l'arbre i
					
					for (int j=1;j<=n;j++){
						if(list[i]==list[j] && i!=j){
							sprintf(filename,"%s%d","outtree",list[i]);
						
							/* Ouvertute du fichier qui contiendra tous les arbres de chaque partition */
							ofstream myfile;
							myfile.open (filename, ios::out | ios::app);
							myfile <<monTableau[j-1];
							myfile << "\n";
							myfile.close();
						}
					}
					
					//Inférer le nouvel arbre consensus majoritaire du cluster source (list[i])
					sprintf(filename,"%s%d","outtree",list[i]);
					ifstream fichieri(filename, ios::in);  // on ouvre en lecture
					
					if(fichieri){		
						// Ouvertute du fichier contenant les parametres du programme Consense
						ofstream myfile;
						myfile.open ("parametresConsense");
						myfile << filename;
						myfile << "\n";
						myfile << "C\n";
						myfile << "C\n";
						myfile << "Y\n";
						myfile.close();
					 
						//appel du logiciel consense
						system("./consense <parametresConsense >outTmp");
						
						//Format consensus tree Ci for only one line
						system("cat outtree|tr \"\n\" \" \">tmp.txt");
						system("cp tmp.txt outtree");

						//Recuperer le string du consensus des arbres du cluster i par la variable Ci
						ifstream fileCi;
						fileCi.open("outtree");
						
						getline(fileCi,str_Ci);
						centroid_k[list[i]]=str_Ci;
						fileCi.close(); 
						system("rm outfile outTmp outtree tmp.txt parametresConsense");
						fichieri.close();
					}
					
					//Calculer la distance intra-groupe pour le cluster source (list[i]) = clusterK_same[list[i]]			
					for (int j=1;j<n;j++){	
						if(list[i]==list[j] && i!=j){
							if(centroid_k[list[i]]!=""){
								ofstream myTrees;
								myTrees.open ("myTrees");
								myTrees << centroid_k[list[i]]<<"\n";
								myTrees << monTableau[j-1]<<"\n";
								myTrees.close();
									
								RF=NouvelleCallRF(16);
								tmp_calc_source+=RF;
							}
						}
					} 
					
					//mise à jour du nombre d'arbres dans le cluster source (list[i])
					nb_cluster_source = nb_cluster_source - 1;
					
					//Cluster destination (k)
					
					//Récupérer tous les arbres du cluster destination (k) avec l'arbre i
					for (int j=1;j<=n;j++){
						if(k==list[j] || i==j){
							sprintf(filename,"%s%d","outtree",k);
						
							/* Ouvertute du fichier qui contiendra tous les arbres de chaque partition */
							ofstream myfile;
							myfile.open (filename, ios::out | ios::app);
							myfile <<monTableau[j-1];
							myfile << "\n";
							myfile.close();
						}
					}
					
					//Inférer le nouvel arbre consensus majoritaire du cluster destination (k)
					sprintf(filename,"%s%d","outtree",k);
					ifstream fichierk(filename, ios::in);  // on ouvre en lecture

					if(fichierk){	
						// Ouvertute du fichier contenant les parametres du programme Consense
						ofstream myfile;
						myfile.open ("parametresConsense");
						myfile << filename;
						myfile << "\n";
						myfile << "C\n";
						myfile << "C\n";
						myfile << "Y\n";
						myfile.close();
					 
						//appel du logiciel consense
						system("./consense <parametresConsense >outTmp");
						
						//Format consensus tree Ci for only one line
						system("cat outtree|tr \"\n\" \" \">tmp.txt");
						system("cp tmp.txt outtree");

						//Recuperer le string du consensus des arbres du cluster i par la variable Ci
						ifstream fileCi;
						fileCi.open("outtree");
						
						getline(fileCi,str_Ci);
						centroid_k[k]=str_Ci;
						fileCi.close(); 
						system("rm outfile outTmp outtree tmp.txt parametresConsense");
						fichierk.close();
					}

					system("rm outtree* myTrees");
					
					//Calculer la distance intra-groupe pour le cluster destination (k) = clusterK_same[k]
					for (int j=1;j<n;j++){	
						if(k==list[j] || i==j){
							if(centroid_k[k]!=""){
								ofstream myTrees;
								myTrees.open ("myTrees");
								myTrees << centroid_k[k]<<"\n";
								myTrees << monTableau[j-1]<<"\n";
								myTrees.close();
									
								RF=NouvelleCallRF(16);
								tmp_calc_dest+=RF;
							}
						}
					} 
					
					//mise à jour du nombre d'arbres dans le cluster destination (k)
					nb_cluster_dest = nb_cluster_dest + 1;
					
					FO_new = FO_new + tmp_calc_dest + tmp_calc_source;	
					
					if(FO_new<FO_old){
						Dref=FO_new;		
						kref=k;
						
						//A VOIR SI UTILE
						new_k = k;
						old_k = list[i];	
						
						//mise à jour de nk_CH[]
						nk_CH[k] = nb_cluster_dest;
						nk_CH[list[i]] = nb_cluster_source;					
						
						//mise à jour de la distance intra-groupe des deux clusters modifiés 
						clusterK_same[list[i]] = tmp_calc_source;
						clusterK_same[k] = tmp_calc_dest;
						
						//mise à jour la liste de distribution des elements
						list[i] = k;	
												
						//mise à jour de la fonction objective FO_old
						FO_old = FO_new;				

					}
									
				}
			}		
		}
		SSE=SSE+Dref;         //SSE=SSE+Dref		 
		howmany[kref]++;	//howmany(kref)=howmany(kref)+1
		
	}
	
    return FO_new;
    
}//  end FO_consense


//Fonctions pour la variante K-medoid
void Centroids_approx(int &n,double** mat,double** xbar,int* list,int &kk, vector<int> &centroid_k_pos)
{  
	//vider les vecteurs des infos des centroids approx
	centroid_k_pos.clear();
	
	int k=0, j=0, i=0, itemp=0;
	double dist_intra[kk+2][n+1];

	//initialisation des variables
	for(k=0; k<=kk; k++){
		for(i=0; i<=n; i++){
			dist_intra[k][i]=1000000000;
		}
	}
	
	//Calculer pour chque arbre sa distance intra groupe avec tous les autres points de ce meme groupe
	int noCluster = 0;
	for(i=1; i<=n; i++){
		noCluster = list[i];
		dist_intra[noCluster][i] = 0.0;
		for(j=1; j<=n; j++){
			if(list[j]==noCluster){
				dist_intra[noCluster][i]+=mat[i-1][j-1];
			}
		}
	}
	
	//initialisation des vecteurs des centroids
	for(k=0; k<=kk; k++){
		centroid_k_pos.push_back(0);
	}
	
	//chercher l'arbres de chque groupe qui minimise les distance
	//pour etre considere comme arbre centroid de ce groupe
	double min_dist = 1000000000;
	for(k=1; k<=kk; k++){
		min_dist = 1000000000;
		for(j=1; j<=n; j++){
			if(dist_intra[k][j]<min_dist){
				xbar[k][j]=dist_intra[k][j];
				centroid_k_pos[k]=j;
			}
			
		}
		
	}
	
	return;
}			//  end
//************************End of Centroids_approx


double DistanceCH_approx_Kmedoid(int &n,int &kk,double** mat,vector<int> &centroid_k_pos,double FO_new,int indice_arbre_cons,int* list) {
	
	double SSB=0.0;
	double SSW=FO_new;
	double distance_total = 0.0;
	int k_cluster = kk;
	int *nk_CH = new int [kk+1];

	//initialisation des variables
	for(int k=0;k<=kk; k++){
		nk_CH[k]=0;
	}
	
	//Compter le nombre d'éléments par cluster
	for(int k=1;k<=n; k++){
		nk_CH[list[k]]++;
	}
	
	//compute SSB
	for(int k=1;k<=kk; k++){
		SSB+=(nk_CH[k]*mat[centroid_k_pos[k]-1][indice_arbre_cons-1]);
	}
	
 	if((fabs(SSW)>0.000001) && (k_cluster>1)){
		distance_total=(SSB/SSW)*((n-k_cluster)/(k_cluster-1.0));
	}
	else if(fabs(SSW)<=0.000001 && (k_cluster>1)){
		distance_total=10000000.0*SSB*((n-k_cluster)/(k_cluster-1.0));
	}
	
	return distance_total;	
    
}//  end ************************End of Distances DistanceCH_approx_Kmedoid


double FO_approx_Kmedoid(int &n,int &kmax,double** mat,double* Dvec,int* list,int* howmany,double &SSE,int &kk,vector<int> &centroid_k_pos, double ** n_identique)
{
	double *clusterK_same = new double [kk+1];
	int *nk_CH = new int [kk+1];
	int cluster_k=0;
	double RF = 0.0;
	double Dref=0,D1=0;		//Real*8 Dref,D1,SSE,Dvec(kk),weight(pmax)
	int	kref=0;		//Integer list(nmax),howmany(kk),kref
	//Integer ishort(pmax)
	// Compute squared distances to group centroids. Assign objects to nearest one
	SSE=0;		//SSE=0.0

	int new_k = 0;
	int old_k = 0;
	int nb_cluster_dest = 0;   
	int nb_cluster_source = 0;  
	double FO_old = 0.0;
	double FO_new = 0.0;
	double tmp_calc = 0.0;
	double tmp_calc_dest = 0.0;
	double tmp_calc_source = 0.0;
	double dist_intra[kk+1][n+1];
	double dist_intra_tmp[kk+1][n+1];
	
	double min_dist = 1000000000;
	int k_source = 0;

	//initialisation des distances intra groupe
	for(int k=0; k<=kk; k++){
		for(int i=0; i<=n; i++){
			dist_intra[k][i]=1000000000.0;
			dist_intra_tmp[k][i]=1000000000.0;
		}
	}
	
	//initialisation des variables
	for(int k=0;k<=kk; k++){
		nk_CH[k]=0;
		clusterK_same[k]=0.0;
	}
	
	//Compter le nombre d'éléments par cluster
	for(int k=1;k<=n; k++){
		nk_CH[list[k]]++;
	}
	
	//compute for each cluster initially, SSW value (intra groupe distance)
	//compute SSW
	for (int i=1;i<n;i++){
		clusterK_same[list[i]]+=mat[i-1][centroid_k_pos[list[i]]-1];
	} 
	
	for (int k=1;k<=kk;k++){
		if(nk_CH[k]>1){
			FO_old += clusterK_same[k];
		}
	}
	
	//les distances intra-groupe (au départ) de chaque arbre vous chaque groupe
	for(int i=1; i<=n; i++){
		dist_intra_tmp[list[i]][i] = 0.0;
		for(int j=1; j<=n; j++){
			if(list[j]==list[i]){
				dist_intra_tmp[list[i]][i]+=mat[i-1][j-1];
			}
		}
	}
	
	//copier la matrice de distance intra
	for(int k_cluster=0; k_cluster<=kk; k_cluster++){
		for(int i_element=0; i_element<=n; i_element++){
			dist_intra[k_cluster][i_element]=dist_intra_tmp[k_cluster][i_element];
		}
	}
  
	for (int i=1;i<=n; i++)			//do 20 i=1,n
	{
		k_source = list[i];
		if(nk_CH[list[i]]>1){
			for (int k=1;k<=kk;k++)			//do 12 k=1,kk
			{				
				//Calcul de la distance RF de chaque point i
				// et assignation du point i au bon cluster
				// Compute a RF distance to the centroid k

				//test si le point i n'appartenait pas initiallement à k 
				//ET que le point i n'est pas le centroid de centroid_k_pos[i] ---  && i!=centroid_k_pos[k]
				if(k_source!=k && list[i]!=k){
					nb_cluster_dest = nk_CH[k];
					nb_cluster_source = nk_CH[list[i]];
					
					FO_new = FO_old - clusterK_same[k] - clusterK_same[list[i]];	

					
					//compute Function objective
					
					//Pour le cluster de destination	
					//Mettre a jour la distance intra de l'element i avec le cluster de destination (k)
					dist_intra_tmp[k][i] = 0.0;
					for(int j=1; j<=n; j++){
						if(list[j]==k){
							dist_intra_tmp[k][i]+=mat[i-1][j-1];
							dist_intra_tmp[k][j]+=mat[j-1][i-1];
						}
					}
					

					//chercher le nouveau pseudo-cendroid du cluster destination (k)
					//chercher l'arbres de chque groupe qui minimise les distance
					//pour etre considere comme arbre centroid de ce groupe
					
					min_dist = 1000000000;
					for(int j=1; j<=n; j++){
						if(dist_intra_tmp[k][j]<min_dist){
							min_dist = dist_intra_tmp[k][j];
							tmp_calc_dest=dist_intra_tmp[k][j];
							centroid_k_pos[k]=j;
						}						
					}			
					
					nb_cluster_dest +=1;
					FO_new = FO_new + tmp_calc_dest;
					
					//Pour le cluster source
					//Mettre a jour la distance intra de l'element i avec le cluster source (list[i]) a 10000000000.0
					dist_intra_tmp[list[i]][i] = 10000000000.0;
					for(int j=1; j<=n; j++){
						if(list[j]==list[i] && i!=j){
							dist_intra_tmp[list[i]][j]-=mat[j-1][i-1];
						}
					}
					
					
					//chercher le nouveau pseudo-cendroid du cluster source (list[i])
					//chercher l'arbres de chque groupe qui minimise les distance
					//pour etre considere comme arbre centroid de ce groupe
					
					min_dist = 1000000000;
					for(int j=1; j<=n; j++){
						if(dist_intra_tmp[list[i]][j]<min_dist){
							min_dist = dist_intra_tmp[list[i]][j];
							tmp_calc_source=dist_intra_tmp[list[i]][j];
							centroid_k_pos[list[i]]=j;
						}						
					}			
					
					//chercher le nouveau pseudo-cendroid du cluster source (list[i])
					nb_cluster_source -=1;
					FO_new = FO_new + tmp_calc_source;
					
					
					//Dvec[k]=D1;
					if(FO_new<FO_old){
						Dref=FO_new;		
						kref=k;
						
						//A VOIR SI UTILE
						new_k = k;
						old_k = list[i];	
						
						//mise à jour de nk_CH[]
						nk_CH[k] = nb_cluster_dest;
						nk_CH[list[i]] = nb_cluster_source;					
						
						//mise à jour de la distance intra-groupe des deux clusters modifiés 
						clusterK_same[list[i]] = tmp_calc_source;
						clusterK_same[k] = tmp_calc_dest;
						
						//mise à jour la liste de distribution des elements
						list[i] = k;	
												
						//mise à jour de la fonction objective FO_old
						FO_old = FO_new;
						
						//Mettre a jour les distances intra groupes
						for(int k_cluster=0; k_cluster<=kk; k_cluster++){
							for(int i_element=0; i_element<=n; i_element++){
								dist_intra[k_cluster][i_element]=dist_intra_tmp[k_cluster][i_element];
							}
						}						

					}else{
						//mise à jour de la fonction objective FO_old
						FO_new = FO_old;
						
						//Reinitialiser distances intra groupes du tmp
						for(int k_cluster=0; k_cluster<=kk; k_cluster++){
							for(int i_element=0; i_element<=n; i_element++){
								dist_intra_tmp[k_cluster][i_element]=dist_intra[k_cluster][i_element];
							}
						}
					}
									
				}
			}		
		}
		SSE=SSE+Dref;         //SSE=SSE+Dref		 
		howmany[kref]++;	//howmany(kref)=howmany(kref)+1
		
	}
    return FO_new;
    
}//  end FO_approx_Kmedoid


////////////////////
//***********CompSST
////////////////////
void CompSST(int &n,int &nmax,int &p,int &pmax,double** mat,double* weight,int* ishort,double &SST)
{
	double	sx=0,sx2=0,var=0,temp=0,dfln=0;	 //Real*8 mat(nmax,pmax),weight(pmax),sx,sx2,var,temp,dfln,SST
	int j=0, i=0;
      dfln=n;		//dfln=dfloat(n)
      SST=0.0;				//SST=0.0

	  for (j=1;j<=p;j++)		// do 22 j=1,p
	  {
		  sx=0.0;		//sx=0.0
		  sx2=0.0;		//sx2=0.0
		  for (i=1;i<=p;i++)		// do 20 i=1,n
		  {    
			  temp=mat[i-1][j-1];
			  sx=sx+temp;				//sx=sx+temp
			  sx2=sx2+temp*temp;		//20 sx2=sx2+temp*temp
		  }
		  var=sx2-(sx*sx/dfln);			//var=sx2-(sx*sx/dfln) 
		  SST=SST+var*weight[ishort[j]];		//22 SST=SST+var*weight(ishort(j))
	  }
      return;
}				//      end
//************************End of CompSST


//////////////////////
//***********Distances
//////////////////////
      
void Distances(int &n,int &nmax,int &p,int &pmax,int &kmax,double** mat,double** xbar,double* Dvec,int* list,int* howmany,double &SSE,int &kk,double* weight,int* ishort,int &idebug, double &W, double ** n_identique)
{

	  double Dref=0,D1=0;		//Real*8 Dref,D1,SSE,Dvec(kmax),weight(pmax)
      int	kref=0, k=0,i=0;		//Integer list(nmax),howmany(kmax),kref
      //Integer ishort(pmax)
        // Compute squared distances to group centroids. Assign objects to nearest one
      SSE=0;		//SSE=0.0
		

	  for (k=1;k<=kk;k++)		// do 10 k=1,kk
	  {
		  howmany[k]=0;			//10 howmany(k)=0
	  }

     
// 
      for (i=1;i<=n; i++)			//do 20 i=1,n
	  {

		  for (k=1;k<=kk;k++)			//do 12 k=1,kk
		  {
		  
			//Calcul de la distance eucledian de chaque point i
			// et assignation du point i au bon cluster
			// Compute a SQUARED (weighted) Euclidean distance to the centroid k
			Euclidean(i,k,nmax,kmax,p,pmax,mat,xbar,weight,ishort,D1);	
			Dvec[k]=D1;
			if(k==1) 			
			{
			  Dref=D1;		
			  kref=1;		
			}
			else		
			{
			  if(D1<Dref)		
			  {
				  Dref=D1;		
				  kref=k;		
			  }		
			}			
		}			  
		  
		  
		  
		SSE=SSE+Dref;         //SSE=SSE+Dref		 
		howmany[kref]++;	//howmany(kref)=howmany(kref)+1

		/* nbNi[kref]=howmany[kref];	//howmany(kref)=howmany(kref)+1 */

		list[i]=kref;		//  20 list(i)=kref

	  }
				   
      return;
    
}//  end

double DistanceCH_old(int &n,int &kmax,double** mat,int* list,double** Ww)
{
	
	double SSB=0.0;
	double SSW=0.0;
	double dist_all = 0.0;
	
	double distance_total = 0.0;
	int cluster_k=0;
	double *clusterK = new double [kmax+1];
	double *clusterK_same = new double [kmax+1];
	int *nk_CH = new int [kmax+1];
	int k_cluster = 0;
	
	for(int k=1;k<=kmax; k++){
		nk_CH[k]=0;
		clusterK[k]=0.0;
		clusterK_same[k]=0.0;
	}
	
	
	for(int k=1;k<=kmax; k++){
		nk_CH[list[k]]++;
	}
	
	for(int k=1;k<=kmax; k++){
		if(nk_CH[k]!=0){
			k_cluster++;
		}
	}
	
	double RF;
	
	//compute dist_all
	for (int i=1;i<n;i++){
		for (int j=i+1;j<=n;j++){
			RF = mat[i-1][j-1];
			dist_all += RF*Ww[i-1][j-1];
		}
	}
	if(n!=0){
		dist_all = dist_all/n;
	}
		
	//compute SSW
	for (int i=1;i<n;i++){
		cluster_k=list[i];
		for (int j=i+1;j<=n;j++){
			if (list[j]==cluster_k){
				RF = mat[i-1][j-1];
				clusterK_same[cluster_k]+=(RF*Ww[i-1][j-1]);
			}
		}
	}
	
	for(int k=1;k<=kmax; k++){
		if(nk_CH[k]>1){
			clusterK_same[k]*=(1.0/(nk_CH[k]));
			SSW+=clusterK_same[k];
		}
	}
	
	//compute SSB
	SSB = dist_all - SSW;
	
	
 	if((fabs(SSW)>0.000001) && (k_cluster>1)){
		distance_total=(SSB/SSW)*((n-k_cluster)/(k_cluster-1.0));
	}
	else if(fabs(SSW)<=0.000001 && (k_cluster>1)){
		distance_total=10000000.0*SSB*((n-k_cluster)/(k_cluster-1.0));
	}

	delete [] clusterK;
	delete [] clusterK_same;
	delete [] nk_CH;
	
	return distance_total;
      
	
    
}//  end ************************End of Distances DistanceCH_old


void Distances_approx(int &n,int &nmax,int &p,int &pmax,int &kmax,double** mat,double** xbar,double* Dvec,int* list,int* howmany,double &SSE,int &kk,double* weight,int* ishort,int &idebug, double &W,vector <string> monTableau,vector<string> &centroid_k,vector<int> &centroid_k_pos)
{

	  double Dref=0,D1=0;		//Real*8 Dref,D1,SSE,Dvec(kmax),weight(pmax)
      int	kref=0;		//Integer list(nmax),howmany(kmax),kref
      //Integer ishort(pmax)
      // Compute squared distances to group centroids. Assign objects to nearest one
      SSE=0;		//SSE=0.0
 
      for (int i=1;i<=n; i++)			//do 20 i=1,n
	  {

		  for (int k=1;k<=kk;k++)			//do 12 k=1,kk
		  {
		  
			//Calcul de la distance RF de chaque point i
			// et assignation du point i au bon cluster
			// Compute a RF distance to the centroid k

			//test si le centroid est vide
			if(centroid_k[k]!=""){
				RF_tree2tree(D1,monTableau[i-1],centroid_k[k]);	
				
				//compute SSW
				Dvec[k]=D1;
				
				if(k==1) 			
				{
				  Dref=D1;		
				  kref=1;		
				}
				else		
				{
				  if(D1<Dref)		
				  {
					  Dref=D1;		
					  kref=k;		
				  }		
				}	
			}
		}	
		  
		SSE=SSE+Dref;         //SSE=SSE+Dref		 
		howmany[kref]++;	//howmany(kref)=howmany(kref)+1


		list[i]=kref;		//  20 list(i)=kref

	  }

      return;
    
}//  end



double DistanceSilhouette(int &n,int &nmax,int &p,int &pmax,int &kmax,double** mat,double** xbar, int* howmany, int &kk,int* list,double* weight,int* ishort) {
	
    //--Note: we suppose that we have a 
    
      double Dref=0,D1=0;		//Real*8 Dref,D1,SSE,Dvec(kmax),weight(pmax)
      int kref=0;		//Integer list(nmax),howmany(kmax),kref
      double* ai =new double[n+1]; //Distance of element i with all elements in the same cluster than i
      double* bi =new double[n+1]; 
      double* si =new double[n+1];
      double* sk =new double[kk+1];    
      double* dist_cluster=new double[kmax+1]; //--distance of i to other cluster (tmp)
	  int cluster_i=0;
	  double distance_i=0.0;
	  int cluster_j=0;
	  int nk = 0;
	  
	  int *nk_silhouette = new int [kmax+1];
	
	for(int k=1;k<=kmax; k++){
		nk_silhouette[k]=0;
	}
	
	//Compute the number of element for each cluster and
	//stocke them on nk_silhouette[]
	//nk_silhouette[] this array start at 1
	for(int k=1;k<=kmax; k++){
		nk_silhouette[list[k]]++;
	}	  
      for (int k=1; k<=kk;k++){
		sk[k]=0.0;
      }
	  
      for (int i=1;i<=n; i++)		
	  {
          //--Calculate a[i] -> calculate the          
          //--Calculate the distance of i to all other point in its cluster
          cluster_i=list[i];
          distance_i=0.0;
          
		  for (int j=1;j<=n;j++){
              if (list[j]==cluster_i) {
                  distance_i+=mat[i-1][j-1];
              }
          }
          
		  nk = nk_silhouette[cluster_i];
   		   if(nk>0){
			   ai[i]=(distance_i/nk); 
		   }else{
               ai[i]=0.0; //--Alone element in cluster             
           }

          //--Calculate b[i] distance of i to the mean distance of each other cluster
          bi[i]=100000000.0;
          
          for (int k=1;k<=kk;k++){
			dist_cluster[k]=0.0;
		  }
		  
          for (int j=1;j<=n;j++) {             
              cluster_j=list[j];

              if (cluster_j!=cluster_i) {
                  dist_cluster[cluster_j]+=mat[i-1][j-1];
              }
          }
		  
          //--Mean distance
           for (int k=1;k<=kk;k++) {
              if (nk_silhouette[k]>0){
				dist_cluster[k]=dist_cluster[k]/(nk_silhouette[k]);
			  }
          }
		  
		  
        //--b[i] is the minimum
        for (int k=1;k<=kk;k++){
            if (k!=cluster_i && dist_cluster[k]<bi[i]){
				/*if (dist_cluster[k]==0){
					bi[i]=0.001;
				}else{*/	
					bi[i]=dist_cluster[k];
				//}
			}			
        }
		
		/* cout<<"bi["<<i<<"] = "<<bi[i]<<endl;   */
		
		if(ai[i]==bi[i]){
			si[i]=0.0;
		}else if(ai[i]<bi[i]){
			if(bi[i]!=0){
				si[i]=1.0-(ai[i]/bi[i]);
			}
		}else if(ai[i]>bi[i]){
			if(ai[i]!=0){
				si[i]=(bi[i]/ai[i])-1.0;
			}
		}
		
		sk[cluster_i]+=si[i];

      }
	  
      //--Average the sk (note: nk_silhouette can be 0)
      for (int k=1;k<=kk;k++){
		if(nk_silhouette[k]>0){
			sk[k]=sk[k]/nk_silhouette[k];
		}
	  }
      
       double C=0.0;
        for (int k= 1;k <= kk; k++){
			C+=sk[k];
        }
		if(kk!=0){
			C = C/(kk*1.0);
		}

      //--Clean up
	  delete [] dist_cluster;
	  delete [] nk_silhouette;
	  return C;
      

    
}//  end ************************End of Distances

double FO_silhouette(int &n,int &kmax,double** mat,double* Dvec,int* list,int* howmany,double &SSE,int &kk,vector<int> &centroid_k_pos)
{
	   //--Note: we suppose that we have a 
    
      double D1=0;		
      int kref=0;
      double* dist_cluster=new double[kmax+1]; //--distance of i to other cluster (tmp)
	  double min_distance = 10000000000.0;
	  
      for (int i=1;i<n; i++)		
	  {    
		  for (int j=1;j<=n;j++) {
              dist_cluster[list[j]]+=mat[i-1][j-1];
          }
          
		  for (int k=1;k<=kk;k++) {
			if(howmany[k]>0){
				dist_cluster[k] = dist_cluster[k]/howmany[k];
			}
		  }
		  
		  //trouver le min des distance du point i pour l'affecter
		  min_distance = 10000000000.0;
		  
		  for (int k=1;k<=kk;k++) {
			if(min_distance<dist_cluster[k]){
				min_distance = dist_cluster[k];
				list[i] = k;
			}
		  }
		  

      }
        
      //--Clean up
      delete [] dist_cluster;
      
	  return D1;
	
}//  end ************************End of FO_silhouette

double DistanceCH_supertree(int &n,int &kmax,double** mat,int* list,double** Ww,double FO_new,double facteur, double ** tree_cluster_leaves, double ** n_identique) {
	
	double SSB=0.0;
	double SSW=0.0;
	double dist_all = 0.0;
	double distance_total = 0.0;
	int cluster_k=0;
	int *nk_CH = new int [kmax+1];

	int k_cluster = 0;
	
	for(int k=1;k<=kmax; k++){
		nk_CH[k]=0;
	}
	
	for(int k=1;k<=kmax; k++){
		nk_CH[list[k]]++;
	}
	
	for(int k=1;k<=kmax; k++){
		if(nk_CH[k]!=0){
			k_cluster++;
		}
	}
	
	//compute dist_all	
	int Nc = 0;
	double dist_by_el = 0.0;
	for (int i=1;i<=n;i++){
		 Nc = 0;
		 dist_by_el = 0.0;
		for (int j=1;j<=n;j++){
			if(n_identique[i-1][j-1]>3){
				dist_by_el += mat[i-1][j-1];
				Nc++;
			}
		}
		if(Nc>=2){
			if(Nc!=n){
				dist_all+=((dist_by_el/(Nc-1))*(n-(Nc-1.0)));
			}else{
				dist_all+=dist_by_el;
			}
		}else{
			dist_all+=n-1;
		}
	} 
	if(n!=0){
		dist_all = dist_all/(2.0*n);
	}
	
	//compute SSW
	SSW = FO_new;
	
	//compute SSB
	SSB = dist_all-SSW;
 	if((fabs(SSW)>0.000001) && (k_cluster>1)){
		distance_total=(SSB/SSW)*((n-k_cluster)/(k_cluster-1.0));
	}
	else if(fabs(SSW)<=0.000001  && (k_cluster>1)){
		distance_total=10000000.0*SSB*((n-k_cluster)/(k_cluster-1.0));
	} 

	delete [] nk_CH;
	return distance_total;
    
}//  end ************************End of Distances DistanceCH_supertree


double FO_super_tree(int &n,int &kmax,double** mat,double* Dvec,int* list,int* howmany,double &SSE,int &kk,vector <string> monTableau, double ** n_identique, double ** tree_cluster_leaves)
{
	double *clusterK_same = new double [kmax+1];
	int *nk_CH = new int [kmax+1];
	int cluster_k=0;
	double RF = 0.0;
	double Dref=0;		
	int	kref=0;		//Integer list(nmax),howmany(kmax),kref
	// Integer ishort(pmax)
	// Compute squared distances to group centroids. Assign objects to nearest one
	SSE=0;		//SSE=0.0

	int k_source = 0;
	int new_k = 0;
	int old_k = 0;
	int nb_cluster_dest = 0;   
	int nb_cluster_source = 0;  
	double FO_old = 0.0;
	double FO_new = 0.0;
	double tmp_calc_dest = 0.0;
	double tmp_calc_source = 0.0;

	int * p = new int [n+1]; //nb de couples d'arbres n'ayant pas de feuilles communes ou inférieur à 3 d'un même cluster
	int **mat_adjacence = new int *[n];
	//Boucle interne de la recherche du pivot de matrice d'adjecence (celui ne contient pas 2)
	int aret_de_recherche = 0;
	int depart_graphe_i = 0;
	int depart_graphe_j = 0;
	double nb_element_dest_i = 0;
	double RF_element_dest_i = 0;
	int nb_el_dest=0;
	int trouve = 0;	
	int k_cluster = 0;
	double tmp_dist_supp = 0;
	
	int * visited = new int[n];//indique si les elements ont ete visites
	int nb_visited = 0; //compte le nombre d'element visite 
	
	for(int i=0;i<n;i++){
		mat_adjacence[i]=new int [n];
	}
	
	//initialisation des variables
	for(int k=1;k<=kmax; k++){
		nk_CH[k]=0;
		p[k] = 0;
		clusterK_same[k]=0.0;
	}
	
	//Initialisation du tableau tree_cluster_leaves
	for(int i=0; i<n; i++){
		tree_cluster_leaves[i][0]=i+1.0;
		tree_cluster_leaves[i][1]=list[i+1];
		tree_cluster_leaves[i][2]=0.0;
		tree_cluster_leaves[i][3]=0.0;
	}
	
	//Compter le nombre d'éléments par cluster
	for(int k=1;k<=kmax; k++){
		nk_CH[list[k]]++;
	}
	
	//mise à jour du nombre d'element par cluster
	for(int k=1;k<=kmax; k++){
		if(nk_CH[k]!=0){
			k_cluster++;
		}
		howmany[k]=nk_CH[k];
	}
	
	//compute for each cluster initially, SSW value (intra groupe distance)
	//compute SSW
	for (int i=1;i<=n;i++){
		cluster_k=list[i];
		for (int j=1;j<=n;j++){
			if (list[j]==cluster_k){
				if(n_identique[i-1][j-1]>3){
					RF = mat[i-1][j-1];
					clusterK_same[cluster_k]+=RF;
					tree_cluster_leaves[i-1][2]+=1;
					tree_cluster_leaves[i-1][3]+=RF;
				}else{
					//compte le nombre de paire d'arbres i et j ayant un nombre de feuilles chevauchants inf à 3
					p[cluster_k]++;
				}
			}
		}
	} 
	
	for (int k=1;k<=kk;k++){
		//tmp_dist_supp = 0;
		for(int nb_p_arbre = 0; nb_p_arbre<n; nb_p_arbre++){
			if(tree_cluster_leaves[nb_p_arbre][1]==double(k)){
				if(tree_cluster_leaves[nb_p_arbre][2]>=2){
					if(nk_CH[k]!=tree_cluster_leaves[nb_p_arbre][2]){
						tree_cluster_leaves[nb_p_arbre][3]+=((tree_cluster_leaves[nb_p_arbre][3]/(tree_cluster_leaves[nb_p_arbre][2]-1))*(nk_CH[k]-(tree_cluster_leaves[nb_p_arbre][2]-1)));
						clusterK_same[k]+=tree_cluster_leaves[nb_p_arbre][3];
					}
				}
			}
		}
		if(nk_CH[k]>1){
			FO_old += (clusterK_same[k]/(2.0*nk_CH[k]));
		}
	}

	Dref=FO_old;
	
	
	for (int i=1;i<=n; i++){
		  
		//Pour le cluster source
		k_source = list[i];
		nb_cluster_source = nk_CH[k_source];
		nb_cluster_source -=1;
		
		//on considere que l'element est dans un autre groupe (fictif son ancien cluster+100)
		list[i]=list[i]+100;		
		
		//Inferer la matrice d'adjacence du cluster k_source
		buildMatAdjecence(k_source,n,list,n_identique,mat_adjacence);
		
		//visited is initialized to zero
		for(int ii=0;ii<n;ii++)
			visited[ii]=0;
		
		//Va permettre de compter le nombre de point connecté
		nb_visited=0;

		aret_de_recherche = 0;
		depart_graphe_i = 0;
		while(depart_graphe_i<n && aret_de_recherche==0){
			depart_graphe_j = 0;
			while(depart_graphe_j<n && aret_de_recherche==0){
				if(mat_adjacence[depart_graphe_i][depart_graphe_j]!=2){
					aret_de_recherche = 1;
				}
				depart_graphe_j++;
			}
			depart_graphe_i++;
		}
		
		//Regarder si le graphe du cluster k_source est connexe sans l'element i
		DFS(depart_graphe_i-1,visited,mat_adjacence,n,nb_visited);		
		list[i]=list[i]-100;
		//Si le cluster source est tjs connexe en l'abs de l'arbre i
		if(nb_cluster_source==nb_visited){
			if(nk_CH[list[i]]>1){
				for (int k=1;k<=kk;k++){
					/* cout<<"Cluster # "<<k_source<<endl; */
					//Calcul de la distance RF de chaque point i
					// et assignation du point i au bon cluster
					// Compute a RF distance to the centroid k

					//test si le point i n'appartenait pas initiallement à k 
					//k_source!=k pour éviter de revifier un élément qui a changé de cluster avec son cluster d'origine
					if(list[i]!=k && k_source!=k && old_k!=k){
						/* cout<<"Vers le cluster "<<k<<endl; */
						nb_cluster_dest = nk_CH[k];
						nb_el_dest=0;
						trouve = 0;
						while(nb_el_dest<n && trouve==0){
							if(list[nb_el_dest+1]==k){
								if(n_identique[nb_el_dest][i-1]>3){
									trouve=1;
								}
							}
							nb_el_dest++;
						}
						//si l'element du cluster source a au moins un lien avec 
						//les élements du cluster de destination							
						if(trouve!=0){

							//retirer les anciens calculs de la FO_new
							if(nk_CH[k]>=1 && nk_CH[list[i]]>=1){
								FO_new = FO_old - (clusterK_same[k]/(2.0*nk_CH[k])) - (clusterK_same[list[i]]/(2.0*nk_CH[list[i]]));	
							}else if(nk_CH[list[i]]>=1){
								FO_new = FO_old - (clusterK_same[list[i]]/(2.0*nk_CH[list[i]]));
							}else if(nk_CH[k]>=1){
								FO_new = FO_old - (clusterK_same[k]/(2.0*nk_CH[k]));
							}else{
								FO_new = FO_old;
							}
							
							//compute la distance intra
							//Pour le cluster de destination
							nb_element_dest_i = 0;
							RF_element_dest_i = 0;
							tmp_calc_dest = 0;
							p[k]=0;
							
							//calcul de la nouvelle distance RF de l'element i 
							//avec les autres elements du cluster k
							for (int i_dest=1;i_dest<=n;i_dest++){
								if (list[i_dest]==k){
									if(n_identique[i-1][i_dest-1]>3){
										RF_element_dest_i+=(2.0*mat[i-1][i_dest-1]);
										nb_element_dest_i++;
									}else{
										//compte le nombre de paire d'arbres i et j ayant un nombre de feuilles chevauchants inf à 3
										p[k]++;
									}
								}
							} 

							nb_cluster_dest +=1;
							if((nb_element_dest_i+1)>=2){
								if((nb_element_dest_i+1)!=nb_cluster_dest){
									RF_element_dest_i+=((RF_element_dest_i/(nb_element_dest_i))*(nb_cluster_dest-(nb_element_dest_i)));
									tmp_calc_dest+=(RF_element_dest_i);
									tmp_calc_dest+=clusterK_same[k];
								}else{
									tmp_calc_dest += (RF_element_dest_i);
									tmp_calc_dest += clusterK_same[k];
								}
							}
							tmp_calc_dest = arrondir(tmp_calc_dest,3);
							if(nb_cluster_dest>0){
								FO_new = FO_new + (tmp_calc_dest/(2.0*nb_cluster_dest));
							}
							//compute la distance intra
							//Pour le cluster source
							
 							tmp_calc_source = clusterK_same[k_source];
							for(int j=1;j<=n; j++){
								if(list[j]==k_source){
									if(n_identique[i-1][j-1]>3){
										tmp_calc_source -= (2.0*mat[i-1][j-1]);
									}
								}
							}	 							
							tmp_calc_source = arrondir(tmp_calc_source,3);
							
							if(nb_cluster_source>0){
								FO_new = FO_new + (tmp_calc_source/(2.0*nb_cluster_source));
							} 
							
							if(FO_new<FO_old){
								Dref=FO_new;		
								kref=k;
								k_source = k;
								
								//mise à jour de nk_CH[]
								nk_CH[k] = nb_cluster_dest;
								nk_CH[list[i]] = nb_cluster_source;
								
								howmany[k] = nb_cluster_dest;
								howmany[list[i]] = nb_cluster_source;

								//mise à jour de la distance intra de lelement i
								tree_cluster_leaves[i-1][1] = k;
								tree_cluster_leaves[i-1][2] = nb_element_dest_i;
								tree_cluster_leaves[i-1][3] = RF_element_dest_i;								
								
								//mise à jour de la distance intra-groupe des deux clusters modifiés 
								clusterK_same[list[i]] = tmp_calc_source;
								clusterK_same[k] = tmp_calc_dest;
								//mise à jour de la fonction objective FO_old
								FO_old = FO_new;
								//mise à jour la liste de distribution des elements
								old_k = list[i];
								list[i] = k;	
								new_k = k;
							}
						}				
					}
					
				}			
			}
		}
		SSE=SSE+Dref;         //SSE=SSE+Dref		 
		/* howmany[kref]++;	//howmany(kref)=howmany(kref)+1 */
	}	
	
	delete [] clusterK_same;
	delete [] nk_CH;
	delete [] p;
	delete [] visited;
	
	for (int i=0;i<n;i++)
	{
		delete [] mat_adjacence[i];
	}

	delete [] mat_adjacence;
	
    return Dref;
    
}//  end FO_super_tree

double arrondir(double num,int digits){
	return floor(num*pow(10,digits)+0.5)/pow(10,digits);   
}

void DFS(int i, int *visited, int **mat_adjacence, int n, int &nb_visited)
{
    int j;
    visited[i]=1;
	nb_visited++;
	
    for(j=0;j<n;j++){
       if(visited[j]==0&&mat_adjacence[i][j]==1){			
            DFS(j,visited,mat_adjacence,n,nb_visited);
		}
	}
}

void buildMatAdjacenceFusion(int k1, int k2, int n, int *list, double **n_identique, int **mat_adjacence)
{
	//Parcours des elements i 
	for(int i=0;i<n;i++){
		//Parcours des elements j
		for(int j=0;j<n;j++){
			//on verifie si les elements i et j appartiennent aux clusters k1 ou k2
			if((k1==list[i+1] && k1==list[j+1]) || (k2==list[i+1] && k2==list[j+1]) || (k1==list[i+1] && k2==list[j+1]) || (k2==list[i+1] && k1==list[j+1])){
				if(n_identique[i][j]<=3){
					mat_adjacence[i][j]=0;
				}else{
					mat_adjacence[i][j]=1;
				}
			//sinon mettre des 2 comme masque
			//mettre des 2 aux elements n'appartenant pas au cluster en cours k
			}else{
				mat_adjacence[i][j]=2;
			}
		}
	}
}

void buildMatAdjecence(int k, int n, int *list, double **n_identique, int **mat_adjacence)
{
	//Parcours des elements i 
	for(int i=0;i<n;i++){
		//Parcours des elements j
		for(int j=0;j<n;j++){
			//on verifie si les elements i et j appartiennent au meme cluster
			if(k==list[i+1] && k==list[j+1]){
				if(n_identique[i][j]<=3){
					mat_adjacence[i][j]=0;
				}else{
					mat_adjacence[i][j]=1;
				}
			//sinon mettre des 2 comme masque
			//mettre des 2 aux elements n'appartenant pas au cluster en cours k
			}else{
				mat_adjacence[i][j]=2;
			}
		}
	}
}

double DistanceCH(int &n,int &kmax,double** mat,int* list,double** Ww,double FO_new,double facteur) {
	
	double SSB=0.0;
	double SSW=0.0;
	double dist_all = 0.0;
	
	double distance_total = 0.0;
	int cluster_k=0;
	int *nk_CH = new int [kmax+1];
	int k_cluster = 0;
	
	for(int k=1;k<=kmax; k++){
		nk_CH[k]=0;
	}
	
	
	for(int k=1;k<=kmax; k++){
		nk_CH[list[k]]++;
	}
	
	for(int k=1;k<=kmax; k++){
		if(nk_CH[k]!=0){
			k_cluster++;
		}
	}
	
	double RF;
	
	//compute dist_all
	for (int i=1;i<n;i++){
		for (int j=i+1;j<=n;j++){
			RF = mat[i-1][j-1];
			dist_all += RF;
		}
	}
	
	dist_all = dist_all*facteur;
	
	//compute SSW
	SSW = FO_new;
	
	//compute SSB
	SSB = dist_all - SSW;
	
 	if((fabs(SSW)>0.000001) && (k_cluster>1)){
		distance_total=(SSB/SSW)*((n-k_cluster)/(k_cluster-1.0));
	}
	else if(fabs(SSW)<=0.000001  && (k_cluster>1)){
		distance_total=10000000.0*SSB*((n-k_cluster)/(k_cluster-1.0));
	}

	delete [] nk_CH;

	return distance_total;
      
	
    
}//  end ************************End of Distances DistanceCH

double FO_borne_supp(int &n,int &kmax,double** mat,double* Dvec,int* list,int* howmany,double &SSE,int &kk,vector <string> monTableau)
{
	double *clusterK_same = new double [kmax+1];
	int *nk_CH = new int [kmax+1];
	int cluster_k=0;
	double RF = 0.0;
	double Dref=0,D1=0;		//Real*8 Dref,D1,SSE,Dvec(kmax),weight(pmax)
	int	kref=0;		//Integer list(nmax),howmany(kmax),kref
	//Integer ishort(pmax)
	// Compute squared distances to group centroids. Assign objects to nearest one
	SSE=0;		//SSE=0.0

	int k_source = 0;
	int new_k = 0;
	int old_k = 0;
	int nb_cluster_dest = 0;   
	int nb_cluster_source = 0;  
	double FO_old = 0.0;
	double FO_new = 0.0;
	double tmp_calc = 0.0;
	double tmp_calc_dest = 0.0;
	double tmp_calc_source = 0.0;
	
	//initialisation des variables
	for(int k=1;k<=kmax; k++){
		nk_CH[k]=0;
		clusterK_same[k]=0.0;
	}
	
	//Compter le nombre d'éléments par cluster
	for(int k=1;k<=kmax; k++){
		nk_CH[list[k]]++;
	}
	
	//compute for each cluster initially, SSW value (intra groupe distance)
	//compute SSW
	for (int i=1;i<n;i++){
		cluster_k=list[i];
		for (int j=i+1;j<=n;j++){
			if (list[j]==cluster_k){
				RF = mat[i-1][j-1];
				clusterK_same[cluster_k]+=RF;
			}
		}
	} 
	
	for (int k=1;k<=kk;k++){
		if(nk_CH[k]>1){
			FO_old += ((clusterK_same[k]*2.0)/(nk_CH[k]));
		}
	}
	Dref = FO_old;
	
	for (int i=1;i<=n; i++)			//do 20 i=1,n
	{
		k_source = list[i];
		if(nk_CH[list[i]]>1){
			for (int k=1;k<=kk;k++)			//do 12 k=1,kk
			{		
				//Calcul de la distance RF de chaque point i
				// et assignation du point i au bon cluster
				// Compute a RF distance to the centroid k

				//test si le point i n'appartenait pas initiallement à k 
				//k_source!=k pour éviter de revifier un élément qui a changé de cluster avec son cluster d'origine
				if(list[i]!=k && k_source!=k){
					nb_cluster_dest = nk_CH[k];
					nb_cluster_source = nk_CH[list[i]];
					
					if(nk_CH[k]>1 && nk_CH[list[i]]>1){
						FO_new = FO_old - ((clusterK_same[k]*2.0)/(nk_CH[k])) - ((clusterK_same[list[i]]*2.0)/(nk_CH[list[i]]));	
					}else if(nk_CH[list[i]]>1){
						FO_new = FO_old - ((clusterK_same[list[i]]*2.0)/(nk_CH[list[i]]));
					}else if(nk_CH[k]>1){
						FO_new = FO_old - ((clusterK_same[k]*2.0)/(nk_CH[k]));
					}else{
						FO_new = FO_old;
					}
					
					//compute Function objective
					//Pour le cluster de destination
					tmp_calc_dest = clusterK_same[k];
					for(int j=1;j<=n; j++){
						if(list[j]==k){
							tmp_calc_dest += mat[i-1][j-1];
						}
					}
					nb_cluster_dest +=1;
					tmp_calc_dest = arrondir(tmp_calc_dest,3);
					if(nb_cluster_dest>1){
						FO_new = FO_new + ((tmp_calc_dest*2.0)/(nb_cluster_dest));
					}
					
					//Pour le cluster source
					tmp_calc_source = clusterK_same[list[i]];
					for(int j=1;j<=n; j++){
						if(list[j]==list[i]){
							tmp_calc_source -= mat[i-1][j-1];
						}
					}
					nb_cluster_source -=1;
					
					tmp_calc_source = arrondir(tmp_calc_source,3);
					if(nb_cluster_source>1){
						FO_new = FO_new + ((tmp_calc_source*2.0)/(nb_cluster_source));
					}
					
					if(FO_new<=FO_old){
						Dref=FO_new;		
						kref=k;
						
						//mise à jour de nk_CH[]
						nk_CH[k] = nb_cluster_dest;
						nk_CH[list[i]] = nb_cluster_source;				
						
						//mise à jour de la distance intra-groupe des deux clusters modifiés 
						clusterK_same[list[i]] = tmp_calc_source;
						clusterK_same[k] = tmp_calc_dest;
						
						//mise à jour de la fonction objective FO_old
						FO_old = FO_new;
						
						//mise à jour la liste de distribution des elements
						list[i] = k;	
						
						//A VOIR SI UTILE
						new_k = k;
						old_k = list[i];
					}else{
						Dref=FO_old;
					}
									
				}
				
			}		
		}
		SSE=SSE+Dref;         //SSE=SSE+Dref		 
		howmany[kref]++;	//howmany(kref)=howmany(kref)+1
		
	}	
	
	delete [] clusterK_same;
	delete [] nk_CH;
	
    return Dref;
    
}//  end FO_borne_supp


double FO_borne_avg(int &n,int &kmax,double** mat,double* Dvec,int* list,int* howmany,double &SSE,int &kk,vector <string> monTableau, double ** n_identique)
{
	double *clusterK_same = new double [kmax+1];
	int *nk_CH = new int [kmax+1];
	int cluster_k=0;
	double RF = 0.0;
	double Dref=0,D1=0;		//Real*8 Dref,D1,SSE,Dvec(kmax),weight(pmax)
	int	kref=0;		//Integer list(nmax),howmany(kmax),kref
	//Integer ishort(pmax)
	// Compute squared distances to group centroids. Assign objects to nearest one
	SSE=0;		//SSE=0.0

	int k_source = 0;
	int new_k = 0;
	int old_k = 0;
	int nb_cluster_dest = 0;   
	int nb_cluster_source = 0;  
	double FO_old = 0.0;
	double FO_new = 0.0;
	double FO_old_supp = 0.0;
	double FO_new_supp = 0.0;
	double FO_old_inf = 0.0;
	double FO_new_inf = 0.0;	
	double tmp_calc = 0.0;
	double tmp_calc_dest = 0.0;
	double tmp_calc_source = 0.0;
	
	//initialisation des variables
	for(int k=1;k<=kmax; k++){
		nk_CH[k]=0;
		clusterK_same[k]=0.0;
	}
	
	//Compter le nombre d'éléments par cluster
	for(int k=1;k<=kmax; k++){
		nk_CH[list[k]]++;
	}
	
	//compute for each cluster initially, SSW value (intra groupe distance)
	//compute SSW
	for (int i=1;i<n;i++){
		cluster_k=list[i];
		for (int j=i+1;j<=n;j++){
			if (list[j]==cluster_k){
				RF = mat[i-1][j-1];
				clusterK_same[cluster_k]+=RF;
			}
		}
	} 
	
	for (int k=1;k<=kk;k++){
		if(nk_CH[k]>1){
			FO_old_supp += (clusterK_same[k]*(2.0/nk_CH[k]));
			FO_old_inf += (clusterK_same[k]/(nk_CH[k]-1.0));
		}
	}
  
	for (int i=1;i<=n; i++)			//do 20 i=1,n
	{
		k_source = list[i];
		if(nk_CH[list[i]]>1){
			for (int k=1;k<=kk;k++)			//do 12 k=1,kk
			{		
				//Calcul de la distance RF de chaque point i
				// et assignation du point i au bon cluster
				// Compute a RF distance to the centroid k

				//test si le point i n'appartenait pas initiallement à k 
				//k_source!=k pour éviter de revifier un élément qui a changé de cluster avec son cluster d'origine
				if(list[i]!=k && k_source!=k){
					nb_cluster_dest = nk_CH[k];
					nb_cluster_source = nk_CH[list[i]];
					
					//calcul de SSW supp et inf
					if(nk_CH[k]>1 && nk_CH[list[i]]>1){
						FO_new_supp = FO_old_supp - (clusterK_same[k]*(2.0/nk_CH[k])) - (clusterK_same[list[i]]*(2.0/nk_CH[list[i]]));
						FO_new_inf = FO_old_inf - (clusterK_same[k]/(nk_CH[k]-1.0)) - (clusterK_same[list[i]]/(nk_CH[list[i]]-1.0));
					}else if(nk_CH[list[i]]>1){
						FO_new_supp = FO_old_supp - (clusterK_same[list[i]]*(2.0/nk_CH[list[i]]));
						FO_new_inf = FO_old_inf - (clusterK_same[list[i]]/(nk_CH[list[i]]-1.0));
					}else if(nk_CH[k]>1){
						FO_new_supp = FO_old_supp - (clusterK_same[k]*(2.0/nk_CH[k]));
						FO_new_inf = FO_old_inf - (clusterK_same[k]/(nk_CH[k]-1.0));
					}else{
						FO_new_supp = FO_old_supp;
						FO_new_inf = FO_old_inf;
					}
										
					//compute Function objective
					//Pour le cluster de destination
					tmp_calc_dest = clusterK_same[k];
					for(int j=1;j<=n; j++){
						if(list[j]==k){
							tmp_calc_dest += mat[i-1][j-1];
						}
					}
					nb_cluster_dest +=1;
					tmp_calc_dest = arrondir(tmp_calc_dest,3);
					if(nb_cluster_dest>1){
						FO_new_supp = FO_new_supp + (tmp_calc_dest*(2.0/nb_cluster_dest));
						FO_new_inf = FO_new_inf + (tmp_calc_dest/(nb_cluster_dest-1.0));
					}
					
					//Pour le cluster source
					tmp_calc_source = clusterK_same[list[i]];
					for(int j=1;j<=n; j++){
						if(list[j]==list[i]){
							tmp_calc_source -= mat[i-1][j-1];
						}
					}
					nb_cluster_source -=1;
					tmp_calc_source = arrondir(tmp_calc_source,3);
					if(nb_cluster_source>1){
						FO_new_supp = FO_new_supp + (tmp_calc_source*(2.0/nb_cluster_source));
						FO_new_inf = FO_new_inf + (tmp_calc_source/(nb_cluster_source-1.0));
					}
					
					FO_new = (FO_new_supp + FO_new_inf)/2.0;
					FO_old = (FO_old_supp + FO_old_inf)/2.0;
					
					if(FO_new<=FO_old){
						Dref=FO_new;		
						kref=k;
						
						//mise à jour de nk_CH[]
						nk_CH[k] = nb_cluster_dest;
						nk_CH[list[i]] = nb_cluster_source;					
						
						//mise à jour de la distance intra-groupe des deux clusters modifiés 
						clusterK_same[list[i]] = tmp_calc_source;
						clusterK_same[k] = tmp_calc_dest;
						
						//mise à jour de la fonction objective FO_old
						FO_old_supp = FO_new_supp;
						FO_old_inf = FO_new_inf;
						
						//mise à jour la liste de distribution des elements
						list[i] = k;	
						
						//A VOIR SI UTILE
						new_k = k;
						old_k = list[i];
					}else{
						Dref=FO_old;
					}
									
				}
			}		
		}
		SSE=SSE+Dref;         //SSE=SSE+Dref		 
		howmany[kref]++;	//howmany(kref)=howmany(kref)+1
		
	}
	
	delete [] clusterK_same;
	delete [] nk_CH;
    return Dref;
    
}//  end FO_borne_avg


double FO_borne_inf(int &n,int &kmax,double** mat,double* Dvec,int* list,int* howmany,double &SSE,int &kk,vector <string> monTableau, double ** n_identique)
{
	double *clusterK_same = new double [kmax+1];
	int *nk_CH = new int [kmax+1];
	int cluster_k=0;
	double RF = 0.0;
	double Dref=0,D1=0;		//Real*8 Dref,D1,SSE,Dvec(kmax),weight(pmax)
	int	kref=0;		//Integer list(nmax),howmany(kmax),kref
	//Integer ishort(pmax)
	// Compute squared distances to group centroids. Assign objects to nearest one
	SSE=0;		//SSE=0.0

	int k_source = 0;
	int new_k = 0;
	int old_k = 0;
	int nb_cluster_dest = 0;   
	int nb_cluster_source = 0;  
	double FO_old = 0.0;
	double FO_new = 0.0;
	double tmp_calc = 0.0;
	double tmp_calc_dest = 0.0;
	double tmp_calc_source = 0.0;
	
	//initialisation des variables
	for(int k=1;k<=kmax; k++){
		nk_CH[k]=0;
		clusterK_same[k]=0.0;
	}
	
	//Compter le nombre d'éléments par cluster
	for(int k=1;k<=kmax; k++){
		nk_CH[list[k]]++;
	}
	
	//compute for each cluster initially, SSW value (intra groupe distance)
	//compute SSW
	for (int i=1;i<n;i++){
		cluster_k=list[i];
		for (int j=i+1;j<=n;j++){
			if (list[j]==cluster_k){
				RF = mat[i-1][j-1];
				clusterK_same[cluster_k]+=RF;
			}
		}
	} 
	
	for (int k=1;k<=kk;k++){
		if(nk_CH[k]>1){
			FO_old += (clusterK_same[k]/(nk_CH[k]-1.0));
		}
	}
	
	for (int i=1;i<=n; i++)			//do 20 i=1,n
	{
		k_source = list[i];
		if(nk_CH[list[i]]>1){
			for (int k=1;k<=kk;k++)			//do 12 k=1,kk
			{		
				//Calcul de la distance RF de chaque point i
				// et assignation du point i au bon cluster
				// Compute a RF distance to the centroid k

				//test si le point i n'appartenait pas initiallement à k 
				//k_source!=k pour éviter de revifier un élément qui a changé de cluster avec son cluster d'origine
				if(list[i]!=k && k_source!=k){
					nb_cluster_dest = nk_CH[k];
					nb_cluster_source = nk_CH[list[i]];
					
					if(nk_CH[k]>1 && nk_CH[list[i]]>1){
						FO_new = FO_old - (clusterK_same[k]/(nk_CH[k]-1.0)) - (clusterK_same[list[i]]/(nk_CH[list[i]]-1.0));	
					}else if(nk_CH[list[i]]>1){
						FO_new = FO_old - (clusterK_same[list[i]]/(nk_CH[list[i]]-1.0));
					}else if(nk_CH[k]>1){
						FO_new = FO_old - (clusterK_same[k]/(nk_CH[k]-1.0));
					}else{
						FO_new = FO_old;
					}
					
					//compute Function objective
					//Pour le cluster de destination
					tmp_calc_dest = clusterK_same[k];
					for(int j=1;j<=n; j++){
						if(list[j]==k){
							tmp_calc_dest += mat[i-1][j-1];
						}
					}
					nb_cluster_dest +=1;
					tmp_calc_dest = arrondir(tmp_calc_dest,3);
					if(nb_cluster_dest>1){
						FO_new = FO_new + (tmp_calc_dest/(nb_cluster_dest-1.0));
					}
					
					//Pour le cluster source
					tmp_calc_source = clusterK_same[list[i]];
					for(int j=1;j<=n; j++){
						if(list[j]==list[i]){
							tmp_calc_source -= mat[i-1][j-1];
						}
					}
					nb_cluster_source -=1;
					
					tmp_calc_source = arrondir(tmp_calc_source,3);
					if(nb_cluster_source>1){
						FO_new = FO_new + (tmp_calc_source/(nb_cluster_source-1.0));
					}
					
					if(FO_new<=FO_old){
						Dref=FO_new;		
						kref=k;
						
						//mise à jour de nk_CH[]
						nk_CH[k] = nb_cluster_dest;
						nk_CH[list[i]] = nb_cluster_source;						
						
						//mise à jour de la distance intra-groupe des deux clusters modifiés 
						clusterK_same[list[i]] = tmp_calc_source;
						clusterK_same[k] = tmp_calc_dest;
						
						//mise à jour de la fonction objective FO_old
						FO_old = FO_new;
						
						//mise à jour la liste de distribution des elements
						list[i] = k;
						
						//A VOIR SI UTILE
						new_k = k;
						old_k = list[i];
					}else{
						Dref=FO_old;
					}	
									
				}
				
			}		
		}
		SSE=SSE+Dref;         //SSE=SSE+Dref		 
		howmany[kref]++;	//howmany(kref)=howmany(kref)+1
		
	}
	
	
	
	delete [] clusterK_same;
	delete [] nk_CH;
	
    return Dref;
    
}//  end FO_borne_inf



//euclidean intParam=13

double FO_euclidean(int &n,int &kmax,double** mat,double* Dvec,int* list,int* howmany,double &SSE,int &kk,vector <string> monTableau, double ** n_identique)
{
	double *clusterK_same = new double [kmax+1];
	int *nk_CH = new int [kmax+1];
	int cluster_k=0;
	double RF = 0.0;
	double Dref=0,D1=0;		//Real*8 Dref,D1,SSE,Dvec(kmax),weight(pmax)
	int	kref=0;		//Integer list(nmax),howmany(kmax),kref
	//Integer ishort(pmax)
	// Compute squared distances to group centroids. Assign objects to nearest one
	SSE=0;		//SSE=0.0

	int k_source = 0;
	int new_k = 0;
	int old_k = 0;
	int nb_cluster_dest = 0;   
	int nb_cluster_source = 0;  
	double FO_old = 0.0;
	double FO_new = 0.0;
	double tmp_calc = 0.0;
	double tmp_calc_dest = 0.0;
	double tmp_calc_source = 0.0;
	
	//initialisation des variables
	for(int k=1;k<=kmax; k++){
		nk_CH[k]=0;
		clusterK_same[k]=0.0;
	}
	
	//Compter le nombre d'éléments par cluster
	for(int k=1;k<=kmax; k++){
		nk_CH[list[k]]++;
	}
	
	//compute for each cluster initially, SSW value (intra groupe distance)
	//compute SSW
	for (int i=1;i<n;i++){
		cluster_k=list[i];
		for (int j=i+1;j<=n;j++){
			if (list[j]==cluster_k){
				RF = mat[i-1][j-1];
				clusterK_same[cluster_k]+=RF;
			}
		}
	} 
	
	for (int k=1;k<=kk;k++){
		if(nk_CH[k]>1){
			FO_old += (clusterK_same[k]/nk_CH[k]);
		}
	}
	
	for (int i=1;i<=n; i++)			//do 20 i=1,n
	{
		k_source = list[i];
		if(nk_CH[list[i]]>1){
			for (int k=1;k<=kk;k++)			//do 12 k=1,kk
			{		
				//Calcul de la distance RF de chaque point i
				// et assignation du point i au bon cluster
				// Compute a RF distance to the centroid k

				//test si le point i n'appartenait pas initiallement à k 
				//k_source!=k pour éviter de revifier un élément qui a changé de cluster avec son cluster d'origine
				if(list[i]!=k && k_source!=k){
					nb_cluster_dest = nk_CH[k];
					nb_cluster_source = nk_CH[list[i]];
					
					if(nk_CH[k]>1 && nk_CH[list[i]]>1){
						FO_new = FO_old - (clusterK_same[k]/nk_CH[k]) - (clusterK_same[list[i]]/nk_CH[list[i]]);	
					}else if(nk_CH[list[i]]>1){
						FO_new = FO_old - (clusterK_same[list[i]]/nk_CH[list[i]]);
					}else if(nk_CH[k]>1){
						FO_new = FO_old - (clusterK_same[k]/nk_CH[k]);
					}else{
						FO_new = FO_old;
					}
					
					//compute Function objective
					//Pour le cluster de destination
					tmp_calc_dest = clusterK_same[k];
					for(int j=1;j<=n; j++){
						if(list[j]==k){
							tmp_calc_dest += mat[i-1][j-1];
						}
					}
					nb_cluster_dest +=1;
					tmp_calc_dest = arrondir(tmp_calc_dest,3);
					if(nb_cluster_dest>1){
						FO_new = FO_new + (tmp_calc_dest/nb_cluster_dest);
					}
					
					//Pour le cluster source
					tmp_calc_source = clusterK_same[list[i]];
					for(int j=1;j<=n; j++){
						if(list[j]==list[i]){
							tmp_calc_source -= mat[i-1][j-1];
						}
					}
					nb_cluster_source -=1;
					
					tmp_calc_source = arrondir(tmp_calc_source,3);
					if(nb_cluster_source>1){
						FO_new = FO_new + (tmp_calc_source/nb_cluster_source);
					}
					
					if(FO_new<=FO_old){
						Dref=FO_new;		
						kref=k;
						
						//mise à jour de nk_CH[]
						nk_CH[k] = nb_cluster_dest;
						nk_CH[list[i]] = nb_cluster_source;			
						
						//mise à jour de la distance intra-groupe des deux clusters modifiés 
						clusterK_same[list[i]] = tmp_calc_source;
						clusterK_same[k] = tmp_calc_dest;
						
						//mise à jour de la fonction objective FO_old
						FO_old = FO_new;
						
						//mise à jour la liste de distribution des elements
						list[i] = k;	
						
						//A VOIR SI UTILE
						new_k = k;
						old_k = list[i];
					}else{
						Dref=FO_old;
					}
									
				}
				
			}		
		}
		SSE=SSE+Dref;         //SSE=SSE+Dref		 
		howmany[kref]++;	//howmany(kref)=howmany(kref)+1
		
	}
	
	delete [] clusterK_same;
	delete [] nk_CH;
	
    return Dref;
    
}//  end FO_euclidean



double FO_gap_statistique(int &n,int &kmax,double** mat,double* Dvec,int* list,int* howmany,double &SSE,int &kk,vector <string> monTableau)
{
	double *clusterK_same = new double [kmax+1]; //la somme de tous les RF des arbres dans le meme cluster k
	int *nk_gap = new int [kmax+1]; //nb d'arbres dans le cluster k
	int cluster_k=0; 
	double RF = 0.0;
	double Dref=0,D1=0;		//Real*8 Dref,D1,SSE,Dvec(kmax),weight(pmax)
	int	kref=0;		//Integer list(nmax),howmany(kmax),kref
	//Integer ishort(pmax)
	// Compute squared distances to group centroids. Assign objects to nearest one
	SSE=0;		//SSE=0.0

	int k_source = 0;
	int new_k = 0;
	int old_k = 0;
	int nb_cluster_dest = 0;   
	int nb_cluster_source = 0;  
	double FO_old = 0.0;
	double FO_new = 0.0;
	double tmp_calc = 0.0;
	double tmp_calc_dest = 0.0;
	double tmp_calc_source = 0.0;
	int k_cluster = 0;
	
	//initialisation des variables
	for(int k=1;k<=kmax; k++){
		nk_gap[k]=0;
		clusterK_same[k]=0.0;
	}
	
	//Compter le nombre d'éléments par cluster
	for(int k=1;k<=kmax; k++){
		nk_gap[list[k]]++;
	}
	
	//compte le nombre de cluster
	for(int k=1;k<=kmax; k++){
		if(nk_gap[k]!=0){
			k_cluster++;
		}
	}
	
	//compute for each cluster initially, SSW value (intra groupe distance)
	//compute SSW
	for (int i=1;i<n;i++){
		cluster_k=list[i];
		for (int j=i+1;j<=n;j++){
			if (list[j]==cluster_k){
				RF = 2*mat[i-1][j-1];
				clusterK_same[cluster_k]+=RF;
			}
		}
	} 
	
	for (int k=1;k<=k_cluster;k++){
		FO_old += clusterK_same[k];
	}
	Dref = FO_old;
	
	for (int i=1;i<=n; i++)
	{
		k_source = list[i];
		
		if(nk_gap[list[i]]>1){
			for (int k=1;k<=kk;k++){		
				//Calcul de la distance RF de chaque point i
				// et assignation du point i au bon cluster
				// Compute a RF distance to the centroid k

				//test si le point i n'appartenait pas initiallement à k 
				//k_source!=k pour éviter de revifier un élément qui a changé de cluster avec son cluster d'origine
				if(list[i]!=k && k_source!=k){
					
										
					FO_new = FO_old - clusterK_same[k] - clusterK_same[list[i]];
					
					//compute Function objective
					//Pour le cluster de destination
					tmp_calc_dest = clusterK_same[k];
					for(int j=1;j<=n; j++){
						if(list[j]==k){
							tmp_calc_dest += (2.0*mat[i-1][j-1]);
						}
					}
					
					nb_cluster_dest = nk_gap[k];
					nb_cluster_dest +=1;
					
					tmp_calc_dest = arrondir(tmp_calc_dest,3);
					FO_new = FO_new + tmp_calc_dest;
												
					//Pour le cluster source
					nb_cluster_source = nk_gap[list[i]];
					nb_cluster_source -=1;	
					tmp_calc_source = clusterK_same[list[i]];
					for(int j=1;j<=n; j++){
						if(list[j]==list[i]){
							tmp_calc_source -= (2.0*mat[i-1][j-1]);
						}
					}
					tmp_calc_source = arrondir(tmp_calc_source,3);
					FO_new = FO_new + tmp_calc_source;
					
					if(FO_new<=FO_old){
						Dref=FO_new;		
						kref=k;
						
						//mise à jour de nk_gap[]
						nk_gap[k] = nb_cluster_dest;
						nk_gap[list[i]] = nb_cluster_source;				
						
						//mise à jour de la distance intra-groupe des deux clusters modifiés 
						clusterK_same[list[i]] = tmp_calc_source;
						clusterK_same[k] = tmp_calc_dest;
						
						howmany[k] = nb_cluster_dest;
						howmany[list[i]] = nb_cluster_source;
						
						SSE=SSE+Dref;         //SSE=SSE+Dref			
						
						//mise à jour de la fonction objective FO_old
						FO_old = FO_new;
							
						//mise à jour la liste de distribution des elements
						old_k = list[i];
						list[i] = k;	
						
						//A VOIR SI UTILE
						new_k = k;						
					}else{
						Dref=FO_old;
					}
									
				}
				
			}		
		} 		 
		
	}	
	delete [] clusterK_same;
	delete [] nk_gap;

    return Dref;
    
}//  end FO_gap_statistique

double DistanceGAP(int &n,int &kmax,double** mat,int* list,double** Ww,double FO_new,double facteur, double ** n_identique) {
	
	//déclaration et initialisation de variables
	double gap=0.0;
	double exp_star = 0.0;
	double VK=0.0;
	double tmp = 0.0;
	
	int *nk_gap = new int [kmax+1];
	int k_cluster = 0;
	double n_avg_identique = 0.0;
	
	for(int k=1;k<=kmax; k++){
		nk_gap[k]=0;
	}
	
	//compte le nombre d'élément par cluster
	for(int k=1;k<=kmax; k++){
		nk_gap[list[k]]++;
	}
	
	//compte le nombre de cluster
	for(int k=1;k<=kmax; k++){
		if(nk_gap[k]!=0){
			k_cluster++;
		}
	}

	
	
	int *n_avg_identique_k = new int [k_cluster+1];
	for(int k=0;k<=k_cluster; k++){
		n_avg_identique_k[k] = 0;
	}
	double RF;
	
	//compute VK
	VK = FO_new;
	
	//compute the average of number of leaves identical in the same cluster --- of all clusters ??
	for(int i=0; i<n; i++){
		for(int j=0;j<n; j++){
			if(list[i+1]==list[j+1]){
				n_avg_identique_k[list[i+1]] += n_identique[i][j];
			}	
		}
	}
	
	for(int i=1; i<=k_cluster; i++){
		if(nk_gap[i]!=0){
			n_avg_identique_k[i] = n_avg_identique_k[i]/(nk_gap[i]*nk_gap[i]);
		}
	}
	
	//compute exp_star
	for(int i=1; i<=k_cluster; i++){
		n_avg_identique += n_avg_identique_k[i];
	}
	
	if(k_cluster!=0){
		n_avg_identique = n_avg_identique/k_cluster;
	}
		
	tmp = (n*n_avg_identique)/(12.0);
	if(tmp!=0 && k_cluster!=0){
		exp_star = log(tmp) - (2.0/n_avg_identique)*log(k_cluster);
	}else if (tmp!=0){
		exp_star = /* log(tmp) -  */(2.0/n_avg_identique)*log(k_cluster);
	}else if(k_cluster!=0){
		exp_star = log(tmp) - (2.0/n_avg_identique)/* *log(k_cluster) */;
	}else{
		exp_star = /* log(tmp) -  */(2.0/n_avg_identique)/* *log(k_cluster) */;
	}
	//compute gap
	if(VK==0){
		gap = exp_star /* - log(0.0000000000001) */;
	}else{
		gap = exp_star - log(VK);
	}
	
	delete [] nk_gap;
	delete [] n_avg_identique_k;

	return gap;     
	
    
}//  end ************************End of Distances DistanceGAP

double DistanceW(int &n,int &kmax,double** mat,int* list,double** Ww,double FO_new,double facteur, double ** n_identique) 
{
	
	//déclaration et initialisation de variables
	double critere_w=0.0;
	double delta_W=FO_new;
	
	int *nk_W = new int [kmax+1];
	int k_cluster = 0;
	
	for(int k=1;k<=kmax; k++){
		nk_W[k]=0;
	}
	
	//compte le nombre d'élément par cluster
	for(int k=1;k<=kmax; k++){
		nk_W[list[k]]++;
	}
	
	//compte le nombre de cluster
	for(int k=1;k<=kmax; k++){
		if(nk_W[k]!=0){
			k_cluster++;
		}
	}

	if((n-k_cluster)!=0){
		critere_w = (1.0/(n-k_cluster))*delta_W;
	}
	
	delete [] nk_W;

	return critere_w;     
	
    
}//  end ************************End of Distances W


double DistanceW2(int &n,int &kmax,double** mat, int* list, double** Ww) {
	
	double dist_intragroupe = 0.0;
	
	double distance_total = 100000000.0;
	int cluster_k=0;
	double *clusterK_same = new double [kmax+1];
	int *nk_W = new int [kmax+1];
	int k_cluster = 0;


	for(int k=1;k<=kmax; k++){
		nk_W[k]=0;
		clusterK_same[k]=0.0;
	}
	

	for(int k=1;k<=kmax; k++){
		nk_W[list[k]]++;
	}
	
	for(int k=1;k<=kmax; k++){
		if(nk_W[k]!=0){
			k_cluster++;
		}
	}
	
	double RF;
	//compute dist_intragroupe
	for (int i=1;i<n;i++){
		cluster_k=list[i];
		for (int j=i+1;j<=n;j++){
			if (list[j]==cluster_k){
				RF = mat[i-1][j-1];
				clusterK_same[cluster_k]+=(RF*Ww[i-1][j-1]);
			}
		}
	}
	
	for(int k=1;k<=kmax; k++){
		if(nk_W[k]>=1.0){
			clusterK_same[k]*=(1.0/(nk_W[k]*nk_W[k]));
			clusterK_same[k]=clusterK_same[k]+1.0;
			if(clusterK_same[k]<0.0001){
				clusterK_same[k]=0;
			}else{
				clusterK_same[k]=log(clusterK_same[k]);
			}
			clusterK_same[k]=nk_W[k]*clusterK_same[k];
			dist_intragroupe+=clusterK_same[k];
		}
	}

	
 	if(k_cluster!=n){
		distance_total=(dist_intragroupe);
	}
		
	delete [] clusterK_same;
	delete [] nk_W;
	
    return distance_total;
    
}//  end ************************End of Distances W2


double FO_W(int &n,int &kmax,double** mat,double* Dvec,int* list,int* howmany,double &SSE,int &kk,vector <string> monTableau)
{
	double *clusterK_same = new double [kmax+1]; //la somme de tous les RF des arbres dans le meme cluster k
	int *nk_W = new int [kmax+1]; //nb d'arbres dans le cluster k
	int cluster_k=0; 
	double RF = 0.0;
	double Dref=0,D1=0;		//Real*8 Dref,D1,SSE,Dvec(kmax),weight(pmax)
	int	kref=0;		//Integer list(nmax),howmany(kmax),kref
	//Integer ishort(pmax)
	// Compute squared distances to group centroids. Assign objects to nearest one
	SSE=0;		//SSE=0.0

	int k_source = 0;
	int new_k = 0;
	int old_k = 0;
	int nb_cluster_dest = 0;   
	int nb_cluster_source = 0;  
	double FO_old = 0.0;
	double FO_new = 0.0;
	double tmp_calc = 0.0;
	double tmp_calc_dest = 0.0;
	double tmp_calc_source = 0.0;
	int k_cluster = 0;
	
	//initialisation des variables
	for(int k=1;k<=kmax; k++){
		nk_W[k]=0;
		clusterK_same[k]=0.0;
	}
	
	//Compter le nombre d'éléments par cluster
	for(int k=1;k<=kmax; k++){
		nk_W[list[k]]++;
	}
	
	//compte le nombre de cluster
	for(int k=1;k<=kmax; k++){
		if(nk_W[k]!=0){
			k_cluster++;
		}
	}
	
	//compute for each cluster initially, SSW value (intra groupe distance)
	//compute SSW
	for (int i=1;i<n;i++){
		cluster_k=list[i];
		for (int j=i+1;j<=n;j++){
			if (list[j]==cluster_k){
				RF = mat[i-1][j-1];
				clusterK_same[cluster_k]+=RF;
			}
		}
	} 
	
	for (int k=1;k<=k_cluster;k++){
		if(nk_W[k]>1){
			FO_old += ((2.0/(nk_W[k]*(nk_W[k]-1)))*clusterK_same[k]);
		}else if(nk_W[k]==1){
			FO_old +=0;
		}
	}
	Dref = FO_old;
	
	for (int i=1;i<=n; i++)
	{
		k_source = list[i];
		
		if(nk_W[list[i]]>1){
			for (int k=1;k<=kk;k++){		
				//Calcul de la distance RF de chaque point i
				// et assignation du point i au bon cluster
				// Compute a RF distance to the centroid k

				//test si le point i n'appartenait pas initiallement à k 
				//k_source!=k pour éviter de revifier un élément qui a changé de cluster avec son cluster d'origine
				if(list[i]!=k && k_source!=k){
					
					if(nk_W[k]>1 && nk_W[list[i]]>1){
						FO_new = FO_old - ((clusterK_same[k]*2.0)/(nk_W[k]*(nk_W[k]-1.0))) - ((2.0*clusterK_same[list[i]])/(nk_W[list[i]]*(nk_W[list[i]]-1.0)));	
					}else if(nk_W[list[i]]>1){
						FO_new = FO_old - ((2.0*clusterK_same[list[i]])/(nk_W[list[i]]*(nk_W[list[i]]-1.0)));
					}else if(nk_W[k]>1){
						FO_new = FO_old - ((clusterK_same[k]*2.0)/(nk_W[k]*(nk_W[k]-1.0)));
					}else{
						FO_new = FO_old;
					}
					
					//compute Function objective
					//Pour le cluster de destination
					tmp_calc_dest = clusterK_same[k];
					for(int j=1;j<=n; j++){
						if(list[j]==k){
							tmp_calc_dest += mat[i-1][j-1];
						}
					}
					
					nb_cluster_dest = nk_W[k];
					nb_cluster_dest +=1;
					tmp_calc_dest = arrondir(tmp_calc_dest,3);
					if(nb_cluster_dest>1){
						FO_new = FO_new + ((tmp_calc_dest*2.0)/(nb_cluster_dest*(nb_cluster_dest-1.0)));
					}
												
					//Pour le cluster source
					nb_cluster_source = nk_W[list[i]];
					nb_cluster_source -=1;	
					tmp_calc_source = clusterK_same[list[i]];
					for(int j=1;j<=n; j++){
						if(list[j]==list[i]){
							tmp_calc_source -= mat[i-1][j-1];
						}
					}
					tmp_calc_source = arrondir(tmp_calc_source,3);
					if(nb_cluster_source>1){
						FO_new = FO_new + ((tmp_calc_source*2.0)/(nb_cluster_source*(nb_cluster_source-1.0)));
					}
					if(FO_new<=FO_old){
						Dref=FO_new;		
						kref=k;
						
						//mise à jour de nk_W[]
						nk_W[k] = nb_cluster_dest;
						nk_W[list[i]] = nb_cluster_source;				
						
						//mise à jour de la distance intra-groupe des deux clusters modifiés 
						clusterK_same[list[i]] = tmp_calc_source;
						clusterK_same[k] = tmp_calc_dest;
						
						howmany[k] = nb_cluster_dest;
						howmany[list[i]] = nb_cluster_source;
						
						SSE=SSE+Dref;         //SSE=SSE+Dref			
						
						//mise à jour de la fonction objective FO_old
						FO_old = FO_new;
							
						//mise à jour la liste de distribution des elements
						old_k = list[i];
						list[i] = k;	
						
						//A VOIR SI UTILE
						new_k = k;						
					}else{
						Dref=FO_old;
					}
									
				}
				
			}		
		} 		 
		
	}	
	delete [] clusterK_same;
	delete [] nk_W;

    return Dref;
    
}//  end FO_W

double DistanceBH(int &n,int &nmax,int &p,int &pmax,int &kmax,double** mat,double** xbar, int* howmany, int &kk,int* list,double* weight,int* ishort,map<int,string>  mapIndicesTreesFinal) {
	
    
	std::map<int,string>::iterator it = mapIndicesTreesFinal.begin();
	int indice = 0;
	double distance_k = 0;
	double distance_total = 0;
	int cluster_k;
	int *clusters_k = new int [kmax+1];
	double *clusterK = new double [kmax+1];
	int Nk=0;
	
	for(int k=1;k<=kmax; k++){
		clusters_k[k]=0;
		clusterK[k]=0;
	}
	
	distance_k=0.0;

	for (int i=1;i<n;i++){
		cluster_k=list[i];
		//pour savoir si ce cluster a ete traite
		if(clusters_k[cluster_k]==0){
			Nk++;
			clusters_k[cluster_k]=1;
		}
		for (int j=i+1;j<=n;j++){
			if (list[j]==cluster_k) {
				distance_k+=mat[i-1][j-1];
				clusterK[cluster_k]+=mat[i-1][j-1];
			}
		}
		
	}
	
	for (int k=1;k<n;k++){
		if(howmany[k]>1){
			distance_total+=clusterK[k]*(2.0/(howmany[k]*(howmany[k]-1.0)));
		}
	}
	
   distance_total/=Nk;
   
   	delete [] clusters_k;
	delete [] clusterK;
	
   return distance_total;
      

    
}//  end ************************End of Distances BH


int fact (int valeur) {
  int resultatFact = 1;
  for(int i = 1; i <= valeur; i++) {
    resultatFact *= i;
  }
  return resultatFact;
}

void listePartition(int Strouve [],map<int,string>  mapIndicesTreesFinal,vector <string> indicesTrees) {
	
    
	std::map<int,string>::iterator it = mapIndicesTreesFinal.begin();
	int indiceGroupe=1;
	
	for ( it = mapIndicesTreesFinal.begin(); it != mapIndicesTreesFinal.end(); it++)
	{
		string ListeIndicesGroupe = it->second;
		int nbTabl = Split(indicesTrees, ListeIndicesGroupe, ',');
		
		for(int i(1); i<indicesTrees.size(); ++i){
			//conversion un string en int
			istringstream iss(indicesTrees[i]);
			// convertir en un int
			int nombre;
			iss >> nombre; // nombre vaut indicesTrees[i]
			
			Strouve[nombre-1]=indiceGroupe;
		}
		indiceGroupe++;
	}

    
}//  end ************************End listePartition

//////////////////////
//***********Euclidean
//////////////////////

// Compute a SQUARED (weighted) Euclidean distance to the centroid
void Euclidean(int &i,int &k,int &nmax,int &kmax,int &p,int &pmax,double** mat,double** xbar,double* weight,int* ishort,double &D1)
{
// Compute a SQUARED (weighted) Euclidean distance to the centroid
	int j=0;
      D1=0.0;			//D1=0.0
	  for (j=1;j<=p;j++)		// do 10 j=1,p
	  {
		  D1=D1+weight[ishort[j]]*((mat[i-1][ishort[j]-1]-xbar[k][j])*(mat[i-1][ishort[j]-1]-xbar[k][j]));		//   10 D1=D1+weight(ishort(j))*((mat(i,ishort(j))-xbar(k,j))**2)
	 }
      return;
}	//end
//************************End of Euclidean

// Compute a RF distance to the centroid
void RF_tree2tree(double &D1,string tree_i_k,string centro_k)
{
	D1=0.0;			//D1=0.0
	double distances[4];
	for(int jk = 0; jk<4; jk++){
		distances[jk] = 0.0;
	}

	main_hgt(tree_i_k,centro_k,distances);
	D1=distances[0];
    return;
}	//end
//************************End of Euclidean

void delete_space(string &str){
	/*string chaineSansEspace;
	 for (unsigned i=0; i<str.length(); i++) {
	   if(str.at(i) != " ")
		  chaineSansEspace += str[i];
	}
	str = chaineSansEspace;*/
	return;
}
void dist(int &i,int &k,double** mat,double &D1,int &n,int* list){
	
	 D1=0.0;	
	int cluster_k = list[i];
	for (int j=i+1;j<=n;j++){
		if (list[j]==cluster_k) {
			D1+=mat[i-1][j-1];
		}
	}
    return;
}	//end
//************************End of dist

//////////////////////
//***********Euclidean
//////////////////////

// Compute a Eucledian distance between point i and j 
double Euclidean_point(int &i, int&j, int &p,int &pmax,double** mat,double* weight,int* ishort)
{

    if (i==j) return 0.0;
    
          double sumSq=0.0;			//D1=0.0
	  for (int a=1;a<=p;a++)		// do 10 j=1,p
	  {
		  sumSq+=pow((weight[ishort[a]]*mat[i-1][ishort[a]-1]-weight[ishort[a]]*mat[j-1][ishort[a]-1]),2);		                  
	  }
      return sqrt(sumSq);
}	//end
//************************End of Euclidean


//////////////////////
//***********Permute
//////////////////////

void Permute(int &iseed,int &n,int &nmax,int *iordre)
{
// This subroutine permutes the first 'n' values of vector 'iordre' at random
// in an equiprobable way. This property has been checked through intensive
// simulations.
      int i=0,j=0,itemp=0,km1=0,m=0;			//Integer iseed,n,i,j,itemp,km1,m
      //Integer iordre(nmax)
      m=n;		//m=n
      km1=n-1;	//km1=n-1

	  for (i=1;i<=km1;i++)			//      do 10 i=1,km1
	  {
m8:		// j = 1 + rand(iseed)*m;		//    8    j = 1 + rand(iseed)*m
		  j=1+(rand()/32767.0)*m;

         if(j>m) goto m8;			//if(j.gt.m) goto 8
         itemp = iordre[m];			//itemp = iordre(m)
         iordre[m] = iordre[j];		//iordre(m) = iordre(j)
         iordre[j] = itemp;			//iordre(j) = itemp
         m = m-1;					//m = m-1
	  }			//   10 continue
	  
	  
      return;
}			//end
//************************End of Permute


//////////////////////
//***********ReadData1
//////////////////////

void ReadData1(int &n,int &nmax,int &p,int &pmax,double** mat/* ,double* coord */,int* ishort,double* weight,double* colsum,int &ntran, char* nameb, int N)
{
     int p1=0,p2=0;
     double rowsum=0,tsum=0;
     
     int i=0, j=0, nlines=0;
       int   nmat=2; // (orientation of data))
	 int iflag=0;//iflag=0
      
	  //Read matrix parameters 
	  n = N;
	  p = N;	
      //printf("\nData:\nn:%d p:%d\n", n,p);
      
      if(n>nmax)    //if(n.gt.nmax) Stop 
	  {	
		  printf ("Too many objects. Use a sample of objects or recompile program to increase nmax.");				//     +'Too many objects. Use a sample of objects or recompile program.'
		  exit(1);
	  }
		       
     if(p>pmax)    // if(p.gt.pmax) Stop 'Too many variables. Recompile the program.'
        {
                printf ("Too many variables. Use a sample of objects or recompile program to increase pmax.");				//     +'Too many objects. Use a sample of objects or recompile program.'
                exit(1);
        }

	
   switch (nmat)   //goto (10,14,18) nmat
   {
   case 1:
	   {
		  
		   break;//goto 22;
	   }
   case 2:
	   {
		   
		   break;//goto 22;
	   }

//        To read the file of QTC variables:
   case 3:
	   {
		   printf ("How many position (p1) and QTC variables (p2)?\n");			
		   printf ("E.g., p1 = 0 or 1 or 2; p2 = 166.  Answer with 2 numbers:");      
		   scanf("%d %d",&p1, &p2);		
 
		   if(p2>pmax)   
		   {
			 printf ("Too many variables. Recompile the program to increase pmax.");			
			 exit(1);
		   }

		   p=p2;
		   printf ("\n"); 
 
		   break;
	   }
   } 
   
   
   //fclose(Input1); 
		  
	  for (j=1;j<=p;j++)		
	  {
		  ishort[j]=j;		
		  weight[j]=1.0;			
	  }
   
       strcpy(nameb,"stat.csv");
	  if((Output4 = fopen(nameb,"a"))==NULL)
	  {
		  printf("\n%s: result file open failed...",nameb);
		  exit(1);
	  } 
      
}//end
//************************End of ReadData1


//////////////////////
//***********Standard
//////////////////////    
void Standard(int &n,int &nmax,int &kmax,int &p,int &pmax,double** mat,double** sx,double** sx2,double* vect,double** var,double* weight,int* ishort,int &istand)
{
	  double temp=0, dfln=0, dflnm1=0;
	  int nlines=0, j=0,i=0;
	  switch (istand)         //goto (10,30,50) istand
	  {
	  case 1:
		  {
                  // (1) Standardize: y(i,j)'=(y(i,j)-Mean(y))/StDev(y)
			  dfln=static_cast<float>(n);		//10 dfln=dfloat(n)
			  dflnm1=static_cast<float>(n-1);	//dflnm1=dfloat(n-1)
			  for (j=1;j<=p;j++)			//do j=1,p
			  {
				  sx[1][j]=0.0;			//sx(1,j)=0.0
				  sx2[1][j]=0.0;			//sx2(1,j)=0.0
			  }			//enddo

			  for (i=1;i<=n;i++)		// do i=1,n
			  {
				  for (j=1;j<=p;j++)		//do j=1,p
				  {
					  temp=mat[i-1][ishort[j]-1];		//temp=mat(i,ishort(j))
					  sx[1][j]=sx[1][j]+temp;			//sx(1,j)=sx(1,j)+temp
					  sx2[1][j]=sx2[1][j]+temp*temp;	//sx2(1,j)=sx2(1,j)+temp*temp
				  }		//enddo
			  }			//enddo
			  
			  for (j=1;j<=p;j++)			//do j=1,p
			  {
				  vect[j]=sx[1][j]/dfln;		//vect(j)=sx(1,j)/dfln
				  var[1][j]=sx2[1][j]-(sx[1][j]*sx[1][j]/dfln);		//var(1,j)=sx2(1,j)-(sx(1,j)*sx(1,j)/dfln)
				  var[1][j]=sqrt(var[1][j]/(dflnm1));				//var(1,j)=dsqrt(var(1,j)/dflnm1)
			  }			//enddo
                        //     Standardize the data      
			  for (i=1;i<=n;i++)			//do i=1,n
			  {
				  for (j=1;j<=p;j++)		//do j=1,p
				  {
					  mat[i-1][ishort[j]-1]=(mat[i-1][ishort[j]-1]-vect[j])/var[1][j];		//mat(i,ishort(j))=(mat(i,ishort(j))-vect(j))/var(1,j)
				  }		//            enddo
			  }		//         enddo

		  break;			//goto 24
		  }
	  case 2:
		  {
// 
// (2) Range: y' =  y(i,j)/yMax
// 'sx(1,j)' is used to store the maximum values of the variables.
   
			  for (j=1;j<=p;j++)			//30 do j=1,p
			  {
				  sx[1][j]=mat[1-1][j-1];			//sx(1,j)=mat(1,j)
			  }		//enddo
			  
			  for (i=2;i<=n;i++)			//do i=2,n
			  {
				  for (j=1;j<=p;j++)		//do j=1,p
				  {
					  temp=mat[i-1][ishort[j]-1];			//temp=mat(i,ishort(j))
					  if(temp>sx[1][j]) sx[1][j]=temp;		//if(temp.gt.sx(1,j)) sx(1,j)=temp
				  }			//enddo
			  }				//enddo

			  for (i=1;i<=n;i++)			//do i=1,n
			  {
				  for (j=1;j<=p;j++)		//do j=1,p
				  {
					  mat[i-1][ishort[j]-1]=mat[i-1][ishort[j]-1]/sx[1][j];		//mat(i,ishort(j))=mat(i,ishort(j))/sx(1,j)
				  }		//enddo
			  }			//enddo
		  break;			//goto 24
		  }
	  case 3:
		  {
// 
// (3) Range: y' = (y(i,j)-yMin)/(yMax-yMin)
// Matrix 'sx' is used to store the minimum and maximum values of the variables.
			  for (j=1;j<=p;j++)			//50 do j=1,p
			  {
				  sx[1][j]=mat[1-1][j-1];		//sx(1,j)=mat(1,j)
				  sx[2][j]=mat[1-1][j-1];		//sx(2,j)=mat(1,j)
			  }			//enddo

			  for (i=2;i<=n;i++)			//do i=2,n
			  {
				  for (j=1;j<=p;j++)		//do j=1,p
				  {
					  temp=mat[i-1][ishort[j]-1];		//temp=mat(i,ishort(j))
					  if(temp<sx[1][j]) sx[1][j]=temp;			//if(temp.lt.sx(1,j)) sx(1,j)=temp
					  if(temp>sx[2][j]) sx[2][j]=temp;			//if(temp.gt.sx(2,j)) sx(2,j)=temp
				  }			//enddo
			  }			//enddo
			  
			  for (i=1;i<=n;i++)			//do i=1,n
			  {
				  for (j=1;j<=p;j++)		//do j=1,p
				  {
					  mat[i-1][ishort[j]-1]=(mat[i-1][ishort[j]-1]-sx[1][j])/(sx[2][j]-sx[1][j]);		//mat(i,ishort(j))=(mat(i,ishort(j))-sx(1,j))/(sx(2,j)-sx(1,j))
            	  }		//enddo
			  }			//enddo
		  break;			//goto 24
		  }
	  }
// 
// Print a few lines of data
   
m24:	  printf ("Print how many lines of transformed data? (None: type 0)");		//24 write(*,*) 'Print how many lines of transformed data? (None: type 0)'
		  scanf	("%d", &nlines);			//read(*,*) nlines
		  
		  if((nlines<0)||(nlines>n)) goto m24;			//		  if((nlines.lt.0).or.(nlines.gt.n)) goto 24
		  if(nlines>0)									//if(nlines.gt.0) then
		  {
			  for (i=1;i<=nlines;i++)		// do 26 i=1,nlines
			  {
				  for (j=1;j<=p;j++)			//26    write(*,101) (mat(i,ishort(j)), j=1,p)
				  {
					  printf("%lf ", mat[i-1][ishort[j]-1]);
				  }
			  }
		  }			//endif
      
		  printf ("\n"); // write(*,*)
      return;

}//   end
//************************End of Standard


//////////////////////
//***********Transform
//////////////////////
void Transform(int &n,int &nmax,int &p,int &pmax,double** mat,double* weight,double* colsum,int &ntran)
		//      Subroutine Transform(n,nmax,p,pmax,mat,weight,colsum,ntran)
{      
    //Integer n,nmax,p,pmax
      double rowsum=0,tsum=0;		//Real*8 mat(nmax,pmax),weight(pmax),colsum(pmax),rowsum,tsum
	  int iflag=0,i=0,j=0;
// Transformation of species data before K-means partitioning?
      printf ("\n"); // write(*,*)
m2:   printf ("Transform species abundance data, in order to base\n"); // write(*,*) 'Transform species abundance data, in order to base'
      printf ("the analysis on a different distance measure?\n"); // write(*,*) 'the analysis on a different distance measure?'
      printf ("(see Legendre & Gallagher, manuscript)\n"); // write(*,*) '(see Legendre & Gallagher, manuscript)'
      printf ("\n"); // write(*,*)
      printf ("(0) No transformation (Euclidean distance preserved)\n"); // write(*,*) '(0) No transformation (Euclidean distance preserved)'
      printf ("(1) Chord distance\n"); // write(*,*) '(1) Chord distance'
      printf ("(2) Chi-square metric\n"); // write(*,*) '(2) Chi-square metric'
      printf ("(3) Chi-square distance\n"); // write(*,*) '(3) Chi-square distance'
      printf ("(4) Distance between species profiles\n"); // write(*,*) '(4) Distance between species profiles'
      printf ("(5) Hellinger distance\n"); // write(*,*) '(5) Hellinger distance'

	  scanf	("%d",&ntran);			//read(5,*) ntran


      if((ntran<0)||(ntran>5)) goto m2;		//      if((ntran.lt.0).or.(ntran.gt.5)) goto 2

      printf ("\n"); // write(*,*)
// Computing the transformed variables.
// The code is written in such a way as to avoid computing in advance
// a vector 'rowsum(i)' of length 'n'

      if(ntran==0) return;			//if(ntran.eq.0) return


      if((ntran==2)||(ntran==3)) 		//      if((ntran.eq.2).or.(ntran.eq.3)) then
	  {
		  tsum=0.0;			//tsum=0.0
		  for (j=1;j<=p;j++)		//do 14 j=1,p
		  {
			  colsum[j]=0.0;		//colsum(j)=0.0
			  for (i=1;i<=n;i++)		//do 14 i=1,n
			  {
				  tsum=tsum+mat[i-1][j-1];				//tsum=tsum+mat(i,j)
				  colsum[j]=colsum[j]+mat[i-1][j-1];		//   14    colsum(j)=colsum(j)+mat(i,j)
			  }
		  }

		  for (j=1;j<=p;j++)					//do 16 j=1,p
		  {
			  if(colsum[j]==0.0) weight[j]=0.0;		//if(colsum(j).eq.0.0) weight(j)=0.0
		  }			//   16    continue
	  }				//endif


	  switch (ntran)         //      goto (20,30,40,50,60) ntran
	  {
	  case 1:
		  {
// Chord distance
			  for (i=1;i<=n;i++)		//20 do 26 i=1,n
			  {
				  rowsum=0.0;			//rowsum=0.0
				  for (j=1;j<=p;j++)	//do 22 j=1,p
				  {
					  rowsum=rowsum+mat[i-1][j-1]*mat[i-1][j-1];			//22 rowsum=rowsum+mat(i,j)**2
				  }
				  
				  if(rowsum==0.0)							//if(rowsum.eq.0.0) then
				  {
					  printf ("Line	%d has total species abundance = 0", i);			//write(*,100) i
					  iflag=1;		//iflag=1
					  break;			//goto 26
				  }		//    endif
				  rowsum=sqrt(rowsum);		//rowsum=dsqrt(rowsum)
				
				  for (j=1;j<=p;j++)			//do 24 j=1,p
				  {
					  mat[i-1][j-1]=mat[i-1][j-1]/rowsum;			//24 mat(i,j)=mat(i,j)/rowsum
				  }
			  }			//26 continue
			  break;		//      goto 68
		  }

	  case 2:
		  {
// Chi-square metric
			  for (i=1;i<=n;i++)				//30 do 36 i=1,n
			  {
				  rowsum=0.0;					//rowsum=0.0
				  for (j=1;j<=p;j++)			//do 32 j=1,p
				  {
					  rowsum=rowsum+mat[i-1][j-1];	//32 rowsum=rowsum+mat(i,j)
				  }
				  
				  if(rowsum==0.0)				//if(rowsum.eq.0.0) then
				  {
					  printf ("Line	%d has total species abundance = 0", i);			//write(*,100) i
					  iflag=1;		//iflag=1
					  break;			//goto 36
				  }		//         endif

				  for (j=1;j<=p;j++)			//do 34 j=1,p
				  {
					  if(colsum[j]==0.0) mat[i-1][j-1]=mat[i-1][j-1]/(rowsum*sqrt(colsum[j]));	//if(colsum(j).ne.0.0) mat(i,j)=mat(i,j)/(rowsum*dsqrt(colsum(j)))
				  }			//   34 continue
			  }				//   36 continue
			  break;		//      goto 68
		  }
	  case 3:
		  {
// Chi-square distance
			  for (i=1;i<=n;i++)				//   40 do 46 i=1,n
			  {
				  rowsum=0.0;					//     rowsum=0.0
				  for (j=1;j<=p;j++)			//     do 42 j=1,p
				  {
					  rowsum=rowsum+mat[i-1][j-1];	//   42 rowsum=rowsum+mat(i,j)
				  }
				  
				  if(rowsum==0.0) 			//      if(rowsum.eq.0.0) then
				  {
					  printf ("Line	%d has total species abundance = 0", i);			// write(*,100) i
					  iflag=1;		//iflag=1
					  break;			//goto 46
				  }			//endif

				  for (j=1;j<=p;j++)			//do 44 j=1,p
				  {
					  if(colsum[j]==0.0) mat[i-1][j-1]=sqrt(tsum)*mat[i-1][j-1]/(rowsum*sqrt(colsum[j]));		// if(colsum(j).ne.0.0) mat(i,j)=dsqrt(tsum)*mat(i,j)/(rowsum*dsqrt(colsum(j)))
				  }		//   44 continue
			  }			//   46 continue
			  break;		//      goto 68
		  }
	  case 4:
		  {
// Distance between species profiles
			  for (i=1;i<=n;i++)			//50 do 56 i=1,n
			  {
				  rowsum=0.0;					//   rowsum=0.0
				  for (j=1;j<=p;j++)			//do 52 j=1,p
				  {
					  rowsum=rowsum+mat[i-1][j-1];	// 52 rowsum=rowsum+mat(i,j)
				  }

				  if(rowsum==0.0)				// if(rowsum.eq.0.0) then
				  {
					  printf ("Line	%d has total species abundance = 0", i);			// write(*,100) i
					  iflag=1;		//iflag=1
					  break;			//goto 56
				  }			//         endif

				  for (j=1;j<=p;j++)      //do 54 j=1,p
				  {
					  mat[i-1][j-1]=mat[i-1][j-1]/rowsum;			//54 mat(i,j)=mat(i,j)/rowsum
				  }
			  }		//   56 continue
			  break;		//   goto 68
		  }
	  case 5:
		  {
// Hellinger distance
			  for (i=1;i<=n;i++)			// 60 do 66 i=1,n
			  {
				  rowsum=0.0;					//    rowsum=0.0
				  for (j=1;j<=p;j++)			//do 62 j=1,p
				  {
					  rowsum=rowsum+mat[i-1][j-1];	//  62 rowsum=rowsum+mat(i,j)
				  }

				  if(rowsum==0.0)				//  if(rowsum.eq.0.0) then
				  {
					  printf ("Line	%d has total species abundance = 0", i);			//   write(*,100) i
					  iflag=1;		//iflag=1
					  break;			//   goto 66
				  }			//           endif
  
				  for (j=1;j<=p;j++)			// do 64 j=1,p
				  {
					  mat[i-1][j-1]=sqrt(mat[i-1][j-1]/rowsum);		//64 mat(i,j)=dsqrt(mat(i,j)/rowsum)
				  }
			  }							//66 continue
			  break;		//      68 continue
		  }
		  
}//???
		  

		  if(iflag==1) //if(iflag.eq.1) stop 'This is incompatible with the selected transformation.'
		  {
			  printf ("This is incompatible with the selected transformation.");
			  exit (1);
		  }

      return;			//return
      
}//************************End of Transform


int Split(vector<string>& vecteur, string chaine, char separateur)
{
	vecteur.clear();

	string::size_type stTemp = chaine.find(separateur);
	
	while(stTemp != string::npos)
	{
		vecteur.push_back(chaine.substr(0, stTemp));
		chaine = chaine.substr(stTemp + 1);
		stTemp = chaine.find(separateur);
	}
	
	vecteur.push_back(chaine);

	return vecteur.size();
}

void convert_list(int *&list, int n, vector<string> &centroid_k,vector<int> &centroid_k_pos, int &kk){
	map<int,int> indice_convert;
	int nouveau_indice = 0;
	string centroid_k_tmp = "";
	int centroid_k_pos_tmp = 0;
	
	for(int i=1; i<=n; i++){
		if(indice_convert.find(list[i])==indice_convert.end()){
			nouveau_indice++;
			indice_convert[list[i]] = nouveau_indice;
			
			//swap centroid
			centroid_k_tmp = centroid_k[list[i]];
			centroid_k[list[i]] = centroid_k[nouveau_indice];
			centroid_k[nouveau_indice] = centroid_k_tmp;
			
			//swap indice des centroids
			centroid_k_pos_tmp = centroid_k_pos[list[i]];
			centroid_k_pos[list[i]] = centroid_k_pos[nouveau_indice];
			centroid_k_pos[nouveau_indice] = centroid_k_pos_tmp;
		}
	}
	
	
	//mettre à jour le kk
	kk = nouveau_indice;
	
	for(int k_supp=kk+1; k_supp<centroid_k.size(); k_supp++){
	
		//supp des cendroid au dela de kk+1
		centroid_k.pop_back();
		
		//supp les indices des cendroids au dela de kk+1
		centroid_k.pop_back();
	
	}
	
	for(int i=1; i<=n; i++){
		list[i] = indice_convert[list[i]];
	}
}

double NouvelleCallRF (double nb_esp_identique){
	
	//Ancienne version
	/* system("./rf myTrees rf_output.txt rf_tmp.txt  rf_output.tmp rf_matrices.txt"); */
	
	//nouvelle version
	system("./rf myTrees rf_output.txt rf_tmp.txt rf_matrices.txt");
	
	//Recuperer le string du consensus des arbres du cluster i par la variable Ci
	ifstream fileCi;
	fileCi.open("rf_output.txt", ios::in);
	string RF_new = "";
	
	if(fileCi){
		// si l'ouverture a fonctionné
		string ligne;
		int sixieme = 0;
		// tant que l'on peut mettre la ligne dans "contenu"
		while(getline(fileCi, ligne)){
			sixieme ++;
			if(sixieme==6){
				size_t pos = ligne.find(" = ");
				RF_new = ligne.substr(pos+3); 
			}
		}

		fileCi.close();
	}else{
		cerr << "Impossible d'ouvrir le fichier !" << endl;
	}
	
	system("rm rf_*");
	//normalization de la distance RF par 2n - 6
	double RF_norm = atof(RF_new.c_str());
	RF_norm = RF_norm/(2*nb_esp_identique-6);
	return RF_norm;
}
///* 
// * Code from Found in http://cran.r-project.org/web/packages/cluster/
// * Original code by Francois Romain 
// */

