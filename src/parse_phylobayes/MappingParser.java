package parse_phylobayes;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.StringReader;
import java.text.DecimalFormat;

/**
 * This class reads multiple .map files and calculates / exports substitution lengths for each type of
 * nucleotide substitution by summarizing multiple sites (from .map files). 
 * There are N sites in an alignment. 
 * There are C trees / steps in a .map files (single site).
 * This class works for single site with 12 different types of substitutions (ungrouped) or 6 types of grouped substitutions,
 * or triplet site with 4 different types of substitutions (Scenario 1) or 9 different types of substituions (Scenario 2).
 * See more details in TripletSiteParser.java.
 * @author Hui-Jie
 *
 */
public class MappingParser {
	
	/** Number of sites*/
	private int N;
	/** Number of iterations / trees in a file*/
	private int C;
	/** Filename prefix */
	private String name;
    /** Number of branches */
	private int branchNum;
	/** Store theta (mu in the manuscript) */
	private double[][][] theta;
	/** Store theta (mu) for groups of changes. i.e. grouping by strand symmetry*/
	private double[][][] theta_ss;
	/** Store theta_ss average over all iterations */
	private double[][] theta_ss_bar;	
	/** Store variance of theta */
	private double[][][] theta_var;
	/** Store variance of theta_ss*/
	private double[][][] theta_ss_var;
    /** Store proportion of states*/
	private double[][][] propStates;
	/** Store time of states */
	private double[][][] timeStates;
	/** Store branch lengths */
	private double[][] br;
	/** Store number of changes */
	private int[][][] numberOfChanges;
	/** Store changes in each substitution group */
	private int[][][] changesInGroups;
	/** Store the tree structure to produce newick format for Multidivtime output file */
	private Tree tree[];
	/** Store the filename of the outgroup */
	private String outgroupFileName;
	/** The number of types of changes. 12 for SingleSiteParser, 4 for TripletParser */
	private int numTypeChanges;
	/** The number of types of states. 4 for SingleSiteParser, 3 for TripletParser */
	private int numTypeStates;
	/** Store the option */
	private int option;
	/** Store number of groups of changes. 6 for SingleSiteParser with strand symmetry
	 *                                     9 for TripletParser with strand symmetry
	 */
	private int numGroup;
	
	/**
	 * Constructor
	 * @param N
	 * @param C
	 * @param name
	 * @param option 1: SingleSiteParser, 2: TripletParser
	 * @param outgroupFileName
	 */
	public MappingParser(int N, int C, String name, int option, String outgroupFileName) {
		this.N = N;
		this.C = C;
		this.name = name;
		this.option = option;
		this.outgroupFileName = outgroupFileName;
		
		switch(option) {
		case 1: //SingleSiteParser
			this.numTypeChanges = 12;
			this.numTypeStates = 4;
			this.numGroup = 6;
			this.tree = new Tree[1];
			break;
		case 2: //TripletSiteParser
			this.numTypeChanges = 4;
			this.numTypeStates = 3;
			this.numGroup = 9;
			this.tree = new Tree[3];
			break;
		default:
			System.out.println("Error!!");
	}
		System.out.println(new File("").getPath().toString());
		
		String inputMap = name+"_0.map";

		try {
			// returns the ClassLoader object associated with this Class
	        //ClassLoader cLoader = this.getClass().getClassLoader();
	        // input stream
			//InputStream inStream = cLoader.getResourceAsStream("./parse_phylobayes/"+inputMap);
			//InputStream inStream = MappingParser.class.getClassLoader().getResourceAsStream(inputMap);
			InputStream inStream = this.getClass().getResourceAsStream(new File("../" + inputMap).getPath().toString());
			BufferedReader r = new BufferedReader(new InputStreamReader(inStream));
			//need to figure out the number of branches first so that 
			//i can declare the size of the array to store info
			String line = r.readLine();
			StringReader str = new StringReader(line);
            TreeParser tp = new TreeParser(str, outgroupFileName);
            //store the tree structure here
            //note that this tree shares the same node/branch numbering and ancestral with all other trees
            this.tree[0] = tp.tokenize();
            this.branchNum = tree[0].getNumBranches();
            r.close();
            inStream.close();
		} catch (IOException e) {
            e.printStackTrace();
        } 
		
		//12 = # of substitution types for single nucleotide / adapt it so that it works for triplet
		theta = new double[numTypeChanges][branchNum][C];
		theta_ss = new double[numGroup][branchNum][C];
		theta_ss_bar = new double[numGroup][branchNum];
		theta_var = new double[numTypeChanges][branchNum][C];
		theta_ss_var = new double[numGroup][branchNum][C];
		try {
			setUp();
		} catch (IOException e) {
			e.printStackTrace();
		}
		calculateTheta();
		calculateThetaSS();
		calculateThetaSSbar();
	}
	
	/**
	 * Get number of branches
     * @return branchNum
	 */
	public int getBranchNum() {
		return branchNum;
	}
	
	/**
	 * Compute theta (mu in the manuscript)
	 * theta = {theta_g} g = 1, ..., # of types of changes
	 * theta_g = {theta_lg} l = 1, 2, ..., # of branches
	 * theta_lg = {theta_lgc} c = 1, 2, ..., C = # of iterations in MCMC
	 * theta_lgc = # of substitutions of type g (a -> b) on branch l in iteration c (sum over all sites) 
	 *             divided by the proportion of time that state is a on branch l in iteration c (sum over all sites).
	 * @throws IOException 
	 *             
	 */
	public void setUp() throws IOException {
		
		// 12 = # of substitution types for single nucleotide / adapt it so that it works for triplet
		// 4 = # of nucleotide types for single site / adapt it so that it works for triplet
		numberOfChanges = new int[numTypeChanges][branchNum][C];
		changesInGroups = new int[numGroup][branchNum][C]; 
		timeStates = new double[numTypeStates][branchNum][C];
		br = new double[branchNum][C];
		propStates = new double[numTypeStates][branchNum][C];
				
		if (option == 1) { //SingleSiteParser, parse 1 file at a time
			String inputMap = null;
			File f = null;
			//process N files (sites)
			for (int i = 0; i < N; i++) {
				inputMap = name+"_"+i+".map";
				//InputStream inStream = MappingParser.class.getClassLoader().getResourceAsStream(inputMap);
				InputStream inStream = this.getClass().getResourceAsStream(new File("../" + inputMap).getPath().toString());
				BufferedReader r = new BufferedReader(new InputStreamReader(inStream));
				System.out.println("i="+i);
				
				String line = null;
				StringReader str = null;
				for (int j = 0; j < C; j++) {
					//read first tree
					try {
						line = r.readLine();
						str = new StringReader(line);
		                TreeParser tp = new TreeParser(str, outgroupFileName);
		                tree[0] = tp.tokenize();
		                
		                SiteParser parse = new SingleSiteParser(tree[0]);

		                for (int l = 0; l < branchNum; l++) {
		                	br[l][j] = parse.getBranchLengths()[l]; //br is the same for all sites in the same iteration
		                	// 4 = # of nucleotide types for single site / adapt it so that it works for triplet
		                	for (int k = 0; k < numTypeStates; k++) {
		                		numberOfChanges[k][l][j] += parse.getNumberOfChanges()[k][l];
		                		timeStates[k][l][j] += parse.getTimeOfStates()[k][l];
		                		propStates[k][l][j] += parse.getPropStates()[k][l];
		                	}//end k
		                	//12 = # of substitution types for single nucleotide / adapt it so that it works for triplet
		                	for (int k = numTypeStates; k < numTypeChanges; k++) {
		                		numberOfChanges[k][l][j]+=parse.getNumberOfChanges()[k][l];
		                	}//end k			                	
		                	for (int k = 0; k < numGroup; k++) {
		                		changesInGroups[k][l][j]+= parse.getChangesInGroups()[k][l];
		                	}
		                }//end l	
		             		           
		                	//read second tree, discard this tree;
                            //REMOVE THIS LINE IF PHYLOBAYES HAS BEEN CHANGED TO INCLUDE ONLY ONE MAPPING PER MCMC ITERATION
		                	line = r.readLine();
		                	//read "" and discard it.
		                	line = r.readLine();	                
					 } catch (IOException e) {
					    // TODO Auto-generated catch block
					    e.printStackTrace();
					 } 				
				}//end j
				try {
	                inStream.close();
	            } catch (IOException e) {
	                e.printStackTrace();
	            }
				r.close();
		    } //end i
			
			
			
		} else if (option == 2) { //TripletParser, parse 3 files at a time
			String inputMap = null;
			//File f = null;
			//move window from site 1 (index 0) to site N-2 (index N-3)
			for (int i = 0; i < (N-2); i++) {
				
				System.out.println("i="+i);
				String[] line = new String[3];
				StringReader[] str = new StringReader[3];
				BufferedReader[] r = new BufferedReader[3];
				InputStream[] inStream = new InputStream[3];
				for (int j = i; j < (i+3); j++) { //construct 3 bufferedreader
					inputMap = name+"_"+j+".map";
					//f = new File(inputMap);
					//r[j-i] = new BufferedReader(new FileReader(f));
					//inStream[j-i] = MappingParser.class.getClassLoader().getResourceAsStream(inputMap);
					inStream[j-i] = this.getClass().getResourceAsStream(new File("../" + inputMap).getPath().toString());					
					r[j-i] = new BufferedReader(new InputStreamReader(inStream[j-i]));
				}
				
				for (int j = 0; j < C; j++) { //read C trees for each file 
					try {
						for (int k = 0; k < 3; k++) { //repeat for each file (total 3 files at a time)
							//read first tree
							line[k] = r[k].readLine();
							str[k] = new StringReader(line[k]);
			                TreeParser tp = new TreeParser(str[k], outgroupFileName);
			                tree[k] = tp.tokenize();
			              //read second tree, discard this tree
                          //REMOVE THIS LINE IF PHYLOBAYES HAS BEEN CHANGED TO INCLUDE ONLY ONE MAPPING PER MCMC ITERATION
		                	line[k] = r[k].readLine();
		                	//read "" and discard it.
		                	line[k] = r[k].readLine();
						}
						SiteParser parse = new TripletParser(tree);
						
						for (int l = 0; l < branchNum; l++) {
							//br = new double[branchNum][C];
		                	br[l][j] = parse.getBranchLengths()[l]; 
		                	for (int k = 0; k < numTypeStates; k++) {
		                		//numberOfChanges = new int[4][trees[0].getNumBranches()];
		                		numberOfChanges[k][l][j] += parse.getNumberOfChanges()[k][l];
		                		//timeOfStates = new double[3][trees[0].getNumBranches()];
		                		timeStates[k][l][j] += parse.getTimeOfStates()[k][l];
		                		//propStates = new double[3][trees[0].getNumBranches()];
		                		propStates[k][l][j] += parse.getPropStates()[k][l];
		                	}//end k
		                	for (int k = numTypeStates; k < numTypeChanges; k++) {
		                		numberOfChanges[k][l][j]+=parse.getNumberOfChanges()[k][l];
		                	}//end k
		                	for (int k = 0; k < numGroup; k++) {
		                		changesInGroups[k][l][j]+= parse.getChangesInGroups()[k][l];
		                	}
		                }//end l	
						
					 } catch (IOException e) {
					    // TODO Auto-generated catch block
					    e.printStackTrace();
					 } 
					
				}//end j
				
				//close bufferedreader
				for (int k = 0; k < 3; k++) {
					r[k].close();
					inStream[k].close();
				}
				
		    } //end i
		}
			
	}
	
	/**
	 * Calculate theta (mu) and theta_var (variance of mu).
	 * compute theta after gathering info for all sites
	 * theta = # of changes on that branch divided by the proportion of time of the starting state.
	 * theta_var = theta / proportion of the time of the starting state.
	 * 
	 */
	public void calculateTheta() {
		//theta = new double[numTypeChanges][branchNum][C];
		//numberOfChanges = new int[numTypeChanges][branchNum][C];
		//timeStates = new double[numTypeStates][branchNum][C];
		//br = new double[branchNum][C];
		//propStates = new double[numTypeStates][branchNum][C];		
		if(option == 1) {
			for (int i = 0; i < numTypeChanges; i++) {
				for (int j = 0; j < branchNum; j++) {
					for (int k = 0; k < C; k++) {
						theta[i][j][k] = numberOfChanges[i][j][k]/propStates[i/3][j][k]; 
						theta_var[i][j][k] = theta[i][j][j]/propStates[i/3][j][k];
					}
				}
			}
		} else if (option == 2) {
			for (int j = 0; j < branchNum; j++) {
				for (int k = 0; k < C; k++) {
					for (int i = 0; i < 2; i++) {
						if (numberOfChanges[i][j][k] == 0 && (propStates[0][j][k] + propStates[1][j][k] ) == 0) {
							theta[i][j][k] = 0;
							theta_var[i][j][k] = 0;
						} else {
							theta[i][j][k] = numberOfChanges[i][j][k]/(propStates[0][j][k]+propStates[1][j][k]);
							theta_var[i][j][k] = theta[i][j][k] / (propStates[0][j][k]+propStates[1][j][k]);
						}
					}
					for (int i = 2; i < 4; i++) {
						if (numberOfChanges[i][j][k] == 0 && propStates[2][j][k] == 0) {
							theta[i][j][k] = 0;
							theta_var[i][j][k] = 0;
						} else {
							theta[i][j][k] = numberOfChanges[i][j][k]/propStates[2][j][k];
							theta_var[i][j][k] = theta[i][j][k] / propStates[2][j][k];
						}
					}
					
				}
			}
		}
		
	}
	
	/**
	 * calculate theta for groups of substitutions. grouping by strand symmetry.
	 * Option 1: 
	 * 	group 0: G->C, C->G, i.e. type 4 & 7
	 *  group 1: G->T, C->A, i.e. type 5 & 6
	 *  group 2: T->A, A->T, i.e. type 2 & 9
	 *  group 3: T->G, A->C, i.e. type 1 & 10
	 *  group 4: G->A, C->T, i.e. type 3 & 8
	 *  group 5: A->G, T->C, i.e. type 0 & 11 
	 *  state 0: A, state 1: G, state 2: C, state 3: T
	 *  
	 *  Option 2:
	 *  group 0-3: non-CpG transversion
	 *  group 4-5: non-CpG transition
	 *  group 6-7: CpG transversion
	 *  group 8: CpG transition
	 *  state 0: non-CpG C+G, state 1: non-CpG A+T, state 2: CpG
	 */
	public void calculateThetaSS() {
		for (int j = 0; j < branchNum; j++) {
			for (int k = 0; k < C; k++) {
				if (option == 1) {
					theta_ss[0][j][k] = changesInGroups[0][j][k]/(propStates[1][j][k]+propStates[2][j][k]);
					theta_ss[1][j][k] = changesInGroups[1][j][k]/(propStates[1][j][k]+propStates[2][j][k]);
					theta_ss[2][j][k] = changesInGroups[2][j][k]/(propStates[0][j][k]+propStates[3][j][k]);
					theta_ss[3][j][k] = changesInGroups[3][j][k]/(propStates[0][j][k]+propStates[3][j][k]);
					theta_ss[4][j][k] = changesInGroups[4][j][k]/(propStates[1][j][k]+propStates[2][j][k]);
					theta_ss[5][j][k] = changesInGroups[5][j][k]/(propStates[0][j][k]+propStates[3][j][k]);	
					theta_ss_var[0][j][k] = theta_ss[0][j][k]/(propStates[1][j][k]+propStates[2][j][k]);
					theta_ss_var[1][j][k] = theta_ss[1][j][k]/(propStates[1][j][k]+propStates[2][j][k]);
					theta_ss_var[2][j][k] = theta_ss[2][j][k]/(propStates[0][j][k]+propStates[3][j][k]);
					theta_ss_var[3][j][k] = theta_ss[3][j][k]/(propStates[0][j][k]+propStates[3][j][k]);
					theta_ss_var[4][j][k] = theta_ss[4][j][k]/(propStates[1][j][k]+propStates[2][j][k]);
					theta_ss_var[5][j][k] = theta_ss[5][j][k]/(propStates[0][j][k]+propStates[3][j][k]);	
				} else if (option == 2) {
					for (int i = 0; i < 4; i++) {
						if (changesInGroups[i][j][k] == 0 && propStates[i/2][j][k] == 0) {
							theta_ss[i][j][k] = 0;
							theta_ss_var[i][j][k] = 0;
						} else {
							theta_ss[i][j][k] = changesInGroups[i][j][k]/propStates[i/2][j][k];
							theta_ss_var[i][j][k] = theta_ss[i][j][k]/propStates[i/2][j][k];
						}
					}
					
					for (int i = 4; i< 6; i++) {
						if (changesInGroups[i][j][k] == 0 && propStates[i/5][j][k] == 0) {
							theta_ss[i][j][k] = 0;
							theta_ss_var[i][j][k] = 0;
						} else {
							theta_ss[i][j][k] = changesInGroups[i][j][k]/propStates[i/5][j][k];
							theta_ss_var[i][j][k] = theta_ss[i][j][k]/propStates[i/5][j][k];
						}
					}
					
					for (int i = 6; i < numGroup; i++) {
						if (changesInGroups[i][j][k] == 0 && propStates[2][j][k] == 0) {
							theta_ss[i][j][k] = 0;
							theta_ss_var[i][j][k] = 0;
						} else {
							theta_ss[i][j][k] = changesInGroups[i][j][k]/propStates[2][j][k];
							theta_ss_var[i][j][k] = theta_ss[i][j][k]/propStates[2][j][k];
						}
					}
					
					
				}
			}
		}
	}
	
	/**
	 * Compute theta_ss average over all iterations
	 */
	public void calculateThetaSSbar() {
		for (int g = 0; g < numGroup; g++) {
			for (int j = 0; j < branchNum; j++ ){
				for (int k = 0; k < C; k++) {
					theta_ss_bar[g][j] += theta_ss[g][j][k];
				}
				theta_ss_bar[g][j] = theta_ss_bar[g][j] / C;
			}
		}
	}
	
	/*
	 * Get the covariance matrix of theta_g
	 * g = the index of the type of substitutions
	 * g has the same order as `SUBSTITUTION_ORDER'
	 * @param g index
	 * @param option 0 for theta, 1 for theta_ss
	 * @return covariance matrix of a certain type
	 */
	public double[][] getCovarianceMatrix(int g, int option) {
		double[][] thetaG;
		double[][] thetaG_var;
		if (option == 0) {
			thetaG = theta[g];
			thetaG_var = theta_var[g];
		} else {
			thetaG = theta_ss[g];
			thetaG_var = theta_ss_var[g];
		}
		
		/**
		 * Compute average of theta or theta_ss
		 * Compute average of theta_var or theta_ss_var
		 */
		double[] theta_bar = new double[branchNum];
		double[] theta_var_bar = new double[branchNum];
		for (int i = 0; i < branchNum; i++) {
			for (int j = 0; j < C; j++) {
				theta_bar[i] += thetaG[i][j];
				theta_var_bar[i] += thetaG_var[i][j];
			}
			theta_bar[i] = theta_bar[i] / C;
			theta_var_bar[i] = theta_var_bar[i] / C;
		}
		
		double[][] cov = new double[branchNum][branchNum];
		
		//diagonal first, variance
		for (int i = 0; i < branchNum; i++) {
				double sum = 0;
				for (int k = 0; k < C; k++) {
					//theta = new double[12][branchNum][C];
					if (option == 0) {
						sum += Math.pow((theta[g][i][k] - theta_bar[i]),2);
					} else {
						sum += Math.pow((theta_ss[g][i][k] - theta_bar[i]),2);
					}
					
			
				}
				// var(theta) = var(E(theta|M_c)) + E(var(theta|M_c))
				cov[i][i] = sum / (C-1) + theta_var_bar[i];
		}
		
		//off diagonal, covariance
		for (int i = 0; i < branchNum; i++) {
			for (int j = 0; j < branchNum; j++) {
				if (i != j) {
					double sum = 0;
					for (int k = 0; k < C; k++) {
						if (option == 0){
							sum += theta[g][i][k] * theta[g][j][k];
						} else {
							sum += theta_ss[g][i][k] * theta_ss[g][j][k];
						}
					}
					cov[i][j] = (sum - theta_bar[i]*theta_bar[j]*C) / (C-1);
				}
			}
		}
		
		return cov;
	}
	
    /**
     * Return the number of changes for each substitution type (scenario 1) on branches for each MCMC iteration.
     * @return numberOfChanges
     */
	public int[][][] getNumberOfChanges() {
		return numberOfChanges;
	}
	
    /**
     * Return the number of changes for each substitution type (scenario 2) on branches for each MCMC iteration.
     * @return numberOfChanges
     */
	public int[][][] getChangesInGroups() {
		return changesInGroups;
	}
	
    /**
	 * Get proportion of time the second position is a CpG / non-CpG site for each branch for each MCMC iteration.
	 * @return propStates
	 */
	public double[][][] getPropStates() {
		return propStates;
	}
	
    /**
	 * Get branch lengths for each branch for each MCMC iteration.
	 * @return br
	 */
	public double[][] getBranchLengths() {
		return br;
	}

    /**
     * Return the total time in each states on each branch for each MCMC iteration
     * @return timeStates
     */
	public double[][][] getTimeStates() {
		return timeStates;
	}
	
    /**
     * Return theta (mu) for each substitution type and each branch for each MCMC iteration
     * @return theta
     */
	public double[][][] getTheta() {
		return theta;
	}
	
    /** 
     * Create output files: o.estb.type, o.estb.group and substitutionLength.txt
     * o.estb.type and o.estb.group can be used as the inputs to Multidivtime.
     *
     * o.estb.type for single nucleotide substitutions are for 12 different substitution types
     * o.estb.group for single nucleotide substitutions are for 6 different substitution types (combine strand symmetric substitution types)
     *
     * o.estb.type for triple site nucleotide substitutions are for 4 different substitution types (non-CpG/CpG transversion/transition)
     * o.estb.group for triplet site nucleotide substitutions are for 9 different substitution types (defined in manuscript)
     *
     * substitutionLength.txt file stores the estimated substitution lengths for each strand symmetric substitution type.
     * This makes it easier to report results.
     **/
	public void printOutput() {
		for (int g = 0; g < numTypeChanges; g++) {
			File output;
			output = new File("o.estb.type"+g);
			String newick = getNewickTree(theta[g]);
			PrintStream print = null;
			
		  	  try{
		  		print = new PrintStream(output);	 
		  		//print tree in newick format with branch lengths averaged over iterations
			  	  print.print(newick);
			  	  print.println();
			  	  print.print(printNodeInfo());
			  	  print.print(printCovarianceMatrix(g,0));
	  	  	  } catch(FileNotFoundException e){
	  	  	  	  System.out.println("Problem creating file!");	  
	  	  	  } finally {
	  	        if (print != null) print.close();
	  	    }
		  	  
		}
		
		for (int g = 0; g < numGroup; g++) {
			File output;
			output = new File("o.estb.group"+g);
			String newick = getNewickTree(theta_ss[g]);
			PrintStream print = null;
			
		  	  try{
		  		print = new PrintStream(output);
		  	//print tree in newick format with branch lengths averaged over iterations
			  	  print.print(newick);
			  	  print.println();
			  	  print.print(printNodeInfo());
			  	  print.print(printCovarianceMatrix(g,1));
		  		
	  	  	  } catch(FileNotFoundException e){
	  	  	  	  System.out.println("Problem creating file!");	  
	  	  	  } finally {
		  	        if (print != null) print.close();
		  	    }
		  	  
		}
		
		File output;
		output = new File("substitutionLength.txt");
		PrintStream print = null;
		try {
			print = new PrintStream(output);
			// print substitution lengths for branch 0 to the last branch to file
			for (int g = 0; g < numGroup; g++) {
				print.println("Substitution length type "+ g);
				for (int j = 0; j < branchNum; j++) {
					print.println(theta_ss_bar[g][j]);
				}
			}
		} catch (FileNotFoundException e) {
			System.out.println("Problem creating substitution length file!");
		} finally {
  	        if (print != null) print.close();
  	    }
			
	}
	
    /**
     * This function print some output to screen. 
     * It prints the number of substitutions per type per branch, the proportion of time of each type per branch,
     * and the sume of theta (mu) for each type per branch to screen.
     **/
	public void printToScreen(){
		for (int g = 0; g < numTypeChanges; g++) {
			for (int j = 0; j < branchNum; j++) {
				int count = 0;
				for (int k = 0; k < C; k++) {
					//numberOfChanges = new int[numTypeChanges][branchNum][C];
					count += numberOfChanges[g][j][k];
				}
				System.out.println("Type "+g+" changes on branch "+j+": "+count);
			}
		}
		
		//output # of substitutions per group per branch
		//9 = 6 groups (4 transversion groups, 2 transition groups) in non-CpG sites 
	    //    + 3 groups (2 transversion groups, 1 transition group) in CpG sites
		//4 groups of transversion    
	    // group 0: G->C, C->G
	    // group 1: G->T, C->A
	    // group 2: T->A, A->T
	    // group 3: T->G, A->C  
		//2 groups of transition
	    //group 0: G->A, C->T
	    //group 1: A->G, T->C 
		for (int g = 0; g < numGroup; g++) {
			for (int j = 0; j < branchNum; j++) {
				int count = 0;
				for (int k = 0; k < C; k++) {
					//numberOfChanges = new int[numTypeChanges][branchNum][C];
					count += changesInGroups[g][j][k];
				}
				System.out.println("Group "+g+" changes on branch "+j+": "+count);
			}
		}
		
		for (int g = 0; g < numTypeStates; g++) {
			for (int j = 0; j < branchNum; j++) {
				double prop = 0.0;
				for (int k = 0; k < C; k++) {
					//propStates = new double[numTypeStates][branchNum][C];
					prop += propStates[g][j][k];
				}
				System.out.println("Proportion of time in state "+g+" on branch "+j+": "+prop);
			}
		}
		//theta = new double[numTypeChanges][branchNum][C];
		
		for (int i = 0; i < numTypeChanges; i++) {
			for (int j  = 0; j < branchNum; j++) {
				double theta_k = 0;
				for (int k = 0; k < C; k++) {
					theta_k += theta[i][j][k];
				}
				System.out.println("Type "+i+" branch "+j+" theta sum: "+theta_k);
			}
		}
		
		for (int i = 0; i < numGroup; i++) {
			for (int j  = 0; j < branchNum; j++) {
				double theta_k = 0;
				for (int k = 0; k < C; k++) {
					theta_k += theta_ss[i][j][k];
				}
				System.out.println("Group "+i+" branch "+j+" theta_ss sum: "+theta_k);
			}
		}
	}

	/**
	 * Store newick tree with theta averaging over C iterations into a String.
	 * theta = new double[12][branchNum][C];
	 * @param theta[g], g = specific types of substitutions
	 */
	public String getNewickTree(double[][] thetaG) {
		//compute theta_g averaging over C
		double[] theta_bar = new double[branchNum];
		for (int i = 0; i < branchNum; i++) {
			for (int j = 0; j < C; j++) {
				theta_bar[i] += thetaG[i][j];
			}
			theta_bar[i] = theta_bar[i] / C;
		}
		
		String newick = inOrderNewick(tree[0].root, theta_bar);		
		
		return newick;
	}
	
	
	/**
	 * Recursive print newick format
	 * @param root cursor to traverse the tree
	 * @return newick tree
	 */
	public String inOrderNewick(TreeNode root, double[] theta) {
		if (!root.isLeaf()) {
	          String output = "";
	          output += "(";
	          output += inOrderNewick(root.firstChild(), theta);
	          output += ",";
	          output += inOrderNewick(root.lastChild(), theta);
	          output += ")";
	          if (root.isRoot()) {
	        	  output += ";";
	          } else {
	        	output += ":" + theta[root.getNodeNum()];  
	          }
	          return output;
	      } else {
	    	  
	          return root.getName() + ":" + theta[root.getNodeNum()];
	      }
	}
	
	/**
	 * Print node information. Used for Multidivtime input.
	 * 1. list of names and tip numbers
	 * 2. list of child1, child2, ..., parent
	 * @return String containing the info
	 */
	public String printNodeInfo() {
		String info = "list of names and tip numbers follows:\n";
		for (int i = 0; i < tree[0].getNumLeaves(); i++) {
			info += tree[0].getNodeByNodeNum(i).getName();
			info += "  " + i + "\n";
		}
		
		info += "list of child1, child2, ..., parent follows:\n";
		for (int i = tree[0].getNumLeaves(); i < tree[0].getTotalNodeCount(); i++) {
			info += tree[0].getNodeByNodeNum(i).getChild(0).getNodeNum() + " ";
			info += tree[0].getNodeByNodeNum(i).getChild(1).getNodeNum() + " ";
			info += i + "\n";
		}
		return info;
	}
	
	/**
	 * Print covariance matrix for type g substitution
	 * @param g index of the type of substitution
	 * @param option 0 for theta, 1 for theta_ss
	 * @return covariance matrix
	 */
	public String printCovarianceMatrix(int g, int option) {
		DecimalFormat formatter = new DecimalFormat("#.#################");
		String cov = "variance-covariance matrix follows:\n";
		double[][] covMatrix = getCovarianceMatrix(g, option);
		for (int i = 0; i < branchNum; i ++) {
			for (int j = 0; j < branchNum; j ++) {
				cov += formatter.format(covMatrix[i][j]) + " ";
			}
			cov += "\n";
		}
		
		return cov;
	}
	
	/**
	 * Main class
	 * @param args[0]: number of sites
	 * @param args[1]: number of iterations
	 * @param args[2]: prefix of filename
	 * @param args[3]: 1 - single site ; 2 - triplet CpG/non-CpG
	 * @param args[4]: outgroup filename, currently will produce incorrect tree topology if an ourgroup file is not given.
	 */
	public static void main(String args[]) {
		int N = Integer.parseInt(args[0]);
		int C = Integer.parseInt(args[1]);
		String name = args[2];
		int option = Integer.parseInt(args[3]);
		if(option == 1 || option == 2) {
			String outgroup = null;
			if (args.length == 5) {
				outgroup = args[4];
			} else {
				outgroup = "";
			}
			MappingParser parse = new MappingParser(N, C, name, option, outgroup);
			//System.out.println("Test!");
			parse.printOutput();
			parse.printToScreen();
		} else {
			System.out.println("Argument error.");
			System.out.println("option 1 : single site.");
			System.out.println("option 2 : CpG dinucleotide.");
			
		}
		
		
	}
	
	
}
