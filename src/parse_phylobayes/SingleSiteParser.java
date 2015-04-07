package parse_phylobayes;


/**
 * This class defines the order of different types of nucleotide substitutions.
 * Order: AG, AC, AT, GA, GC, GT, CA, CG, CT, TA, TG, TC
 * Number of types of nucleotide substitutions = 12.
 * Nucleotide types = {A, G, C, T}.
 * @author Hui-Jie Lee
 *
 */
public class SingleSiteParser implements SiteParser {
	
    /** Input tree **/
	private Tree tree;
    /** Store the number of changes for each type on each branch **/
	private int[][] numberOfChanges; //dim = 12 x (# of branches)
    /** Store the total time of each state (A, G, C, T) on each branch **/
	private double[][] timeOfStates; //dim = 4 x (# of branches)
    /** Store the branch lengths for each branch **/
	private double[] branchLengths; //dim = # of branches
    /** Store the proportion of time in each state on each branch **/
	private double[][] propStates; //dim = 4 x (# of branches)
	
	//												   0    1    2    3    4    5    6    7    8    9    10   11
	public static final String[] SUBSTITUTION_ORDER= {"AG","AC","AT","GA","GC","GT","CA","CG","CT","TA","TG","TC"};
	//                                                0   1   2   3
    public static final String[] NUCLEOTIDE_ORDER = {"A","G","C","T"};
    
    
    /**  
     *   group 0: G->C, C->G
     *   group 1: G->T, C->A
     *   group 2: T->A, A->T
     *   group 3: T->G, A->C
     *   group 4: G->A, C->T
     *   group 5: A->G, T->C
     */                                            0       1     2      3      4      5
    public static final String[] GROUP_ORDER = {"GCCG","GTCA","TAAT","TGAC","GACT","AGTC"};
	
    private int[][] changesInGroups; //dim = 6 x (# of branches)
    
    
	/**
	 * Constructor of SingleSiteParser
	 * @param tree input tree
	 */
	public SingleSiteParser(Tree tree) {
		super();
		this.tree = tree;
		
		numberOfChanges = new int[12][tree.getNumBranches()];
		timeOfStates = new double[4][tree.getNumBranches()];
		branchLengths = new double[tree.getNumBranches()];
		propStates = new double[4][tree.getNumBranches()];
		
		setNumberOfChanges(computeNumberOfChanges());
		setTimeOfStates(computeTimeOfStates());
		setBranchLengths(computeBranchLength());
		setPropStates(computePropState());
		setChangesInGroups(computeChangesInGroups());
	}
	
	/** Return the number of branches in the tree
	 * @return number of branches
	 */
	public int numberOfBranches() {
		return branchLengths.length;
	}
	
	/**
	 * Assign number of changes for each type of changes on each branch.
	 * @return numberOfChanges
	 */
	private int[][] computeNumberOfChanges() {
		int[][] changes = new int[12][tree.getNumBranches()];
		for (int i = 0; i < 12; i++) {
			for (int j = 0; j < tree.getNumBranches(); j++) {
				changes[i][j] = tree.getNumberOfChanges(SUBSTITUTION_ORDER[i].substring(0, 1),SUBSTITUTION_ORDER[i].substring(1, 2))[j];
			}
		}		
		return changes;
	}
	/**
	 * Assign changes in each substitution groups on each branch.
	 *  group 0: G->C, C->G, i.e. type 4 & 7
	 *  group 1: G->T, C->A, i.e. type 5 & 6
	 *  group 2: T->A, A->T, i.e. type 2 & 9
	 *  group 3: T->G, A->C, i.e. type 1 & 10
	 *  group 4: G->A, C->T, i.e. type 3 & 8
	 *  group 5: A->G, T->C, i.e. type 0 & 11 
	 * @return changesInGroups
	 */
	private int[][] computeChangesInGroups() {
		int[][] changes = new int[6][tree.getNumBranches()];
		for (int j = 0; j < tree.getNumBranches(); j++) {
			changes[0][j] = numberOfChanges[4][j] + numberOfChanges[7][j];
			changes[1][j] = numberOfChanges[5][j] + numberOfChanges[6][j];
			changes[2][j] = numberOfChanges[2][j] + numberOfChanges[9][j];
			changes[3][j] = numberOfChanges[1][j] + numberOfChanges[10][j];
			changes[4][j] = numberOfChanges[3][j] + numberOfChanges[8][j];
			changes[5][j] = numberOfChanges[0][j] + numberOfChanges[11][j];
		}
		
		return changes;
	}
	
	/**
	 * Assign times of each nucleotide for each type of changes on each branch.
	 * @return timeOfStates
	 */
	private double[][] computeTimeOfStates() {
		double[][] times = new double[4][tree.getNumBranches()];
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < tree.getNumBranches(); j++) {
				times[i][j] = tree.getTimeOfState(NUCLEOTIDE_ORDER[i])[j];
			}
		}		
		return times;
	}
	
	/**
	 * Assign branch lengths for each branch.
	 * @return branchLengths
	 */
	private double[] computeBranchLength() {
		double[] br = new double[tree.getNumBranches()];
		for (int i = 0; i < tree.getNumBranches(); i++) {
			br[i] = tree.getBranchLength()[i];
		}
		return br;
	}
	
    /** 
     * Compute the proportion of states on each branch.
     * @return prop
     */
	private double[][] computePropState() {
		double[][] prop = new double[4][tree.getNumBranches()];
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < tree.getNumBranches(); j++) {
				prop[i][j] = tree.getPropState(NUCLEOTIDE_ORDER[i])[j];
			}
		}		
		return prop;
	}
	
    /** 
     * Set proportion of states 
     * @param propStates to be set
     */
	public void setPropStates(double[][] propStates) {
		this.propStates = propStates;
	}
	
    /**
     * Return the proportion of states on each branch
     * @return propStates
     */
	public double[][] getPropStates() {
		return propStates;
	}

    /**
     * Return the number of changes for each substitution type (ungrouped) on branches
     * @return numberOfChanges
     */
	public int[][] getNumberOfChanges() {
		return numberOfChanges;
	}
	
    /**
     * Return the number of changes for each substitution type (grouped) on branches
     * @return changesInGroups
     */
	public int[][] getChangesInGroups() {
		return changesInGroups;
	}

    /**
     * Set the number of changes for each substitution type (ungrouped) on branches
     * @param number of changes to be set 
     */
	private void setNumberOfChanges(int[][] numberOfChanges) {
		this.numberOfChanges = numberOfChanges;
	}
	
    /**
     * Set the number of changes for each substitution type (ungrouped) on branches
     * @param number of changes to be set
     */
	private void setChangesInGroups(int[][] changesInGroups) {
		this.changesInGroups = changesInGroups;
	}

    /** 
     * Return the total time in each states on each branch
     * @return timeOfStates
     */
	public double[][] getTimeOfStates() {
		return timeOfStates;
	}

    /**
     * Set the total time in each states on each branch
     * @param timeOfStates
     */
	private void setTimeOfStates(double[][] timeOfStates) {
		this.timeOfStates = timeOfStates;
	}

    /**
     * Return the branch length on each branch
     * @return branchLengths
     */
	public double[] getBranchLengths() {
		return branchLengths;
	}

    /**
     * Set the branch length on each branch
     * @param branchLengths
     */
	private void setBranchLengths(double[] branchLengths) {
		this.branchLengths = branchLengths;
	}
	
    /** 
     * Return the input tree
     * @return tree
     */
	public Tree getTree() {
		return tree;
	}

	
	/**
	 * Return the `branch length' of different types of substitutions.
	 * `branch length' = expected # of substitutions of type g (from state a to b) given the proportion of the time
	 * that the state is a. 
	 * @return theta as a vector
	 */
	/*
	public double[][] getTheta() {
		double[][] theta = new double[12][tree.getNumBranches()];
		for (int i = 0; i <12; i++) {
			for (int j = 0; j < tree.getNumBranches(); j++) {
				theta[i][j] = numberOfChanges[i][j]/propStates[i/3][j];
			}
		}		
		return theta;
	}
	*/
	
	

}
