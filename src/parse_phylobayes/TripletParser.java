package parse_phylobayes;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;

/**
 * Parse three consecutive sites /triplet at a time.
 * There are three states: non-CpG C+G site, non-CpG A+T site, CpG C+G site
 * Denote XX' -> YY' as substitution from X to Y or its complementary, X' to Y'.
 * Two scenarios are consider here.
 * Scenario 1:
 * 4 types of changes: non-CpG transversion, non-CpG transition, CpG transversion, CpG transition.
 * Scenario 2: (report in manuscript)
 * non-CpG site: transversion: GC -> CG (Type 1), GC -> AT (Type 2), AT -> TA (Type 3), AT -> CG (Type4)
 *               transition:   GC -> AT (Type 5), AT -> GC (Type 6)
 * CpG site:     transversion: GC -> CG (Type 7), CG -> AT (Type 8)
 *               transition:   GC -> AT (Type 9)
 * 
 * @author Hui-Jie Lee
 *
 */
public class TripletParser implements SiteParser {
    
    // A->T,C ; T->A,G ; C->A,G ; G->C,T                 0     1     2     3
    public static final String[] TRANSVERSION_ORDER = {"ATC","TAG","CAG","GCT"};
    // A->G ; T->C ; C->T ; G->A                      0    1    2    3
    public static final String[] TRANSITION_ORDER = {"AG","TC","CT","GA"};
    
    /** 
     * 4 groups of transversion
     * group 0: G->C, C->G
     * group 1: G->T, C->A
     * group 2: T->A, A->T
     * group 3: T->G, A->C                                0      1      2      3
     **/
    public static final String[] TRANSVERSION_GROUP = {"GCCG","GTCA","TAAT","TGAC"};
    /**
     * 2 groups of transition
     * group 0: G->A, C->T
     * group 1: A->G, T->C                              0      1
     */
    public static final String[] TRANSITION_GROUP = {"GACT","AGTC"};

	
    /** Store three trees for three sites **/
	private Tree[] trees;
    
    /** Store number of changes in Scenario 1 
     *  dim = 4 x (# of branches)
     **/
	private int[][] numberOfChanges;
    
    /** Store number of changes in Scenario 2
     * dim = 9 x (# of branches)
     **/
    private int[][] changesInGroups;
    
    /** Store the total time spent in each state on branches
     *  dim = 3 x (# of branches)
     **/
	private double[][] timeOfStates;
    
    /** Store the branch lengths for each branch
     *  dim = # of branches
     **/
	private double[] branchLengths;
    
    /** Store the proportion of time in each state on branches
     *  dim = 3 x (# of branches)
     **/
	private double[][] propStates;
    
	/** Store path state as an array of arraylist of String, each String contains 3 characters (triplet). 
	 *  There are branchNum of such arraylist for each branch.
     **/
	private ArrayList<String>[] pathState;
    
	/** Store path time as an array of arraylist of double.
	 * There are branchNum of such arraylist for each branch 
     **/
	private ArrayList<Double>[] pathTime;
    
	/**
	 * Constructor
	 * @param trees: 3 trees for 3 sites
	 */
	public TripletParser(Tree[] trees) {
		this.trees = new Tree[3];
		for (int i = 0; i < 3; i++) {
			this.trees[i] = trees[i];
		}
		
		numberOfChanges = new int[4][trees[0].getNumBranches()];
		timeOfStates = new double[3][trees[0].getNumBranches()];
		branchLengths = new double[trees[0].getNumBranches()];
		propStates = new double[3][trees[0].getNumBranches()];
		changesInGroups = new int[9][trees[0].getNumBranches()];
	
		pathState = (ArrayList<String>[])new ArrayList[trees[0].getNumBranches()];
		pathTime = (ArrayList<Double>[])new ArrayList[trees[0].getNumBranches()]; 

		for(int i = 0; i < trees[0].getNumBranches(); i++) {
			pathState[i] = new ArrayList<String>();
			pathTime[i] = new ArrayList<Double>();
		}
		
		constructWholePath();
	
		setBranchLengths(computeBranchLength());
		setNumberOfChanges(computeNumberOfChanges());
		setTimeOfStates(computeTimeOfStates());
		setPropStates(computePropState());
	}
    
    
	/** Return the number of branches in the tree
	 * @return number of branches
	 */
	public int numberOfBranches() {
		return branchLengths.length;
	}

	
	/**
	 * Wrapper class for constructPath(). Construct path for all branches.
	 */
	public void constructWholePath() {
		for (int i = 0; i < numberOfBranches(); i++) {
			constructPath(i);
		}
	}
	
	/**
	 * Construct the path state and path time by combining three sites/trees for a given branch.
	 * path state: triplet starting from the parent node to current node
	 * path time: time interval between substitution events.
	 *            The last time will be branch length of site 2 minus the time at the last substitution event.
	 * @param branchIndex index of the given branch (the index of the node that ends the branch)
	 */
	public void constructPath(int branchIndex) {
		ArrayList<String>[] state = (ArrayList<String>[])new ArrayList[3];
		ArrayList<Double>[] time = (ArrayList<Double>[])new ArrayList[3];
		HashMap<Double, Integer> map = new HashMap<Double, Integer>();
		
		//create a HashMap: key = time when substitutions occurred, value = site where substitutions occurred
		//do not add the last time in the HashMap b/c they will be the same for all sites.
		//i.e. last time = branch length
		for (int i = 0; i < 3; i ++) {
			state[i] = trees[i].getNodeByNodeNum(branchIndex).getPathState();
			time[i] = cumulativeSum(trees[i].getNodeByNodeNum(branchIndex).getPathTime());
			for (int j = 0; j < time[i].size()-1; j ++) {
				map.put(time[i].get(j), i);
			}
		}
		
		//starting state (triplet)
		int[] index = new int[3]; //track the index of states for three sites
		pathState[branchIndex].add(state[0].get(index[0])+state[1].get(index[1])+state[2].get(index[2]));
		
		//sort the HashMap by key => create a TreeMap (ordered) from HashMap
        Map<Double, Integer> treemap = new TreeMap<Double, Integer>(map);
        
        //iterate the TreeMap
        Set set = treemap.entrySet();
        Iterator iterator = set.iterator();
        double previous = 0.0;
        while(iterator.hasNext()) {
             Map.Entry me = (Map.Entry)iterator.next();
             //System.out.print(me.getKey() + ": ");
             //System.out.println(me.getValue());
             pathTime[branchIndex].add(((Double) me.getKey() - previous)); //store time interval length but not cumulative time
             previous = (Double) me.getKey();
             int position = (Integer) me.getValue();
             index[position] ++; //increment the index of pathState[position]
             char[] lastState = pathState[branchIndex].get(pathState[branchIndex].size()-1).toCharArray();
             lastState[position] = state[position].get(index[position]).charAt(0); //replace character at position
             pathState[branchIndex].add(String.valueOf(lastState));             
        }
        //add last time interval (branch length - last event time)
        pathTime[branchIndex].add(trees[1].getBranchLength()[branchIndex]-previous);
        //add last state. i.e. state at current node
        pathState[branchIndex].add(state[0].get(state[0].size()-1)+state[1].get(state[1].size()-1)+state[2].get(state[2].size()-1));        
		
	}
	
	/**
	 * Return pathState for a given branch
	 * @param branchIndex
	 * @return
	 */
	public ArrayList<String> getPathStateAtBranch(int branchIndex) {
		//constructPath(branchIndex);
		return pathState[branchIndex];
	}
	
	/**
	 * Return pathTime for a given branch
	 * @param branchIndex
	 * @return
	 */
	public ArrayList<Double> getPathTimeAtBranch(int branchIndex) {
		//constructPath(branchIndex);
		return pathTime[branchIndex];
	} 
	
	/**
	 * Return an arraylist which contains the cumulative sum of the given arraylist.
	 * e.g. given {1, 2, 3, 4, 5} would return {1, 3, 6, 10, 15}
	 * @param time
	 * @return cumulative sum
	 */
	public ArrayList<Double> cumulativeSum(ArrayList<Double> time) {
		ArrayList<Double> list = new ArrayList<Double>();
		int length = time.size();
		double sum = 0;
		for (int i = 0; i < length; i ++) {
			sum += time.get(i);
			list.add(sum);
		}
		return list;
	}
	
	
	/**
	 * Assign number of changes for each type of changes on each branch.
	 * There are 4 types of changes: CpG / non-CpG transition/transversion
	 * @return numberOfChanges
	 */
	private int[][] computeNumberOfChanges() {
		int[][] changes = new int[4][numberOfBranches()];
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < numberOfBranches(); j++) {
				changes[i][j] = computeNumberOfChangesAtBranch(j, i);
			}
		}		
		return changes;
	}
	
	/**
	 * Count the number of given type of changes on a given branch
	 * type 0: non-CpG transversion
	 * type 1: non-CpG transition
	 * type 2: CpG transversion
	 * type 3: CpG transition
	 * TRANSVERSION_ORDER = {"ATC","TAG","CAG","GCT"}
	 * TRANSITION_ORDER = {"AG","TC","CT","GA"}
	 * TRANSVERSION_GROUP = {"GCCG","GTCA","TAAT","TGAC"};
	 * TRANSITION_GROUP = {"GACT","AGTC"}
	 * changesInGroups[9][branchNum]: 0-3 for non-CpG transversion, 4-5 for non-CpG transition
	 *                                6-7 for CpG transversion, 8 for CpG transition
	 * @param branchIndex
	 * @param type 
	 * @return number of given type of changes on a given branch
	 */
	private int computeNumberOfChangesAtBranch(int branchIndex, int type) {
		int count = 0;
		switch(type) {
			case 0:
				for (int i = 0; i < pathState[branchIndex].size() - 1; i++) {
					if(pathState[branchIndex].get(i).charAt(1) != pathState[branchIndex].get(i+1).charAt(1)) {
						if(!pathState[branchIndex].get(i).contains("CG")) { //non-CpG site
							if(pathState[branchIndex].get(i).charAt(1) == TRANSVERSION_ORDER[0].charAt(0) && (pathState[branchIndex].get(i+1).charAt(1) == TRANSVERSION_ORDER[0].charAt(1) || pathState[branchIndex].get(i+1).charAt(1) == TRANSVERSION_ORDER[0].charAt(2))) {
								count++;
								//start with `A'
								if(pathState[branchIndex].get(i+1).charAt(1) == TRANSVERSION_GROUP[2].charAt(3)) {
									changesInGroups[2][branchIndex]++; //A->T
								} else if (pathState[branchIndex].get(i+1).charAt(1) == TRANSVERSION_GROUP[3].charAt(3)) {
									changesInGroups[3][branchIndex]++; //A->C
								}
							} else if(pathState[branchIndex].get(i).charAt(1) == TRANSVERSION_ORDER[1].charAt(0) && (pathState[branchIndex].get(i+1).charAt(1) == TRANSVERSION_ORDER[1].charAt(1) || pathState[branchIndex].get(i+1).charAt(1) == TRANSVERSION_ORDER[1].charAt(2))) {
								count++;
								//start with `T'
								if(pathState[branchIndex].get(i+1).charAt(1) == TRANSVERSION_GROUP[2].charAt(1)) {
									changesInGroups[2][branchIndex]++; //T->A
								} else if (pathState[branchIndex].get(i+1).charAt(1) == TRANSVERSION_GROUP[3].charAt(1)) {
									changesInGroups[3][branchIndex]++; //T->G
								}
							} else if(pathState[branchIndex].get(i).charAt(1) == TRANSVERSION_ORDER[2].charAt(0) && (pathState[branchIndex].get(i+1).charAt(1) == TRANSVERSION_ORDER[2].charAt(1) || pathState[branchIndex].get(i+1).charAt(1) == TRANSVERSION_ORDER[2].charAt(2))) {
								count++;
								//start with `C'
								if(pathState[branchIndex].get(i+1).charAt(1) == TRANSVERSION_GROUP[0].charAt(3)) {
									changesInGroups[0][branchIndex]++; //C->G
								} else if (pathState[branchIndex].get(i+1).charAt(1) == TRANSVERSION_GROUP[1].charAt(3)) {
									changesInGroups[1][branchIndex]++; //C->A
								}
							} else if(pathState[branchIndex].get(i).charAt(1) == TRANSVERSION_ORDER[3].charAt(0) && (pathState[branchIndex].get(i+1).charAt(1) == TRANSVERSION_ORDER[3].charAt(1) || pathState[branchIndex].get(i+1).charAt(1) == TRANSVERSION_ORDER[3].charAt(2))) {
								count++;
								//start with `G'
								if(pathState[branchIndex].get(i+1).charAt(1) == TRANSVERSION_GROUP[0].charAt(1)) {
									changesInGroups[0][branchIndex]++; //G->C
								} else if (pathState[branchIndex].get(i+1).charAt(1) == TRANSVERSION_GROUP[1].charAt(1)) {
									changesInGroups[1][branchIndex]++; //G->T
								}
							}
						}
					}
				}
				break;
			case 1:
				for (int i = 0; i < pathState[branchIndex].size() - 1; i++) {
					if(pathState[branchIndex].get(i).charAt(1) != pathState[branchIndex].get(i+1).charAt(1)) {
						if(!pathState[branchIndex].get(i).contains("CG")) { //non-CpG site
							if(pathState[branchIndex].get(i).charAt(1) == TRANSITION_ORDER[0].charAt(0) && pathState[branchIndex].get(i+1).charAt(1) == TRANSITION_ORDER[0].charAt(1)) {
								count++;
								//A->G
								changesInGroups[5][branchIndex]++;
							} else if(pathState[branchIndex].get(i).charAt(1) == TRANSITION_ORDER[1].charAt(0) && pathState[branchIndex].get(i+1).charAt(1) == TRANSITION_ORDER[1].charAt(1)) {
								count++;
								//T->C
								changesInGroups[5][branchIndex]++;
							} else if(pathState[branchIndex].get(i).charAt(1) == TRANSITION_ORDER[2].charAt(0) && pathState[branchIndex].get(i+1).charAt(1) == TRANSITION_ORDER[2].charAt(1)) {
								count++;
								//C->T
								changesInGroups[4][branchIndex]++;
							} else if(pathState[branchIndex].get(i).charAt(1) == TRANSITION_ORDER[3].charAt(0) && pathState[branchIndex].get(i+1).charAt(1) == TRANSITION_ORDER[3].charAt(1)) {
								count++;
								//G->A
								changesInGroups[4][branchIndex]++;
							}
						}
					}
				}
				break;
			case 2:
				for (int i = 0; i < pathState[branchIndex].size() - 1; i++) {
					if(pathState[branchIndex].get(i).charAt(1) != pathState[branchIndex].get(i+1).charAt(1)) {
						if(pathState[branchIndex].get(i).contains("CG")) { //CpG site
							if(pathState[branchIndex].get(i).charAt(1) == TRANSVERSION_ORDER[2].charAt(0) && (pathState[branchIndex].get(i+1).charAt(1) == TRANSVERSION_ORDER[2].charAt(1) || pathState[branchIndex].get(i+1).charAt(1) == TRANSVERSION_ORDER[2].charAt(2))) {
								count++;
								//start with `C'
								if(pathState[branchIndex].get(i+1).charAt(1) == TRANSVERSION_GROUP[0].charAt(3)) {
									changesInGroups[6][branchIndex]++; //C->G
								} else if (pathState[branchIndex].get(i+1).charAt(1) == TRANSVERSION_GROUP[1].charAt(3)) {
									changesInGroups[7][branchIndex]++; //C->A
								}
							} else if(pathState[branchIndex].get(i).charAt(1) == TRANSVERSION_ORDER[3].charAt(0) && (pathState[branchIndex].get(i+1).charAt(1) == TRANSVERSION_ORDER[3].charAt(1) || pathState[branchIndex].get(i+1).charAt(1) == TRANSVERSION_ORDER[3].charAt(2))) {
								count++;
								//start with 'G'
								if(pathState[branchIndex].get(i+1).charAt(1) == TRANSVERSION_GROUP[0].charAt(1)) {
									changesInGroups[6][branchIndex]++; //G->C
								} else if (pathState[branchIndex].get(i+1).charAt(1) == TRANSVERSION_GROUP[1].charAt(1)) {
									changesInGroups[7][branchIndex]++; //G->T
								}
							} 
						}
					}
				}
				break;
			case 3:
				for (int i = 0; i < pathState[branchIndex].size() - 1; i++) {
					if(pathState[branchIndex].get(i).charAt(1) != pathState[branchIndex].get(i+1).charAt(1)) {
						if(pathState[branchIndex].get(i).contains("CG")) { //CpG site
							if(pathState[branchIndex].get(i).charAt(1) == TRANSITION_ORDER[2].charAt(0) && pathState[branchIndex].get(i+1).charAt(1) == TRANSITION_ORDER[2].charAt(1)) {
								count++;
								//C->T
								changesInGroups[8][branchIndex]++;
							} else if(pathState[branchIndex].get(i).charAt(1) == TRANSITION_ORDER[3].charAt(0) && pathState[branchIndex].get(i+1).charAt(1) == TRANSITION_ORDER[3].charAt(1)) {
								count++;
								//G->A
								changesInGroups[8][branchIndex]++;
							}
						}
					}
				}
				break;
		}
		
		
		return count;
	}
	
	/**
	 * Assign times of each state on each branch.
	 * There are three states: non-CpG C+G, non-CpG A+T, CpG C+G
	 * @return timeOfStates
	 */
	private double[][] computeTimeOfStates() {
		double[][] times = new double[3][numberOfBranches()];
		for (int j = 0; j < numberOfBranches(); j++) {
			times[0][j] = computeTimeOfnonCpGCGAtBranch(j); //non-CpG C+G
			times[1][j]	= branchLengths[j]-computeTimeOfCpGAtBranch(j) - computeTimeOfnonCpGCGAtBranch(j); //non-CpG A+T
			times[2][j] = computeTimeOfCpGAtBranch(j); //CpG site
		}				
		return times;
	}
	
	/**
	 * Compute time of CpG site at the second position on a given branch
	 * @param branchIndex
	 * @return
	 */
	private double computeTimeOfCpGAtBranch(int branchIndex) {
		double timeState = 0;
		
		for (int i = 0; i < pathState[branchIndex].size()-1; i++) {
			if(pathState[branchIndex].get(i).contains("CG")) {
				timeState += pathTime[branchIndex].get(i);
			}
		}
		
		return timeState;
	}
	
	/**
	 * Compute time of non-CpG C+G site at the second position on a given branch
	 * @param branchIndex
	 * @return
	 */
	private double computeTimeOfnonCpGCGAtBranch(int branchIndex) {
		double timeState = 0;
		
		for (int i = 0; i < pathState[branchIndex].size()-1; i++) {
			if(!pathState[branchIndex].get(i).contains("CG")) {
				if((pathState[branchIndex].get(i).charAt(1) == 'C') || (pathState[branchIndex].get(i).charAt(1) == 'G')) {
					timeState += pathTime[branchIndex].get(i);
				}
			}
		}
		
		return timeState;
	}
	
	/**
	 * Assign branch lengths for each branch. 
	 * Branch lengths are supposed to be the same for all three sites, but in case of rounding error, I pick site 2.
	 * @return branchLengths
	 */
	private double[] computeBranchLength() {
		double[] br = new double[numberOfBranches()];
		for (int i = 0; i < trees[1].getNumBranches(); i++) {
			br[i] = trees[1].getBranchLength()[i];
		}
		return br;
	}
	
	/**
	 * Compute the proportion of time that the second position is a non-CpG C+G / non-CpG A+T / CpG site.
	 * @return prop
	 */
	private double[][] computePropState() {
		double[][] prop = new double[3][numberOfBranches()];
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < numberOfBranches(); j++) {
				prop[i][j] = timeOfStates[i][j]/branchLengths[j];
			}
		}		
		return prop;
	}
	
	/**
	 * Get proportion of time the second position is a CpG / non-CpG site for each branch. 
	 * dim = 3 x (# of branches)
	 * @return propStates
	 */
	public double[][] getPropStates() {
		return propStates;
	}
	
    /**
     * Set proportion of states
     * @param propStates to be set
     */
	public void setPropStates(double[][] propStates) {
		this.propStates = propStates;
	}

    /**
     * Return the number of changes for each substitution type (scenario 1) on branches
     * @return numberOfChanges
     */
	public int[][] getNumberOfChanges() {
		return numberOfChanges;
	}
	
    /**
     * Return the number of changes for each substitution type (scenario 2) on branches
     * @return numberOfChanges
     */
	public int[][] getChangesInGroups() {
		return changesInGroups;
	}
	
    /**
     * Set the number of changes for each substitution type (scenario 1) on branches
     * @param numberOfChanges
     */
	private void setNumberOfChanges(int[][] numberOfChanges) {
		this.numberOfChanges = numberOfChanges;
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

}
