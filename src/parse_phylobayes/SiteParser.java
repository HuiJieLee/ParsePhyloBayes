package parse_phylobayes;

/**
 * Interface for SingleSiteParse and TripletParser
 * @author Hui-Jie Lee
 *
 */
public interface SiteParser {
    
	/** Return number of branches in the tree **/
	public int numberOfBranches();
	
    /** Return the proportion of time in states on branches **/
	public double[][] getPropStates();
	
    /** Return the number of changes on branches **/
	public int[][] getNumberOfChanges();
	
    /** Return the total time in states on branches **/
	public double[][] getTimeOfStates();
	
    /** Return the branch lengths for each branch **/
	public double[] getBranchLengths();
	
    /** Return the number of changes in each substitution group/type **/
	public int[][] getChangesInGroups();
	
	
	
	
}
