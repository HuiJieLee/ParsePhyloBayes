/**
 * 
 */
package parse_phylobayes;

import java.io.StringReader;
import java.util.Arrays;

import junit.framework.TestCase;

/**
 * @author Hui-Jie Lee
 *
 */
public class SingleSiteParserTest extends TestCase {

	/** A TreeParser to construct a tree */
	private TreeParser parser;
	/** A Tree used to perform tests */
	private Tree tree;
	/** A SingleSiteParser to perform tests */
	private SingleSiteParser site;

	/**
	 * Sets up the SingleSiteParserTest by creating one SingleSiteParser objects. 
	 */
	protected void setUp() throws Exception {
		StringReader str = new StringReader("(((A_A:0.1:C:0.5:C,B_A:0.4:C:0.3:C)_C:0.7:G:0.3:G,(C_A:0.3:A:0.1:T:0.1:T,D_T:0.5:T)_T:0.8:T:0.2:G:0.3:G)_G:0.8:G,((E_C:0.4:C,F_C:0.7:C)_C:0.7:C,G_A:0.4:C:0.6:C)_C:0.9:G:0.3:G)_G;");
        parser = new TreeParser(str,"");
        tree = parser.tokenize();
        site = new SingleSiteParser(tree);
	}
	
	public void testNumberOfBranches() {
		assertEquals(12, site.numberOfBranches());
	}
	
	public void testGetBranchLengths() {
		double[] br = new double[] {0.6,0.7,0.5,0.5,0.4,0.7,1,0.7,1.2,1.3,1,0.8};
		assertTrue(Arrays.equals(br, site.getBranchLengths()));
	}
	
	public void testGetNumberOfChanges() {
		int[][] change = new int[][] {{0,0,0,0,0,0,0,0,0,0,0,0},
									  {0,0,0,0,0,0,0,0,0,0,0,0},
									  {0,0,0,0,0,0,0,0,0,0,0,0},
									  {0,0,0,0,0,0,0,0,0,0,0,0},
									  {0,0,0,0,0,0,0,0,1,0,1,0},
									  {0,0,0,0,0,0,0,0,0,1,0,0},
									  {1,1,0,0,0,0,1,0,0,0,0,0},
									  {0,0,0,0,0,0,0,0,0,0,0,0},
									  {0,0,0,0,0,0,0,0,0,0,0,0},
									  {0,0,1,0,0,0,0,0,0,0,0,0},
									  {0,0,0,0,0,0,0,0,0,0,0,0},
									  {0,0,0,0,0,0,0,0,0,0,0,0}};
		
		for(int i = 0; i < 12; i++) {
			assertTrue(Arrays.equals(change[i],site.getNumberOfChanges()[i]));
		}
	}
		
	public void testGetTimeOfStates() {
		double[][] time = new double[][] {{0,0,0.3,0,0,0,0,0,0,0,0,0},
										  {0,0,0,0,0,0,0,0,1.2,0.5,1,0.8},
										  {0.6,0.7,0,0,0.4,0.7,1,0.7,0,0,0,0},
										  {0,0,0.2,0.5,0,0,0,0,0,0.8,0,0}};
		for (int i = 0; i < 4; i++) {
			assertTrue(Arrays.equals(time[i],site.getTimeOfStates()[i]));
		}
	}
	
	public void testGetPropStates() {
		double[][] prop = new double[][] {{0,0,0.3/0.5,0,0,0,0,0,0,0,0,0},
										  {0,0,0,0,0,0,0,0,1,0.5/1.3,1,1},
										  {1,1,0,0,1,1,1,1,0,0,0,0},
										  {0,0,0.2/0.5,1,0,0,0,0,0,0.8/1.3,0,0}};
		for (int i = 0; i < 4; i++) {
			assertTrue(Arrays.equals(prop[i],site.getPropStates()[i]));
		}
	}
	
	/**
     * Tears down the test fixture.
     * (Called after every test case method.)
     */
    protected void tearDown() {
         parser = null;
         tree = null;
         site = null;
    }

}
