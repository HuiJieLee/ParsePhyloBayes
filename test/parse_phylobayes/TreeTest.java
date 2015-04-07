/**
 * 
 */
package parse_phylobayes;

import java.io.StringReader;
import java.util.ArrayList;
import java.util.Arrays;

import junit.framework.TestCase;

/**
 * Test class for Tree. Tests that the object is constructed correctly
 * and that the state of the object is manipulated properly.
 * @author Hui-Jie Lee
 *
 */
public class TreeTest extends TestCase {

	/** A TreeParser to construct a tree */
	private TreeParser parser;
	/** A Tree used to perform tests */
	private Tree tree;
	
	/** A TreeParser to construct a tree */
	private TreeParser parser1;
	/** A Tree used to perform tests */
	private Tree tree1;
	
	/**
	 * Sets up the TreeTest by creating one Tree objects. 
	 */
	protected void setUp() throws Exception {
		StringReader str = new StringReader("(((A_A:0.1:C:0.5:C,B_A:0.4:C:0.3:C)_C:0.7:G:0.3:G,(C_A:0.3:A:0.1:T:0.1:T,D_T:0.5:T)_T:0.8:T:0.2:G:0.3:G)_G:0.8:G,((E_C:0.4:C,F_C:0.7:C)_C:0.7:C,G_A:0.4:C:0.6:C)_C:0.9:G:0.3:G)_G;");
		StringReader str1 = new StringReader("(((human_A:0.100993:A,orangutan_A:0.0694486:A)_A:0.198984:A,(rhesus_A:0.195176:A,baboon_A:0.233133:A)_A:0.171989:A)_A:0.406658:A,marmoset_A:0.297841:A,bushbaby_G:3.73678e-05:G:0.490021:A)_A;");
        
		parser = new TreeParser(str,"");
        tree = parser.tokenize();
        
        parser1 = new TreeParser(str1,"");
        tree1 = parser1.tokenize();
	}
	
	public void testNumBranches() {
		assertEquals(12, tree.getNumBranches());
		//assertEquals(10, tree1.getNumBranches());
	}
	
	public void testGetPathState(){
		ArrayList<String> states = new ArrayList<String>();
		states.add("C");
		states.add("C");
		states.add("A");
		assertEquals(states, tree.getNodeByNodeNum(0).getPathState());
		assertEquals(states, tree.getNodeByNodeNum(1).getPathState());
		
		states.clear();
		states.add("T");
		states.add("T");
		states.add("A");
		states.add("A");
		//System.out.print(tree.getNodeByNodeNum(2).parent().getState());
		//System.out.print(tree.getNodeByNodeNum(2).getPathState().get(0));
		//System.out.print(tree.getNodeByNodeNum(2).getPathState().get(1));
		//System.out.print(tree.getNodeByNodeNum(2).getState());
		//System.out.print(states.get(0));
		//System.out.print(states.get(1));
		//System.out.print(states.get(2));
		//System.out.print(states.get(3));
		//System.out.print(states.size());
		assertEquals(states, tree.getNodeByNodeNum(2).getPathState());	
	}
	
	public void testGetPathTime(){
		ArrayList<Double> times = new ArrayList<Double>();
		times.add(0.5);
		times.add(0.1);
		assertEquals(times, tree.getNodeByNodeNum(0).getPathTime());
		
		times.clear();
		times.add(0.1);
		times.add(0.1);
		times.add(0.3);
		assertEquals(times, tree.getNodeByNodeNum(2).getPathTime());
		
	}
	
	
	public void testGetNumberOfChanges() {
		int[] change = new int[] {0,0,0,0,0,0,0,0,0,0,0,0};
		assertTrue(Arrays.equals(change, tree.getNumberOfChanges("A","G")));
		assertTrue(Arrays.equals(change, tree.getNumberOfChanges("A","C")));
		assertTrue(Arrays.equals(change, tree.getNumberOfChanges("A","T")));
		assertTrue(Arrays.equals(change, tree.getNumberOfChanges("G","A")));
		change[8] = 1;
		change[10] = 1;
		assertTrue(Arrays.equals(change, tree.getNumberOfChanges("G","C")));
		change[8] = 0;
		change[9] = 1;
		change[10] = 0;
		assertTrue(Arrays.equals(change, tree.getNumberOfChanges("G","T")));
		change[0] = 1;
		change[1] = 1;
		change[6] = 1;
		change[9] = 0;
		assertTrue(Arrays.equals(change, tree.getNumberOfChanges("C","A")));
		change[0] = 0;
		change[1] = 0;
		change[6] = 0;
		assertTrue(Arrays.equals(change, tree.getNumberOfChanges("C","G")));
		assertTrue(Arrays.equals(change, tree.getNumberOfChanges("C","T")));
		change[2] = 1;
		assertTrue(Arrays.equals(change, tree.getNumberOfChanges("T","A")));
		change[2] = 0;
		assertTrue(Arrays.equals(change, tree.getNumberOfChanges("T","G")));
		assertTrue(Arrays.equals(change, tree.getNumberOfChanges("T","C")));			
	}
	
	public void testTimeOfState() {
		double[] time = new double[] {0,0,0.3,0,0,0,0,0,0,0,0,0}; 
		assertTrue(Arrays.equals(time, tree.getTimeOfState("A")));
		time[2] = 0.2;
		time[3] = 0.5;
		time[9] = 0.8;
		assertTrue(Arrays.equals(time, tree.getTimeOfState("T")));
		time[0] = 0.6;
		time[1] = 0.7;
		time[2] = 0;
		time[3] = 0;
		time[4] = 0.4;
		time[5] = 0.7;
		time[6] = 1;
		time[7] = 0.7;
		time[9] = 0;
		assertTrue(Arrays.equals(time, tree.getTimeOfState("C")));
		time = new double[] {0,0,0,0,0,0,0,0,1.2,0.5,1,0.8}; 
		assertTrue(Arrays.equals(time, tree.getTimeOfState("G")));		
	}
	
	public void testGetBranchLength() {
		double[] br = new double[] {0.6,0.7,0.5,0.5,0.4,0.7,1,0.7,1.2,1.3,1,0.8};
		assertTrue(Arrays.equals(br, tree.getBranchLength()));
	}
	
	public void testGetPropOfState() {
		double[] prop = new double[] {0,0,0.3/0.5,0,0,0,0,0,0,0,0,0};
		assertTrue(Arrays.equals(prop, tree.getPropState("A")));
		prop[2] = 0.2/0.5;
		prop[3] = 1;
		prop[9] = 0.8/1.3;
		assertTrue(Arrays.equals(prop, tree.getPropState("T")));
		prop = new double[] {1,1,0,0,1,1,1,1,0,0,0,0};
		assertTrue(Arrays.equals(prop, tree.getPropState("C")));
		prop = new double[] {0,0,0,0,0,0,0,0,1,0.5/1.3,1,1};
		assertTrue(Arrays.equals(prop, tree.getPropState("G")));
	}
	
	/**
     * Tears down the test fixture.
     * (Called after every test case method.)
     */
    protected void tearDown() {
         parser = null;
         tree = null;
    }

}
