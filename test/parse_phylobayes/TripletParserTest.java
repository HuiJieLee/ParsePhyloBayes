/**
 * 
 */
package parse_phylobayes;

import java.io.StringReader;
import java.util.ArrayList;
import java.util.Arrays;

import junit.framework.TestCase;

/**
 * @author Hui-Jie Lee
 *
 */
public class TripletParserTest extends TestCase {
	
	/** A TreeParser to construct a tree */
	private TreeParser parser;
	/** An array of Tree used to perform tests */
	private Tree[] trees;
	/** A TripletParser to perform tests */
	private TripletParser sites;

	/* (non-Javadoc)
	 * @see junit.framework.TestCase#setUp()
	 */
	protected void setUp() throws Exception {
		trees = new Tree[3];
		
		//from gtr_322.map gtr_323.map gtr_324.map
		//outgroup has not been removed yet. handle outgroup in MappingParser.java
		//unrooted trees
		StringReader str = new StringReader("(((human_A:0.101652:A,orangutan_A:0.175401:A)_A:0.649724:A,(rhesus_G:0.00557255:G:0.314872:A,baboon_G:0.136734:G:0.279524:A)_A:0.27683:A)_A:1.42672:A:0.197287:T:0.302921:A:0.677233:T:0.683041:C:0.366039:G:0.300407:A:1.29285:T,marmoset_C:0.0890465:C:1.71246:A:0.0718244:G:0.685267:A:2.21303:T,bushbaby_G:0.146267:G:1.23497:T:0.433369:C:0.174881:G:1.14392:T:0.192338:C:0.350913:G:0.364961:A:0.508608:C:1.1939:G:1.20611:T:0.687221:G:0.548099:T)_T;");
        parser = new TreeParser(str,"");
        trees[0] = parser.tokenize();
        
        str = new StringReader("(((human_G:0.101652:G,orangutan_G:0.175401:G)_G:0.649724:G,(rhesus_G:0.320445:G,baboon_G:0.416258:G)_G:0.27683:G)_G:5.24651:G,marmoset_C:3.52904:C:0.295141:T:0.947449:G,bushbaby_G:8.18555:G)_G;");
        parser = new TreeParser(str,"");
        trees[1] = parser.tokenize();
        
        str = new StringReader("(((human_G:0.101652:G,orangutan_G:0.175401:G)_G:0.649724:G,(rhesus_G:0.320445:G,baboon_G:0.416258:G)_G:0.27683:G)_G:5.24651:G,marmoset_T:2.89905:T:0.630191:C:1.24239:G,bushbaby_G:8.18555:G)_G;");
        parser = new TreeParser(str,"");
        trees[2] = parser.tokenize();
        
        sites = new TripletParser(trees);
	}

	public void testNumberOfBranches() {
		assertEquals(9, sites.numberOfBranches());
	}
	
	public void testGetPathStateAtBranch() {		
		/**
		 * site1: marmoset_C:0.0890465:C:1.71246:A:0.0718244:G:0.685267:A:2.21303:T
		 * {2.21303, 2.898297, 2.970121, 4.682581}
		 * site2: marmoset_C:3.52904:C:0.295141:T:0.947449:G
		 * {0.947449, 1.24259}
         * site3: marmoset_T:2.89905:T:0.630191:C:1.24239:G
         * {1.24239, 1.872581}
		 */
	
		ArrayList<String> state = new ArrayList<String>();
		state.add("TGG"); //start at parent node
		state.add("TTG"); //t = 0.947449, site 2 from G -> T
		state.add("TTC"); //t = 1.24239, site 3 from G -> C 
		state.add("TCC"); //t = 1.24259, site 2 from T -> C
		state.add("TCT"); //t = 1.872581, site 3 from C -> T
		state.add("ACT"); // t = 2.21303, site 1 from T -> A
		state.add("GCT"); // t = 2.898297, site 1 from A -> G
		state.add("ACT"); // t = 2.970121, site 1 from G -> A
		state.add("CCT"); // t= 4.692581, site 1 from A -> C
		state.add("CCT"); // site at current node
		
		//test marmoset branch, index = 4
		assertEquals(state, sites.getPathStateAtBranch(4));		
	}

	public void testGetPathTimeAtBranch(){
		/**
		 * site1: marmoset_C:0.0890465:C:1.71246:A:0.0718244:G:0.685267:A:2.21303:T
		 * {2.21303, 2.898297, 2.970121, 4.682581}
		 * site2: marmoset_C:3.52904:C:0.295141:T:0.947449:G
		 * {0.947449, 1.24259}
         * site3: marmoset_T:2.89905:T:0.630191:C:1.24239:G
         * {1.24239, 1.872581}
		 */
		ArrayList<Double> time = new ArrayList<Double>();
		time.add(0.947449);
		time.add(0.294941);
		time.add(0.0002);
		time.add(0.629991);
		time.add(0.340449);
		time.add(0.685267);
		time.add(0.071824);
		time.add(1.71246);
		time.add(0.0890465);
		
		//test marmoset branch, index = 4
		assertEquals(time.size(),sites.getPathTimeAtBranch(4).size());
		for (int i = 0; i < time.size(); i++) {
			assertTrue(time.get(i)-sites.getPathTimeAtBranch(4).get(i)<0.0001);
		}
		
	}
	
	public void testGetBranchLengths() {
		double[] br = new double[] {0.101652,0.175401,0.320445,0.416258,4.77163,8.18555,0.27683,0.649724,5.24651};
		assertTrue(Arrays.equals(br, sites.getBranchLengths()));
	}
	
	public void testGetNumberOfChanges() {
		int[][] change = new int[][] {
				  {0,0,0,0,1,0,0,0,0},
				  {0,0,0,0,1,0,0,0,0},
				  {0,0,0,0,0,0,0,0,0},
				  {0,0,0,0,0,0,0,0,0}};
		
		for(int i = 0; i < 4; i++) {
			assertTrue(Arrays.equals(change[i],sites.getNumberOfChanges()[i]));
		}
	}
	
	public void testGetTimeOfStates() {
		double[][] time = new double[][] {
				  {0.101652,0.175401,0.320445,0.416258,4.77163,7.051235,0.27683,0.649724,4.563469},
				  {0,0,0,0,0,1.134315,0,0,0.683041}};
		for (int i = 0; i < 2; i++) {
			//System.out.println(i);
			for (int j = 0; j < 9; j++) {
				//System.out.println(j);
				//System.out.println(sites.getTimeOfStates()[i][j]);
				assertTrue(time[i][j] - sites.getTimeOfStates()[i][j]<0.0001);
			}
		}
	}

	public void testGetPropStates() {
		double[][] prop = new double[][] {
				{1,1,1,1,1,7.051235/8.18555,1,1,4.563469/5.24651},
				{0,0,0,0,0,1.134315/8.18555,0,0,0.683041/5.24651}
		};
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 9; j++) {
				assertTrue(prop[i][j] - sites.getPropStates()[i][j]<0.0001);
			}
		}
	}
	
	/* (non-Javadoc)
	 * @see junit.framework.TestCase#tearDown()
	 */
	protected void tearDown() throws Exception {
		super.tearDown();
	}

}
