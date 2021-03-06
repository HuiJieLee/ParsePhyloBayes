/**
 * 
 */
package parse_phylobayes;

import java.util.Arrays;

import junit.framework.TestCase;

/**
 * @author Hui-Jie Lee
 *
 */
public class MappingParserTest extends TestCase {

	/** A MappingParser to perform test */
	private MappingParser parse;
	/** A MappingParser to perform test */
	private MappingParser parse1;

	/**
	 * Sets up the MappingParserTest by creating one MappingParser objects. 
	 */
	protected void setUp() throws Exception {
		parse = new MappingParser(4, 2, "test",1, "");
		parse1 = new MappingParser(4, 2, "neutral",2, "");
	}
	
	public void testGetBranchNum() {
		assertEquals(12, parse.getBranchNum());
		assertEquals(9, parse1.getBranchNum());
	}
	
	public void testGetNumberOfChanges() {
/*		int[][][] num = new int[][][] {
				{{1,0,0,0,0,0,1,0,0,0,0,0},
				 {1,0,0,0,0,0,1,0,0,0,0,0},
				 {0,0,0,0,0,0,0,0,1,1,0,0},
				 {0,0,0,1,0,0,0,0,0,0,0,0},
				 {0,0,0,0,0,0,0,0,0,0,0,0},
				 {0,0,0,0,0,0,0,0,0,0,1,0},
				 {1,0,0,0,0,0,1,0,0,0,1,0},
				 {0,0,0,0,0,0,0,0,0,0,1,0},
				 {0,1,0,0,1,0,1,0,1,0,0,0},
				 {1,0,0,0,0,1,0,0,0,0,0,1},
				 {0,0,0,0,1,0,0,0,0,0,0,0},
				 {0,1,0,0,0,0,1,0,0,0,0,0}					
				},
				{{1,1,0,0,0,0,2,0,0,0,0,0},
				 {1,0,0,0,0,0,1,0,0,0,0,0},
				 {0,0,0,0,0,0,0,0,1,1,0,0},
				 {0,0,0,1,0,0,0,0,0,0,0,0},
				 {0,0,0,0,0,0,0,0,0,0,0,0},
				 {0,0,0,0,0,0,0,0,0,0,1,0},
				 {1,0,0,0,0,0,1,0,0,0,1,0},
				 {0,0,0,0,0,0,0,0,0,0,1,0},
				 {0,1,0,0,1,0,1,0,1,0,0,0},
				 {1,0,0,0,1,1,0,1,0,0,0,1},
				 {0,0,0,0,1,0,0,0,0,0,0,0},
				 {0,1,0,0,0,0,1,0,0,0,0,0}					
				}
		};
		
		int[][][] num1 = new int[][][] {
				{{1,0,0,0,0,0,1,0,0,0,0,0},
				 {1,1,0,0,0,0,2,0,0,0,0,0}	
				},
				{{1,0,0,0,0,0,1,0,0,0,0,0},
			     {1,0,0,0,0,0,1,0,0,0,0,0}				
				},
				{{0,0,0,0,0,0,0,0,1,1,0,0},
				 {0,0,0,0,0,0,0,0,1,1,0,0}	
				},
				{{0,0,0,1,0,0,0,0,0,0,0,0},
				 {0,0,0,1,0,0,0,0,0,0,0,0} 	
				},
				{{0,0,0,0,0,0,0,0,0,0,0,0},
				 {0,0,0,0,0,0,0,0,0,0,0,0}
				},
				{{0,0,0,0,0,0,0,0,0,0,1,0},
				 {0,0,0,0,0,0,0,0,0,0,1,0}		
				},
				{{1,0,0,0,0,0,1,0,0,0,1,0},
				 {1,0,0,0,0,0,1,0,0,0,1,0}	
				},
				{{0,0,0,0,0,0,0,0,0,0,1,0},
				 {0,0,0,0,0,0,0,0,0,0,1,0}	
				},
				{{0,0,0,0,0,0,0,0,0,0,1,0},
				 {0,0,0,0,0,0,0,0,0,0,1,0}		
				},
				{{1,0,0,0,0,1,0,0,0,0,0,1},
				 {1,0,0,0,1,1,0,1,0,0,0,1}	
				},
				{{0,0,0,0,1,0,0,0,0,0,0,0},
				 {0,0,0,0,1,0,0,0,0,0,0,0}		
				},
				{{0,1,0,0,0,0,1,0,0,0,0,0},
				 {0,1,0,0,0,0,1,0,0,0,0,0}	
				}
		};
		int[][][] num2 = new int[][][] {
				{{1,1},{0,1},{0,0},{0,0},{0,0},{0,0},{1,2},{0,0},{0,0},{0,0},{0,0},{0,0}},
				{{1,1},{0,0},{0,0},{0,0},{0,0},{0,0},{1,1},{0,0},{0,0},{0,0},{0,0},{0,0}},
				{{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{1,1},{1,1},{0,0},{0,0}},
				{{0,0},{0,0},{0,0},{1,1},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}},
				{{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}},
				{{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{1,1},{0,0}},
				{{1,1},{0,0},{0,0},{0,0},{0,0},{0,0},{1,1},{0,0},{0,0},{0,0},{1,1},{0,0}},
				{{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{1,1},{0,0}},
				{{0,0},{1,1},{0,0},{0,0},{1,1},{0,0},{1,1},{0,0},{1,1},{0,0},{0,0},{0,0}},
				{{1,1},{0,0},{0,0},{0,0},{0,1},{1,1},{0,0},{0,1},{0,0},{0,0},{0,0},{1,1}},
				{{0,0},{0,0},{0,0},{0,0},{1,1},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}},
				{{0,0},{1,1},{0,0},{0,0},{0,0},{0,0},{1,1},{0,0},{0,0},{0,0},{0,0},{0,0}}
		};
*/		
		//[12][branchNum][C]
		//{1,1} = AG, branch0
		//{1,1} = AG, branch1
		int[][][] num3 = new int[][][]{
				{{1,1},{1,1},{0,0},{0,0},{0,0},{0,0},{1,1},{0,0},{0,0},{1,1},{0,0},{0,0}},
				{{0,1},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{1,1},{0,0},{0,0},{1,1}},
				{{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}},
				{{0,0},{0,0},{0,0},{1,1},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}},
				{{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{1,1},{0,1},{1,1},{0,0}},
				{{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{1,1},{0,0},{0,0}},
				{{1,2},{1,1},{0,0},{0,0},{0,0},{0,0},{1,1},{0,0},{1,1},{0,0},{0,0},{1,1}},
				{{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,1},{0,0},{0,0}},
				{{0,0},{0,0},{1,1},{0,0},{0,0},{0,0},{0,0},{0,0},{1,1},{0,0},{0,0},{0,0}},
				{{0,0},{0,0},{1,1},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}},
				{{0,0},{0,0},{0,0},{0,0},{0,0},{1,1},{1,1},{1,1},{0,0},{0,0},{0,0},{0,0}},
				{{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{1,1},{0,0},{0,0}}
		};
		for(int i = 0; i < 12; i++) {
			for(int j = 0; j < 12; j++){
				assertTrue(Arrays.equals(num3[i][j],parse.getNumberOfChanges()[i][j]));
			}
		}
		
		//numberOfChanges = new int[numTypeChanges][branchNum][C];
		num3 = new int[][][]{
				{{0,0},{0,0},{0,0},{0,0},{2,2},{0,2},{0,0},{0,0},{0,2}},
				{{0,0},{0,0},{0,0},{0,0},{2,1},{0,0},{0,0},{0,0},{0,0}},
				{{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}},
				{{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}}	
		};
		
		for(int i = 0; i < 4; i++) {
			//System.out.print("i="+i);
			for(int j = 0; j < 9; j++){
				//System.out.println("j="+j);
				for(int k = 0; k < 2; k++) {
					//System.out.print("k="+k+" ");
					//System.out.println(parse1.getNumberOfChanges()[i][j][k]);
					assertEquals(num3[i][j][k], parse1.getNumberOfChanges()[i][j][k]);
				}
				//assertTrue(Arrays.equals(num3[i][j],parse1.getNumberOfChanges()[i][j]));
			}
		}
	}
	
	
	public void testPropStates() {
		double[][][] prop = new double[][][] {
				{{0.4/0.6,0.5/0.6},{0.2/0.7,0.4/0.8},{0.3/0.5,0.3/0.5},{0.2/0.5,0.1/0.4},{0.4/0.4,0.5/0.5},{0.7/0.7,0.8/0.8},{0.5/1.0,0.6/1},{0.7/0.7,0.6/0.6},{1.5/1.2,1.3/1.1},{0.4/1.3,0.4/1.2},{1.0/1,1.0/1},{0.4/0.8,0.4/1}},
				{{0.2/0.6,0.3/0.6},{0.5/0.7,0.4/0.8},{0.5/0.5,0.5/0.5},{0.3/0.5,0.3/0.4},{0.4/0.4,0.5/0.5},{1.2/0.7,1.1/0.8},{1.3/1.0,0.9/1},{0.3/0.7,0.2/0.6},{1.2/1.2,1.1/1.1},{1.4/1.3,2.4/1.2},{1.0/1,1.0/1},{0.8/0.8,1.0/1}},
				{{1.2/0.6,1.0/0.6},{1.4/0.7,1.6/0.8},{0.7/0.5,0.8/0.5},{1.0/0.5,0.8/0.4},{0.4/0.4,0.5/0.5},{0.7/0.7,0.8/0.8},{1.0/1.0,1.0/1},{0.7/0.7,0.6/0.6},{0.9/1.2,0.9/1.1},{2.1/1.3,1.0/1.2},{1.0/1,1.0/1},{1.2/0.8,1.6/1}},
				{{0.6/0.6,0.6/0.6},{0.7/0.7,0.8/0.8},{0.5/0.5,0.4/0.5},{0.5/0.5,0.4/0.4},{0.4/0.4,0.5/0.5},{0.2/0.7,0.5/0.8},{1.2/1.0,1.5/1},{1.1/0.7,1.0/0.6},{1.2/1.2,1.1/1.1},{1.3/1.3,1.0/1.2},{1.0/1,1.0/1},{0.8/0.8,1.0/1}}			
		};
	
		//System.out.println(parse.getTimeStates()[2][0][0]);
		//System.out.println(parse.getBranchLengths()[0][0]);
		//System.out.println(parse.getPropStates()[2][0][0]);
		//System.out.println(parse.getPropStates()[2][0][0]);
		//assertEquals(prop[0][0][0], parse.getPropStates()[0][0][0]);
		//assertTrue(Math.abs(prop[0][0][0]-parse.getPropStates()[0][0][0])< 0.01);
		//assertTrue(Math.abs(prop[0][0][1]-parse.getPropStates()[0][0][1])< 0.01);
		//assertTrue(Arrays.equals(prop[0][0],parse.getPropStates()[0][0]));
		
		
		
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 12; j++) {	
				for (int k = 0; k < 2; k++) {
					//System.out.println("i"+i);
					//System.out.println("j"+j);
					//System.out.println("k"+k);
					assertTrue(Math.abs(prop[i][j][k]-parse.getPropStates()[i][j][k])< 0.0001);
				}
			}
		}
		
		
		prop = new double[][][] {
				{{2,2},{2,2},{2,2},{2,2},{2,2},{1+7.051235/8.18555,2},{2,2},{2,2},{1+4.563469/5.24651,2}},
				{{0,0},{0,0},{0,0},{0,0},{0,0},{1.134315/8.18555,0},{0,0},{0,0},{0.683041/5.24651,0}}	
		};
		
		for (int i = 0; i < 2; i++) {
			for (int j = 0; j < 9; j++) {	
				for (int k = 0; k < 2; k++) {
					//System.out.println("i"+i);
					//System.out.println("j"+j);
					//System.out.println("k"+k);
					assertTrue(Math.abs(prop[i][j][k]-parse1.getPropStates()[i][j][k])< 0.0001);
				}
			}
		}
		
	}
	
	public void testGetTheta() {
		//theta = new double[12][branchNum][C];
		double[][][] theta = new double[][][] {
				{{1/(0.4/0.6),1/(0.5/0.6)},{1/(0.2/0.7),1/(0.4/0.8)},{0,0},{0,0},{0,0},{0,0},{1/(0.5/1.0),1/(0.6/1)},{0,0},{0,0},{1/(0.4/1.3),1/(0.4/1.2)},{0,0},{0,0}},
				{{0,1/(0.5/0.6)},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{1/(1.5/1.2),1/(1.3/1.1)},{0,0},{0,0},{1/(0.4/0.8),1/(0.4/1.0)}},
				{{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}},
				{{0,0},{0,0},{0,0},{1/(0.3/0.5),1/(0.3/0.4)},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}},
				{{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{1/(1.2/1.2),1/(1.1/1.1)},{0,1/(2.4/1.2)},{1/(1.0/1),1/(1.0/1)},{0,0}},
				{{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{1/(1.4/1.3),1/(2.4/1.2)},{0,0},{0,0}},
				{{1/(1.2/0.6),2/(1.0/0.6)},{1/(1.4/0.7),1/(1.6/0.8)},{0,0},{0,0},{0,0},{0,0},{1,1},{0,0},{1/(0.9/1.2),1/(0.9/1.1)},{0,0},{0,0},{1/(1.2/0.8),1/(1.6/1)}},
				{{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,1/(1.0/1.2)},{0,0},{0,0}},
				{{0,0},{0,0},{1/(0.7/0.5),1/(0.8/0.5)},{0,0},{0,0},{0,0},{0,0},{0,0},{1/(0.9/1.2),1/(0.9/1.1)},{0,0},{0,0},{0,0}},
				{{0,0},{0,0},{1/(0.5/0.5),1/(0.4/0.5)},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0}},
				{{0,0},{0,0},{0,0},{0,0},{0,0},{1/(0.2/0.7),1/(0.5/0.8)},{1/(1.2/1.0),1/(1.5/1)},{1/(1.1/0.7),1/(1.0/0.6)},{0,0},{0,0},{0,0},{0,0}},
				{{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{0,0},{1/(1.3/1.3),1/(1.0/1.2)},{0,0},{0,0}}
		};
		
		//parse.printOutput();
		parse.printToScreen();
		/*
		double[][] cov = parse.getCovarianceMatrix(0);
		for (int i = 0; i < 12; i ++) {
			for (int j = 0; j < 12; j++) {
				System.out.print(cov[i][j]+" ");
			}
			System.out.println();
		}
		*/
		
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 12; j++) {	
				for (int k = 0; k < 2; k++) {
					//System.out.println("i"+i);
					//System.out.println("j"+j);
					//System.out.println("k"+k);
					assertTrue(Math.abs(theta[i][j][k]-parse.getTheta()[i][j][k])< 0.0001);
				}
			}
		}
	}
	
	/**
     * Tears down the test fixture.
     * (Called after every test case method.)
     */
    protected void tearDown() {
         parse = null;
    }

}
