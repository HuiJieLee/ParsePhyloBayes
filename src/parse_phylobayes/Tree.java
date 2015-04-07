/* This program is modified from Tree.java in the treejuxtaposer package from:
   http://olduvai.sourceforge.net/tj/ */

/*
 Copyright (c) 2002 Compaq Computer Corporation
 
 SOFTWARE RELEASE
 
 Permission is hereby granted, free of charge, to any person obtaining
 a copy of this software and associated documentation files (the
 "Software"), to deal in the Software without restriction, including
 without limitation the rights to use, copy, modify, merge, publish,
 distribute, sublicense, and/or sell copies of the Software, and to
 permit persons to whom the Software is furnished to do so, subject to
 the following conditions:
 
 - Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.
 
 - Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.
 
 - Neither the names of Compaq Research, Compaq Computer Corporation
 nor the names of its contributors may be used to endorse or promote
 products derived from this Software without specific prior written
 permission.
 
 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 IN NO EVENT SHALL COMPAQ COMPUTER CORPORATION BE LIABLE FOR ANY CLAIM,
 DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
 THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */


package parse_phylobayes;
import java.text.Collator;
import java.util.*;


/**
 * Store the tree structure
 * Additional functions are added to the original Tree.java
 * @author Yunhong Zhou, Li Zhang
 * @author Hui-Jie
 *
 */
public class Tree {
	/** The list of nodes of the tree indexed by their keys  */ 
	public ArrayList<TreeNode> nodes; 


	/** 
	 * Assign nodeName to each node
	 */
	private HashMap<String, TreeNode> nodesByName; 
	
	/**
	 * Assign nodeNum to each node
     * Added by Hui-Jie
	 */
	private Map<Integer, TreeNode> nodesByNodeNum;

	/** key should be unique for each tree, set by object that creates trees  */
	private int key;

	/** Leaf counter, for determining grid size, making arrays for tree comparisons */
	private int numLeaves = 0;
	
	/**
	 * Root node of this tree
	 */
	protected TreeNode root = null;
	
	/**
	 * Height of tree, which is also the longest path from the root to some leaf node.
	 */
	private int height = 0;
	
	/**
	 * Default tree constructor.  Nodes are created by parser and added in later.
	 *
	 */
	public Tree() {
		root = new TreeNode();
		nodes = new ArrayList<TreeNode>();
		nodesByName = new HashMap<String, TreeNode>();
		nodesByNodeNum = new HashMap<Integer, TreeNode>();
	}

	/**
	 * Returns the number of interior nodes in this tree.
	 * @return Total number of nodes minus the number of leaves.
	 */
	private int getInteriorCount() { return nodes.size() - numLeaves;}
	/**
	 * Returns the node count, for internal and leaf nodes.
	 * @return Size of the {@link #nodes} array, which contains all nodes.
	 */
	protected int getTotalNodeCount() { return nodes.size();}
	
	/**
	 * Return the number of branches on this tree.
     * Added by Hui-Jie
	 * @return Number of branches.
	 */
	public int getNumBranches() { return nodes.size()-1;}
	
	/**
	 * Return the number of leaves on this tree
     * Added by Hui-Jie
	 * @return number of leaves.
	 */
	public int getNumLeaves() {
		return numLeaves;
	}
	
	
	/**
	 * Accessor for height of tree.  This is also the longest path from the root to some leaf node.
	 * @return value of {@link #height}.
	 */
	public int getHeight() { return height; }

	/**
	 * Returns the node indexed by the given key.
	 * @param key Key of the node to retrieve.
	 * @return Treenode referenced by the given key.
	 */
	public TreeNode getNodeByKey(int key){ if (key >= nodes.size()) return null; return (TreeNode) nodes.get(key);}
	
	
	/**
	 * Returns the node indexed by nodeNum.
     * Added by Hui-Jie
	 * @param nodeNum node numbered same as Multidivtime
	 * @return Treenode referenced by the given nodeNum
	 */
	public TreeNode getNodeByNodeNum(int nodeNum) {
		if (nodeNum >= nodes.size()) return null; return (TreeNode) nodesByNodeNum.get(nodeNum);
	}
	
	/**
	 * Returns the node given by the string.
	 * @param s Name/label of node to retrieve.
	 * @return Treenode referenced by the given name.
	 */
	public TreeNode getNodeByName(String s){ 
		return (TreeNode) nodesByName.get(s);
	}

	/** Left most leaf accessor.  This is the "min leaf"
	 * @return root's left most leaf, which is the smallest indexed leaf node in the tree.
	 */
	public TreeNode getLeftmostLeaf() { return root.getLeftmostLeaf(); }
	/** Root accessor.
	 * @return Value of {@link #root}*/
	public TreeNode getRoot() { return root;}
	public void setRootNode(TreeNode newRoot) { this.root = newRoot; }

	/**
	 * Post processing includes computing size of each node, 
	 * linking nodes in different order, etc.
	 * Sets left and right-most leaves of the tree.
	 * Computes and stores pre- and post-orders for leaves.
	 * Can't do minmax until after linkNodesInPreorder is called 
	 * to set index values!
	 *
	 * @see     TreeNode
	 */
	public void postProcess() {
		preorderPostProcess();
		linkLeaves();
//		System.out.println("progress bar updated: min:" + jpb.getMinimum() + " max:" + jpb.getMaximum() + " value:" + jpb.getValue());
	}

	/**
	 * 
	 * Traverses the tree in pre-order, stores the ordering in the preorderNext field of TreeNodes
	 * Sets node count for the tree.
	 *
	 * @see     TreeNode
	 */
	private void preorderPostProcess()
	{
		// munge names here, names become fully qualified, labels are what names were
		final char separator = '/'; // separator between name fields
		// arbitrary seen by users in search, no parsing on this is required later
		int index = 0;
		height = 1;
		for(TreeNode n = root; n != null; n = n.getPreorderNext())
		{
			//n.key = index++;
			n.setKey(index++);
			nodes.add(n);
			if(n.name != null && n.name.length() > 0) {
				// don't put an empty string in the
				// hash table
				nodesByName.put(n.name, n);
			}
			//n.height = (null != n.parent) ? n.parent.height+1 : 1;
			if(n.parent() != null) { //not root
				n.setHeight(n.parent().getHeight()+1);
			} else {
				n.setHeight(1);
			}
			
			//height = (n.height > height) ? n.height : height;
			if(n.getHeight() > height) {
				height = n.getHeight();
			} 
		}

	}

	/**
	 * Traverse the tree and initialize the {@link #nodesByName} and {@link #nodes} data structures.
	 * Used when modifying the names of nodes as well as initialization.
	 *
	 */
	public void setUpNameLists()
	{
		nodes = new ArrayList<TreeNode>();
		nodesByName = new HashMap<String,TreeNode>();
		final char separator = '/'; // separator between name fields
		for(TreeNode n = root; n != null; n = n.getPreorderNext())
		{
			nodes.add(n);
			if(n.name != null && n.name.length() > 0) {
				// don't put an empty string in the
				// hash table
				nodesByName.put(n.name, n);
			}
			//n.height = (null != n.parent) ? n.parent.height+1 : 1;
			if (n.parent() != null) {
				n.setHeight(n.parent().getHeight()+1);
			} else {
				n.setHeight(1);
			}
			//height = (n.height > height) ? n.height : height;
			if (n.getHeight() > height) {
				height = n.getHeight();
			} 
		}
	}
	

	/**
	 * Wrapper for initiating {@link #linkSubtreeNodesInPreorder(TreeNode)} with the root node.
	 */
	private void linkNodesInPreorder() {

		linkSubtreeNodesInPreorder(root);

	}

	/**
	 * Traverses the subtree rooted at TreeNode n in pre-order, stores the
	 * ordering in the preorderNext field of TreeNodes. 
	 * @param   n the root of the subtree
	 *
	 * @see     TreeNode
	 */
	private void linkSubtreeNodesInPreorder(TreeNode n) {

		if(n.isLeaf()) return;
		for(int i=0; i<n.numberChildren(); i++) {
			linkSubtreeNodesInPreorder(n.getChild(i));
		}

		n.setPreorderNext(n.firstChild());
		for(int i = 0; i < n.numberChildren()-1; i++) {
			//n.getChild(i).rightmostLeaf.preorderNext = n.getChild(i+1);
			n.getChild(i).getRightmostLeaf().setPreorderNext(n.getChild(i+1));
		}
		//n.rightmostLeaf.preorderNext = null;
		n.getRightmostLeaf().setPreorderNext(null);
	}

	/**
	 * 
	 * Links leaves of the tree in pre-order,
	 * check to see whether leaves have distinct names.
	 * If leaves have the same name, add a suffix index separated by " "
	 *
	 * @see     #linkNodesInPreorder()
	 * @see     TreeNode
	 * @see     NameComparator
	 * @param jpb Progress bar.
	 */
	private void linkLeaves() {
		int counter = 0;
		int percentage = 0;
		TreeNode pren = root.getLeftmostLeaf();
		Vector leaves = new Vector<TreeNode>();
		leaves.add(pren);
//		pren.lindex = 0;
		for(TreeNode n = pren.getPreorderNext(); n!=null; n=n.getPreorderNext())
		{
			counter++;
			if(n.isLeaf())
			{
				leaves.add(n);
			}
		}
		numLeaves = leaves.size();

		NameComparator myNameComparator = new NameComparator();
		TreeNode[] sortedLeafArray = (TreeNode[])leaves.toArray(new TreeNode[leaves.size()]);
		Arrays.sort(sortedLeafArray, myNameComparator);
		int index = 0;
		TreeNode curr = sortedLeafArray[0];
		TreeNode next;
		for(int i=0; i<leaves.size()-1; i++){
			next = sortedLeafArray[i+1]; // only 1 index lookup per iteration
			boolean compare = myNameComparator.compare(curr, next) == 0; 
			if (compare || index > 0)
			{
				String name = curr.getName();
				nodesByName.remove(curr); // before all nodes with
				// same name were being ignored in search and comparing two identically named
				// leaves was broken, much fewer differences in trees with many leaves that
				// have the same name (imagine: all index.html occurences being marked as
				// different since numbering convention doesn't string match the original node name)
				curr.setName(name+ " " + index); //sb.toString());
				nodesByName.put(name+ " " + index, curr);//sortedLeafArray[i].getName(), sortedLeafArray[i]); // add the node back with number convention
				if (!compare)
					index = 0;
				else
					index++;
			}
			curr = next;
		}
	}

	/**
	 * Get the leaves under this node.  Used for tree to tree comparison, removing leaf nodes from difference calculations when they only appear in one side of the tree.
	 * This operation is constant time per leaf, since it relies on pre-ordered node links and pointers to extreme leaves.
	 * Time complexity of this function is linear in the number of leaves in the subtree under the node.
	 * @param node Node to get leaves under.  The root node will return all leaves in the tree, leaves return a list of just themselves.
	 * @return List of leaves under this node.
	 */
	public LinkedList getLeaves(TreeNode node)
	{
		LinkedList<TreeNode> leaves = new LinkedList<TreeNode>();
		TreeNode currNode = node.getLeftmostLeaf();
		while (currNode != node.getRightmostLeaf())
		{
			if (!currNode.isLeaf()) // internal node?
				currNode = currNode.getLeftmostLeaf(); // descend to minimal leaf
			leaves.add(currNode);
			currNode = currNode.getPreorderNext();
		}
		leaves.add(node.getRightmostLeaf());
		return leaves;
	}

	/**
	 * Return the key
	 * @return key
	 */
	public int getKey() {
		return key;
	}
    
	/**
	 * Set the key
	 * @param key
	 */
	public void setKey(int key) {
		this.key = key;
	}

	/**
	 * Renumbering nodes by setting nodeNum so that it is indexed the same as Multidivtime.
	 * Set up the {@link #nodesByNodeNum} and {@link #nodes} data structures.
     * Added by Hui-Jie
	 */
	public void setNodeNum() {
		int index = 0;
		for (int i = 0; i < nodes.size(); i++) {
			if(nodes.get(i).isLeaf()) {
				nodes.get(i).setNodeNum(index);
				nodesByNodeNum.put(index, nodes.get(i));
				index++;
			}
		}
		for (int j = (nodes.size()-1); j >= 0; j--) {
			if(!nodes.get(j).isLeaf()) {
				nodes.get(j).setNodeNum(index);
				nodesByNodeNum.put(index, nodes.get(j));
				index++;
			}
		}			
	}
	
	/**
	 * Remove the node which is the most recent common ancestor of the taxa given.
	 * Note: this node is the child of the `root' of the unrooted tree. (root has degree 3) 
	 * Also remove all the nodes above this node.
     * Added by Hui-Jie
	 * @param taxa names of the taxa to be removed 
	 */
	public void removeOutgroup(String[] taxa) {
		if (root.numberChildren() != 3) {
			System.out.println("Input tree is not an unrooted tree.");
			System.out.println("Cannot root the tree by an outgroup");
		} else { //root has three children, one of them is the outgroup and should be remove			
			int index = -1; //use to find the index of the outgroup
			boolean find = false; //find the outgroup or not
			while (!find) {
				index++;
				if (index >= 3) {
					System.out.println("Error in the outgroup file. Please check again.");
					break;
				}
				for (int i = 0; i < taxa.length; i++) {
					if(root.getChild(index).getLeftmostLeaf().getName().equals(taxa[i])) find = true;
					if(root.getChild(index).getRightmostLeaf().getName().equals(taxa[i])) find = true;
				}
			}
			switch(index) {
			case 0: //outgroup is the first child
				root.setPreorderNext(root.getChild(1));
				break;
			case 1: //outgroup is the second child
				root.getChild(0).getRightmostLeaf().setPreorderNext(root.getChild(2));
				root.getChild(0).setPostorderNext(root.getChild(2).getLeftmostLeaf());
				break;
			case 2: //outgroup is the third child
				root.getChild(1).getRightmostLeaf().setPreorderNext(null);
				root.getChild(1).setPostorderNext(root);
				break;
			default: //error in the outgroup file. do nothing
				break;
			}
			
			removeNodes(root.getChild(index));
			root.removeChild(root.getChild(index));
			root.setNumberLeaves();
			root.setExtremeLeaves();
			numLeaves = numLeaves - taxa.length;
			nodesByNodeNum = new HashMap<Integer, TreeNode>();
			setNodeNum();
		}
	}
	
	/**
	 * Recursively remove nodes from the tree
     * Added by Hui-Jie
	 * @param base any nodes above the base will be removed.
	 */
	private void removeNodes(TreeNode base) {
		if (base.isLeaf()) {
			nodesByName.remove(base.getName());
			nodes.remove(base);
		} else { //no name for internal node
			removeNodes(base.firstChild());
			removeNodes(base.lastChild());
			nodes.remove(base);
		}
	}
	
	
	/**
	 * Compute number of state changes on each branch. Assume each state is a character: {A, G, C, T}.
	 * Branches are numbered by the node that ends the branch.
	 * # of branches = # of nodes - 1.
	 * Root has the last index, and it does not correspond to a branch.
     * Added by Hui-Jie
	 * @param from starting state
     * @param to ending state
     * @return an array of number of state changes 
	 */
	public int[] getNumberOfChanges(String from, String to) {
		int numBranch = nodes.size()-1;
		int changes[] = new int[numBranch];
		for (int i = 0; i < numBranch; i++) {
			changes[i] = getNodeByNodeNum(i).getNumberOfChanges(from, to);
		}
		return changes;
	}
	
	/**
	 * Return time duration of a state on each branch. Assume state can be : {A, G, C, T}
	 * Branches are numbered by the node that ends the branch.
	 * # of branches = # of nodes - 1.
	 * Root has the last index, and it does not correspond to a branch.
     * Added by Hui-Jie
	 * @param type state
     * @return an array of time duration 
	 */
	public double[] getTimeOfState(String type){
		int numBranch = nodes.size()-1;
		double times[] = new double[numBranch];
		for (int i = 0; i < numBranch; i++) {
			times[i] = getNodeByNodeNum(i).getTimeOfState(type);
		}
		return times;
	}
	
	/**
	 * Compute the branch length (time) of a branch.
	 * Branches are numbered by the node that ends the branch.
	 * # of branches = # of nodes - 1.
	 * Root has the last index, and it does not correspond to a branch.
     * Added by Hui-Jie
	 * @return an array of branch length
	 */
	public double[] getBranchLength() {
		int numBranch = nodes.size()-1;
		double br[] = new double[numBranch];
		for (int i = 0; i < numBranch; i++) {
			br[i] = getNodeByNodeNum(i).getBranchLength();
		}
		return br;
	}
	
	/**
	 * Compute the proportion of time that a state exists on each branch.
	 * * Branches are numbered by the node that ends the branch.
	 * # of branches = # of nodes - 1.
	 * Root has the last index, and it does not correspond to a branch.
     * Added by Hui-Jie
	 * @param type state
     * @return an array of time proportion
	 */
	public double[] getPropState(String type) {
		int numBranch = nodes.size()-1;
		double prop[] = new double[numBranch];
		for (int i = 0; i < numBranch; i++) {
			prop[i] = getNodeByNodeNum(i).getTimeOfState(type)/getNodeByNodeNum(i).getBranchLength();
		}
		return prop;
	}
 
	
}

/** Comparator class for Strings */
class NameComparator implements Comparator{
	/** collator object used for string comparison. */
	Collator myCollator = Collator.getInstance(Locale.US);

	/** String comparator, uses {@link Collator} comparator. */
	public int compare(Object o1, Object o2){
		String s1 = ((TreeNode)o1).getName();
		String s2 = ((TreeNode)o2).getName();
		return myCollator.compare(s1, s2);
	}

}

