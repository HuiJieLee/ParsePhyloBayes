/* This program is modified from TreenNode.java in the treejuxtaposer package from:
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
import java.util.*;



/**
 * This class respresnets a node of a tree and store info about node.
 * Nodes have fields that store pre- and post-order nodes, 
 * list of children, parent node, node number, left / right most leaf,
 * state of the node (A,C,T,G), state changes (path) to its parent node,
 * and the times of substitution events.
 * @author Tamara Munzner, Li Zhang, Yunhong Zhou
 * @author Hui-Jie
 *
 */
public class TreeNode {
	
	/** Array of child nodes that are attached below this internal node.  Null if this is a leaf. */
	protected ArrayList<TreeNode> children;

	/** key is unique for nodes in one tree.  Keys are pre-ordered (root = 0, depth-traversal ordering). */
	private int key;
	
	/** nodeNum is unique for nodes in one tree. They are ordered as the same way in MultidivTime.*/
	private int nodeNum;

	/** The parent of this node.  This is null for the root node. */
	private TreeNode parent;

	/**
	 * Node name with default "". Most internal nodes have no name and all leaf
	 * nodes have a name.  This becomes the long version of the node name when fully
	 * qualified names are used.
	 */
	protected String name = ""; // the long form in fully qualified names

	/** Distance from this node to the root node. The root is at height 1. */
	private int height;
	
	/** Leftmost (minimum) leaf node under this internal node (or this node for leaves). */
	private TreeNode leftmostLeaf;
	/** Rightmost (maximum) leaf node under this internal node (or this node for leaves). */
	private TreeNode rightmostLeaf;

	/** The number of leaves under this internal node (or 1 for leaves). */
	private int numberLeaves;

	/** The next preorder node. */
	private TreeNode preorderNext = null;

	/** The next postorder node. */
	private TreeNode posorderNext = null;
	
	/** Store the state (ACTG  or codon state) of the node */
	private String state = null;
	
	/** Store the time to parent node */
	private double time;

	public String label;
	
	/** Store the state path to its parent node*/
	public ArrayList<String> pathState; 
	
	/** Store the time path to its parent node */
	public ArrayList<Double> pathTime;

	/**
	 * Default tree node constructor.
	 * Children list initially set to capacity 2 as in most case binary.
	 * 	Used in 2 places: create the root when creating the tree;
	 *  the parser uses this to create nodes attached to the root.
	 */
	public TreeNode() {
		children = new ArrayList<TreeNode>();
		pathState = new ArrayList<String>();
		pathTime = new ArrayList<Double>(); 
	}
	
	/**
	 * Set the name for this node, the name is usually the label drawn with this node.
	 * @param s The new value of {@link #name}, the name for this node.
	 */
	public void setName(String s) {
		name = s;
	}

	/**
	 * Get the number of children under this node.
	 * @return Number of nodes stored in the children array {@link #children}.
	 */
	public int numberChildren() {
		return children.size();
	}

	/**
	 * Get a given child for this node, with range checking and casting.
	 * @param i The child index to get.
	 * @return The i(th) child for this node.
	 */
	public TreeNode getChild(int i) {
		if (i < children.size())
			return (TreeNode) children.get(i);
		else
			return null;
	}

	/**
	 * Tests to determine if this node is a leaf.  Does not work for nodes not in the tree structure.
	 * @return True if this node has no linked children, and therefore is a leaf node for the tree.
	 */
	public boolean isLeaf() {
		return children.isEmpty();
	}

	/**
	 * Tests to determine if this node is the root of its tree. Does not work for nodes not in the tree structure.
	 * @return True if this node has no linked parent, and therefore is the root of the tree.
	 */
	public boolean isRoot() {
		return (null == parent);
	}

	/**
	 * Tests nodes for equality, based on the name of the node.
	 * Can only be used after every node get unique NodeNum.
	 * @param n Second node to test vs. this node.
	 * @return True if the names of both nodes are the same, false otherwise.
	 */
	public boolean equals(TreeNode n) {
		return (nodeNum == n.getNodeNum());
	}

	/**
	 * Add a child to the end of the list of children.  Note there is no remove child method, this is permanent.
	 * Additional processing for linking nodes (setting up pointers and leaf properties, for example) is done later.
	 * @param n New child node for this node.
	 */
	public void addChild(TreeNode n) {
		children.add(n);
		n.parent = this;
	}
	
	public void removeChild(TreeNode n) {
		children.remove(n);
		n.parent = null;
	}
	/**
	 * Get the parent for this node.
	 * @return Value of {@link #parent}.
	 */
	public TreeNode parent() {
		return parent;
	}

	/** Get the first child of this node. Doesn't work with leaf nodes.
	 * @return First child of this internal node.
	 */
	protected TreeNode firstChild() {
		return (TreeNode) children.get(0);
	}

	/** Get the last child of this node. Doesn't work with leaf nodes.
	 * @return Last child of this internal node.
	 */
	public TreeNode lastChild() {
		return (TreeNode) children.get(children.size() - 1);
	}

	/**
	 * Returns the key for this node.
	 * @return The value of {@link #key} for this node.
	 */
	public int getKey() {
		return key;
	}
	
	/**
	 * Set the key for the node
	 * @param 
	 * @return 
	 */
	public void setKey(int key) {
		this.key = key;
	}

	/**
	 * Returns the label for this node, which is {@link #name}.
	 * @return The value of {@link #name} for this node.
	 */
	public String getName() {
		return name;
	}

	/**
	 * Set the extreme leaves for this node.  This is done in leaf->root direction, so all linking can be done in O(n) time.
	 *
	 */
	public void setExtremeLeaves() {
		if (isLeaf()) {
			leftmostLeaf = this;
			rightmostLeaf = this;
			return;
		}
		leftmostLeaf = firstChild().leftmostLeaf;
		rightmostLeaf = lastChild().rightmostLeaf;
	}

	/** root->leaf traversal, depth first in direction of leftmost leaf. */
	public void linkNodesInPreorder() {
		if (isLeaf())
			return;
		preorderNext = firstChild();
		for (int i = 0; i < numberChildren() - 1; i++)
			getChild(i).rightmostLeaf.preorderNext = getChild(i + 1);
		// rightmostLeaf.preorderNext = null; // redundant
	}

	/** Leaf->root traversal, starting at leftmost leaf of tree. */
	public void linkNodesInPostorder() {
		if (isLeaf())
			return;
		// n.posorderNext = null; // redundant
		for (int i = 0; i < numberChildren() - 1; i++)
			getChild(i).posorderNext = getChild(i + 1).leftmostLeaf;
		lastChild().posorderNext = this;
	}

	/**
	 * Sets the number of leaves, must be run on leaves first (pre-order)
	 * 
	 * @return The number of leaves ({@link #numberLeaves}) including the
	 *         current node (leaves = 1)
	 */
	public int setNumberLeaves() {
		numberLeaves = 0;
		if (isLeaf())
			numberLeaves = 1;
		else
			for (int i = 0; i < children.size(); i++)
				numberLeaves += getChild(i).numberLeaves;
		return numberLeaves;
	}

	/**
	 *  Returns the height of the node.
     *  Add by Hui-Jie
	 */
	public int getHeight() {
		return height;
	}
	
	/**
	 *  Set height of the node
     *  Add by Hui-Jie
	 * @param height
	 */
	public void setHeight(int height) {
		this.height = height;
	}
	
	/**
	 *  Return the left most leaf
     *  Add by Hui-Jie
	 */
    public TreeNode getLeftmostLeaf() {
    	return leftmostLeaf;
    }
    
	/**
	 *  Return the right most leaf
     *  Add by Hui-Jie
	 */
    public TreeNode getRightmostLeaf() {
    	return rightmostLeaf;
    }
	
    /**
     *  Return the number of leaves under this internal nodes, may be zero.
     *  Add by Hui-Jie
     */
    public int getNumberLeaves() {
    	return numberLeaves;
    }
	

    /**
     *  Return the next preorder node
     *  Add by Hui-Jie
     */
    public TreeNode getPreorderNext() {
    	return preorderNext;
    }
    
    /**
     *  Setter for the preorder
     *  Add by Hui-Jie
     * @param preorderNext
     */
    public void setPreorderNext(TreeNode next) {
    	this.preorderNext = next;
    }
    
    /**
     *  Setter for the postorder
     *  Add by Hui-Jie
     * @param postorderNext
     */
    public void setPostorderNext(TreeNode next) {
    	this.posorderNext = next;
    }
    
    /**
     *  Return the next postorder node
     *  Add by Hui-Jie
     */
    public TreeNode getPostorderNext() {
    	return posorderNext;
    }
    
    /**
     *  Set the state of the node
     *  Add by Hui-Jie
     * @param String state
     */
    public void setState(String state) {
    	this.state = state;
    }
    
    /**
     *  Return the state of the node
     *  Add by Hui-Jie
     */
    public String getState() {
    	return state;
    }
    
    /**
     *  Set the time
     *  Add by Hui-Jie
     * @param double time
     */
    public void setTime(double time) {
    	this.time = time;
    }
    
    /**
     *  Return the time to its parent node
     *  Add by Hui-Jie
     */
    public double getTime() {
    	return time;
    }
    
    
    /**
     *  Add path state
     *  Add by Hui-Jie
     * @param state
     */
    public void addPathState(String state) {
    	pathState.add(state);
    }
    
    /**
     *  Return size of pathState
     *  Add by Hui-Jie
     *  @return size
     */
    public int getPathStateSize() {
    	return pathState.size();
    }
    
    
    /**
     *  Return size of pathTime
     *  Add by Hui-Jie
     *  @return size
     */
    public int getPathTimeSize() {
    	return pathTime.size();
    }
    
    /**
     *  Add path time
     *  Add by Hui-Jie
     *  @param time
     */
    public void addPathTime(double time) {
    	pathTime.add(time);
    }

    /**
     *  Return the whole path of state changes starting from the parent state to the state of this node.
     *  Add by Hui-Jie
     *  @return state path
     */
    public ArrayList<String> getPathState() {
    	ArrayList<String> states = new ArrayList<String>();
    	states.add(this.parent.getState());
    	Collections.reverse(pathState);
    	states.addAll(pathState);
    	Collections.reverse(pathState);
    	states.add(state);
    	return states;
    }
    
    /**
     *  Return the whole path of times starting from the parent state to the current node.
     *  Add by Hui-Jie
     *  @return time path
     */
    public ArrayList<Double> getPathTime() {
    	ArrayList<Double> times = new ArrayList<Double>();
    	Collections.reverse(pathTime);
    	times.addAll(pathTime);
    	Collections.reverse(pathTime);
    	times.add(time);
    	return times;
    }
    
    /**
     *  Return the node number
     *  Add by Hui-Jie
     * @return node number
     */
	public int getNodeNum() {
		return nodeNum;
	}

    /**
     *  Set node number
     *  Add by Hui-Jie
     */
	public void setNodeNum(int nodeNum) {
		this.nodeNum = nodeNum;
	}
	
	/**
     *  Count number of state changes on a branch. Assume each state is a character: {A, G, C, T}.
     *  Specify which kind of change it is
     *  Add by Hui-Jie
     *  @param from starting state
     *  @param to ending state
     *  @return number of state changes
     */
    public int getNumberOfChanges(String from, String to) {
    		String start = parent.getState();
    		String end = null;
    		int count = 0;
    		
    		if(pathState.size()==0) {
    			if(start.equals(from) && getState().equals(to)) 
    				count = 1;
    		} else {
    			for (int i = pathState.size() - 1; i >= 0; i--) {
        			end = pathState.get(i);
        			if (start.equals(from) && end.equals(to)) {
        				count ++;
        			}
        			start = pathState.get(i);
        		}   	
    			if(start.equals(from) && getState().equals(to)) count++;
    		}
    	return count;
    }
    
    /**
     *  Count time duration for a state on a branch. Assume each state is a character: {A, G, C, T}
     *  Specify which type of state it is.
     *  Add by Hui-Jie
     *  @param state type
     *  @return time duration of that state type
     */
    public double getTimeOfState(String type) {
    	double timeState = 0;
    	String start = parent.getState();
    	//if no path for this node, pathState.size() = 0, skip for loop
    	if (pathState.size() > 0) {
    		for (int i = (pathState.size()-1); i>= 0; i--) {
        		if (start.equals(type)) {
        			timeState += pathTime.get(i);
        		}
        		start = pathState.get(i);
        	}
    	}
    	
    	if(start.equals(type)) {
    		timeState += getTime();
    	}
    	
    	return timeState;
    }
    
    /**
     *  Calculate the branch length (time) of the branch ending with this node.
     *  Sum all path time.
     *  Add by Hui-Jie
     *  @return branch length
     */
    public double getBranchLength() {
    	double br = time;
    	if(pathState.size() != 0) {
    		for (int i = 0; i < pathState.size(); i++) {
    			br += pathTime.get(i);
    		}
    	} 
    	return br;
    }
    
}
