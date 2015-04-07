/* This program is modified from TreenParser.java in the treejuxtaposer package from:
 http://olduvai.sourceforge.net/tj/ */


package parse_phylobayes;
import java.io.*;
import java.util.*;

/**
 * Parse the newick portion of a file
 * Note: this class is still not functioning correctly when outgroupFileName is "". (i.e. no outgroup)
 * @author James
 * @author Hui-Jie
 *
 */
public class TreeParser {
	/** Line (and tree information) termination. */
	private static final char lineTerminator = ';';
	
	private StreamTokenizer tokenizer;
    /**
     * Root node of the tree being parsed.  Must be initialized outside the tokenizer.
     */
    private TreeNode rootNode;
    
    /**
     * Store the number of taxa in the outgroup
     */
    private int outgroupNum;
    
    /**
     * Store the name of the taxon in the outgroup
     */
    private String[] outgroup;
    
    
    /**
     * Constructor
     * Initializes parsing of a tree by creating a tokenizer and setting default
     * properties (such as spacing, quoting characters). 
     *
     * Modified by Hui-Jie
     *
     * @param bufferedReader Buffered reader that could start in the middle of a nexus file or
     * the start of a newick file (basically the beginning of a newick tree, is run
     * for each tree in a nexus file)
     * @param outgroupFileName a string storing the outgroup species name
     */
    public TreeParser(StringReader b, String outgroupFileName)
    {
        tokenizer = new StreamTokenizer(b);
        tokenizer.eolIsSignificant(false); //Determines whether or not ends of line are treated as tokens.
        // Specifies that matching pairs of this character delimit string constants in this tokenizer.
        tokenizer.quoteChar('"'); //
//        tokenizer.quoteChar('\''); // TODO: check quote layering, quoted quotes
        
        //wordChar(,);
        //Specifies that all characters c in the range low <= c <= high are word constituents.
        tokenizer.wordChars('\'', '\''); // quote problem, turn this into a prime symbol?
        // 32 = space
        tokenizer.wordChars('!', '!'); // 33
        // 34 = "
        tokenizer.wordChars('#', '&'); // 35-38
        // 39-41 = '() newick
        tokenizer.wordChars('*', '+'); // 42-43
        // 44 = , newick
        tokenizer.wordChars('-', '/'); // 45-47
        // 48-59 = [0-9]:;
        tokenizer.wordChars('<', '<'); // 60
        // 61 = = nexus
        tokenizer.wordChars('>', '@'); // 62-64
        // 65-90 = [A-Z]
//        tokenizer.wordChars('[', '['); // 91 [ nexus comment character, treat as char
        // 92 = \ (esc, support esc'd spaces)
//      93 = ] nexus comment character
        tokenizer.wordChars('^', '^'); // 94
        // 95 = _ (underscore)
        tokenizer.wordChars('`', '`'); // 96
        // 97-122 = [a-z]
        tokenizer.wordChars('{', '~'); // 123-126
        // 127 = del
        
        if (!outgroupFileName.equals("")) {
        	//outgroup file
            /*
        	File f = new File(outgroupFileName);
            try
            {
                BufferedReader r = new BufferedReader(new FileReader(f));
                String line = r.readLine();
               //first line contains the number of taxa in the outgroup
                this.outgroupNum = Integer.parseInt(line);
                outgroup = new String[outgroupNum];
                for (int i = 0; i < outgroupNum; i++) {
                	line = r.readLine();
                	outgroup[i] = line;
                }
                
            }
            catch (FileNotFoundException e)
            {
                System.out.println("Couldn't find file: " + outgroupFileName);
            } catch (IOException e) {
    			// TODO Auto-generated catch block
    			e.printStackTrace();
    		}
    		*/
        	
        	//InputStream inStream = TreeParser.class.getClassLoader().getResourceAsStream(outgroupFileName);
        	InputStream inStream = this.getClass().getResourceAsStream(new File("../" + outgroupFileName).getPath().toString());
        	BufferedReader r = new BufferedReader(new InputStreamReader(inStream));
        	String line;
			try {
				line = r.readLine();
				//first line contains the number of taxa in the outgroup
	             this.outgroupNum = Integer.parseInt(line);
	             outgroup = new String[outgroupNum];
	             for (int i = 0; i < outgroupNum; i++) {
	             	line = r.readLine();
	             	outgroup[i] = line;
	             }
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			} finally {
	            try {
	                inStream.close();
	            } catch (IOException e) {
	                e.printStackTrace();
	            }
	        }
               
             
        	
        } else { //no outgroup
        	outgroupNum = 0;
        	outgroup = null;
        }
        
        
    }


	/**
     * Adds node at the top of the stack to the tree.  TreeNode is already created based
     * on Newick properties.
     * @param name Name of the node.
     * @param nodeStack Stack of nodes that haven't been added to the tree yet.  Nodes are popped when
     * they have names and all children are processed.
     * @return Newly added treeNode linked into the tree. 
     */
    private TreeNode popAndName(String name, Stack nodeStack)
    {
	    TreeNode topNode = (TreeNode)nodeStack.pop();
	    if (name == null)
	    {
	    	topNode.label = "";
	    	topNode.setName("");
	    }
	    else
	    {
	    	topNode.label = name;
	    	topNode.setName(name);
	    }
	    try //causing troubles when its parent is a pseudo node
	    {
	    	TreeNode parent = (TreeNode) nodeStack.peek();
	    	parent.addChild(topNode);
	    }
	    catch (EmptyStackException e)
	    {
	        if (topNode != rootNode)
	            System.out.println("Parser error on node " + topNode);
	    }
	    topNode.setExtremeLeaves(); // sets leftmost and rightmost leaf, non-recursive
	    topNode.setNumberLeaves(); // sets number of leaves, non-recursive
	    topNode.linkNodesInPreorder();
	    topNode.linkNodesInPostorder();
	    return topNode;
    }

    /**
     * PhyloBayes output tree tokenizer: converts a string (tree as a string) into a tree object.
     * The stream tokenizer should be initialized before calling this function.
     *
     * Modified by Hui-Jie
     *
     * @return Tree parsed from the stream.
     */
    public Tree tokenize()
    {
        final char openBracket = '(', closeBracket = ')', childSeparator = ',',
        	treeTerminator = lineTerminator, quote = '\'', doubleQuote = '"', 
        	infoSeparator = ':', stateSeparator = '_';
        //int progress = 0;
        rootNode = new TreeNode();
        Tree t = new Tree();
        t.setRootNode(rootNode);
        //t.setFileName(streamName);
        Stack<TreeNode> nodeStack = new Stack<TreeNode>();
        nodeStack.push(rootNode);
        int thisToken;
        TreeNode lastNamed = null;
        boolean EOT = false; //end of tree
        boolean nameNext = true;
        boolean stateNext = false;
        boolean isTrueNode = true; //whether lastNamed is a true node or a pseudo node
        int index = 0; //use to index nodes        
	try {
            while (EOT == false &&
                    (thisToken = tokenizer.nextToken()) != StreamTokenizer.TT_EOF)
            {
            switch (thisToken)
            {
                case stateSeparator:
                	if (nameNext) { //pop internal node whose children have been processed
                		lastNamed = popAndName(null, nodeStack);
                	}
                	nameNext = false;
                	stateNext = true;
                	isTrueNode = true;
                	break;
            	case StreamTokenizer.TT_WORD:
            		if(nameNext) { //name next
            			lastNamed = popAndName(tokenizer.sval, nodeStack);
            			nameNext = false;
            			//stateNext = true;
            		} else if (stateNext && isTrueNode) { //state for true node
            			lastNamed.setState(tokenizer.sval); //lastNamed is a true node, store state at lastNamed
            		} else if (stateNext && !isTrueNode) { //state for pseudo node
            			String tempState = tokenizer.sval;
            			//check if it's redundant
            			int tempToken = tokenizer.nextToken();
            			if (tempToken != infoSeparator) { //redundant state
            				nameNext = true;
            				stateNext = false;
            			} else { 
            				/*create pseudo node
            				if (lastNamed != null) {
            					System.err.println("Error: didn't expect last named node exist: " + tokenizer.sval);
            				} else {
            					TreeNode temp = new TreeNode();
            					temp.setState(tempState);
            					temp.setPseudo(true);
            					nodeStack.push(temp);
                				lastNamed = popAndName(null,nodeStack);
            				}
            				*/
            				lastNamed.addPathState(tempState);
            			}
            			//push back tempToken
            			tokenizer.pushBack();
            		} else {
            			System.err.println("Error: didn't expect this name/state here: " + tokenizer.sval);
            		}
            		
            		/*
            	    if (!nameNext)
            	        System.err.println("Error: didn't expect this name here: " + tokenizer.sval);
            	    lastNamed = popAndName(tokenizer.sval, nodeStack);
            		//progress += tokenizer.sval.length();
            		nameNext = false;
            		*/
            		break;
            	case StreamTokenizer.TT_NUMBER:
            		//check if it's in the form of `5e-05' instead of 0.00005
            		double time = tokenizer.nval;
            		int tempToken = tokenizer.nextToken();
            		if (tempToken != infoSeparator) { //`e-05' is a token
            			String temp = time + tokenizer.sval;             			
            			time = Double.parseDouble(temp);
            		} else {
            			tokenizer.pushBack();
            		}
            		            		
            		if (lastNamed != null) {
            			if(lastNamed.getPathStateSize() == 0){
            				lastNamed.setTime(time);
            			} else {
            				lastNamed.addPathTime(time);
            			}          			
            		} else {
            			System.err.println("Error: can't set value " + time + " to a null node");
            		}
            		//lastNamed = null;
            		nameNext = false;
            		stateNext = true;
            		isTrueNode = false;
            		/*
            		if (nameNext) //name contains number
            		    lastNamed = popAndName(tokenizer.sval, nodeStack);
            		else //time to parent node
            		{
            		    if (lastNamed != null)
            		        lastNamed.setTime(tokenizer.nval);
            		    else
            		        System.err.println("Error: can't set value " + tokenizer.nval + " to a null node");
            		    lastNamed = null;
            		}
            		//progress += (new Double(tokenizer.nval).toString()).length();
            		nameNext = false;
            		*/
            		break;
            	case infoSeparator:
            		/*
            	    if (nameNext)
            	        lastNamed = popAndName(null, nodeStack);
            	    //progress += 1;
            	    */
            		nameNext = false;
            		stateNext = true;
            	    break;
            	case treeTerminator:
            	case StreamTokenizer.TT_EOF:
            		if (nameNext)
            	        lastNamed = popAndName(null, nodeStack);
            	    EOT = true;
            	    //progress += 1;
            	    nameNext = false;
            	    break;
            	case openBracket:
            	    nodeStack.push(new TreeNode());
            	    //progress += 1;
            	    nameNext = true;
            	    stateNext = false;
            	    //isTrueNode = true;
            	    break;
            	case closeBracket:
            	    /*
            	     if (nameNext)
            	        lastNamed = popAndName(null, nodeStack);
            	    //progress += 1;
            	     */
            		nameNext = true;
            	    break;
            	case childSeparator:
            		/*
            	    if (nameNext)
            	        lastNamed = popAndName(null, nodeStack);
            	    nodeStack.push(new TreeNode());
            	    //progress += 1;
            	    */
            		nodeStack.push(new TreeNode());
            		nameNext = true;
            		stateNext = false;
            	    break;
            	default:
            	    debugOutput("default " + (char)thisToken);
            		break;
            }
        }
        }
        catch (IOException e) {
        }
        if (!nodeStack.isEmpty())
            System.err.println("Node stack still has " + nodeStack.size() + " things");
        t.postProcess();
        t.setNodeNum(); //need to call it first to ensure `equals' works properly.
        //remove outgroup if t is an unrooted tree
        if(outgroupNum != 0 && t.getRoot().numberChildren() == 3) {
        	t.removeOutgroup(outgroup);
        	//t.setNodeNum();
        }
        
        return t;
    }

    /**
     * Debug printout function.  Avoid using the system calls and use this, and set flag
     * {@link #debugOutput} depending on debugging or not.
     * @param s Display the string, for debugging.
     */
    public void debugOutput(String s)
    {
            System.out.println(s);
    }


	/**
	 * @param args[0] pb_mpi output filename
	 * @param args[1] outgroup filename
	 */
/*    
	public static void main(String[] args) {
        String fileName = args[0];
        long start = System.currentTimeMillis();
        String outgroupFileName = args[1];        
        File f = new File(fileName);
        try
        {
            BufferedReader r = new BufferedReader(new FileReader(f));
            String line = r.readLine();
            while(line != null) { //does not work if there is a blank line between tree! 
            	                 // "" can be read in and it is not a null string object.
            	StringReader str = new StringReader(line);
                TreeParser tp = new TreeParser(str, outgroupFileName);
                Tree t = tp.tokenize();
                line = r.readLine();
            }
        }
        catch (FileNotFoundException e)
        {
            System.out.println("Couldn't find file: " + fileName);
        } catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
        System.out.println("Parsed in " + ((System.currentTimeMillis() - start)/1000.0) + " s");
        System.exit(0);

	}
	*/
    
}
