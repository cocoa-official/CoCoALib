#ifndef CoCoA_JBDatastructure_H
#define CoCoA_JBDatastructure_H

//   Copyright (c)  2015  John Abbott,  Anna M. Bigatti
//   Original author: 2015  Mario Albert

//   This file is part of the source of CoCoALib, the CoCoA Library.
//
//   CoCoALib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   CoCoALib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with CoCoALib.  If not, see <http://www.gnu.org/licenses/>.

#include <list>
#include "CoCoA/SparsePolyOps-RingElem.H"
//#include "CoCoA/SparsePolyRing.H" // from SparsePolyOps-RingElem.H

namespace CoCoA {
  namespace Involutive {

    //forward Declaration
    class JanetIterator;

    /****************************************************************************************************************************/
    /* This class stores the JanetTriple. It contains the myPolynomial
     * p, its ancestor, and its already prolonged variables
     */
    /*****************************************************************************************************************************/
    class JanetTriple {
    public:
      /**
       * The constructor of the JanetTriple.
       * @param pol a RingElem for the polynomial
       * @param anc a PPMonoid for the ancestor
       * The JanetTriple doesn't contain any already prolonged variables
       **/
      JanetTriple(ConstRefRingElem pol, ConstRefPPMonoidElem anc)
          : myPolynomial(pol),
            myAncestor(anc),
            myAlreadyProlongedVars(NumIndets(owner(pol)), false)
      {
      }

      /**
       * The constructor of the JanetTriple.
       * @param pol a RingElem for the polynomial
       * The ancestor is the LPP of polynomial.
       * The JanetTriple doesn't contain any already prolonged variables
       **/
      JanetTriple(ConstRefRingElem pol)
          : myPolynomial(pol),
            myAncestor(LPP(pol)),
            myAlreadyProlongedVars(NumIndets(owner(pol)), false)
      {
      }

      /** The constructor of the JanetTriple
       *	@param pol the new polynomial
       *	@param anc the new ancestor
       *	@param nm the already prolonged variables as std::vector<bool>
       */
      JanetTriple(ConstRefRingElem pol, ConstRefPPMonoidElem anc, std::vector<bool> nm)
          : myPolynomial(pol),
            myAncestor(anc),
            myAlreadyProlongedVars(nm)
      {
      }

      /**
       * This function changes the polynomial of the JanetTriple
       * @param pol the new polynomial
       */
      inline void mySetPol(ConstRefRingElem pol) {
        CoCoA_ASSERT(LPP(myPolynomial) == LPP(pol));
        myPolynomial = pol;
      }

      /**
       * This function changes the ancestor of the JanetTriple
       * @param anc the the ancestor
       */
      inline void mySetAnc(ConstRefPPMonoidElem anc) {
        myAncestor = anc;
      }

      /**
       * This function sets a already prolonged variable
       * @param index is the index of the variable which was already prolonged
       * If the variable is already prolonged, nothing will happen
       */
      inline void myVarIsProlonged(const long& index) {
        myAlreadyProlongedVars[index] = true;
      }

      /**
       * This function sets a already prolonged variable to not prolonged
       * @param index becomes not prolonged
       * If the variable is not prolonged, nothing will happen
       */
      inline void myVarIsNotProlonged(const long& index) {
        myAlreadyProlongedVars[index] = false;
      }

      /**
       * LPP of myPolynomial
       */
      ConstRefPPMonoidElem myGetLppPol() const {
        return static_cast<const SparsePolyRingBase*>(owner(myPolynomial).myRawPtr())->myLPP(raw(myPolynomial));
        // return const_cast<SparsePolyRingBase*>(owner(myPolynomial))->myLPP(raw(myPolynomial));
        // return AsSparsePolyRing(owner(myPolynomial))->myLPP(raw(myPolynomial));
      }

      /**
       * This function returns the current polynomial as RingElem
       */
      inline RingElem& myGetPol() {
        return myPolynomial;
      }

      inline const RingElem& myGetPol() const {
        return myPolynomial;
      }

      /**
       *This function returns a pointer to the current polynomial
       */
      inline RingElem* myGetPolPtr() {
        return &myPolynomial;
      }

      inline const RingElem* myGetPolPtr() const {
        return &myPolynomial;
      }
      /**
       * This function returns the current ancestor as PPMonoidElem
       */
      inline const PPMonoidElem& myGetAnc() const {
        return myAncestor;
      }

      /**
       * This function checks if a variable is already prolonged
       * @param index is the index of the variable which shall be checked
       * returns true if the variable is already prolonged
       */
      inline bool IamProlonged(const long& index) const {
        return myAlreadyProlongedVars[index];
      }

      /**
       * This function sets all variables to not already prolonged
       */
      void myClearProlongedVars();

    private:
      friend struct isLower;
      friend struct IsLowerLPP;
      friend class TQSets;

      /**
       * The polynomial of the Janet-Triple.
       */
      RingElem myPolynomial;

      /**
       * The ancestor of the polynomial.
       */
      PPMonoidElem myAncestor;

      /**
       * The already prolonged variables of the polynomial
       */
      std::vector<bool> myAlreadyProlongedVars;

    };

    typedef std::list<JanetTriple>::iterator nodeData;

    class JanetNodeBase;
    typedef std::shared_ptr<JanetNodeBase> JanetHandle;

    /****************************************************************************************************************************/
    /*
     * The JanetNode are the nodes in the JanetTree. It's the father
     * class of two different node types. In the father-class we store
     * the varDistance and the degDistance to the next node and
     * define some virtual functions for the child-classes (We have to
     * do this because we use a Handle-Class to store the JanetNodes in
     * std::lists.
     */
    class JanetNodeBase {
    public:

      /**
       * Virtual destructor (further informations: C++ Primer)
       */
      virtual ~JanetNodeBase() {
      }

      /**
       * This function returns the degree-distance to the next
       * degree-node
       */
      inline long myGetDisNextDeg() const {
        return myDisNextDeg;
      }

      /**
       * This function returns the variable-distance to the next
       * variable-node
       */
      inline long myGetDisNextVar() const {
        return myDisNextVar;
      }

      /**
       * This function sets the degree-distance to the next degree-node
       * @param dis is the new degree-distance
       */
      inline void mySetDisNextDeg(const long& dis) {
        myDisNextDeg = dis;
      }

      /**
       * This function sets the variable-distance to the next
       * variable-node
       * @param dis is the new variable-distance
       * maybe protected and a new virtual function so that JanetLeafNodeImpl can't
       * set another distance
       */
      inline void mySetDisNextVar(const long& dis) {
        myDisNextVar = dis;
      }

      /**
       * only a virtual function. You need this function in the class
       * JanetInternalNode
       */
      virtual inline std::list<JanetHandle>* myNextVarArm() {
        CoCoA_THROW_ERROR("You don't have the permisson to deal with JanetNodes, maybe there is a mistake in the datastructure",
                    "myNextVarArm");  //??remove in the final version??
        return nullptr;
      }

      /**
       * only a virtual function. You need this function in both child
       * classes
       */
      virtual inline bool IamLeafNode() {
        return false;
      }

      /**
       * only a virtual function. You need this function in the class
       * JanetInternalNode
       * Warnings because I don't use nextArm
       */
      virtual inline void mySetNextVarArm(const std::list<JanetHandle>& /*NextArm*/) {
        CoCoA_THROW_ERROR("You don't have the permisson to deal with JanetNodes, maybe there is a mistake in the datastructure",
                    "mySetNextVarArm");
      }

    protected:
      /**
       * The constructor of the JanetNode. The constructor is proteced
       * because the user shall not use this class
       * @param deg set the distance to the next degree-node
       * @param var set the distance to the next var-node
       * @param newTriple set the JanetTriple in the node
       */
      JanetNodeBase(const long& deg, const long& var)
          : myDisNextDeg(deg),
            myDisNextVar(var) {
      }

      /**
       * the degree-distance to the next degree-node in the
       * JanetTree. If the distance is zero there is no next
       * degree-node in the JanetTree
       */
      long myDisNextDeg;

      /**
       * the variable-distance to the next variable-node in the
       * JanetTree. If the distance is zero there is no next
       * variable-node in the JanetTree
       */
      long myDisNextVar;
    };

    /***********************************************************************************************************************************/
    /** \brief A child-class of JanetNode. In the JanetTree this class
     * represents the leaf nodes
     *
     * The class represents the leaf nodes in the JanetTree. The class
     * contains the variables disNextDeg and disNextVar from the
     * father-class JanetNode and a JanetTriple as nodeData
     */
    class JanetLeafNodeImpl : public JanetNodeBase {
    public:
      /**
       * The constructor of the JanetLeafNode
       * @param newTriple is the triple of the JanetLeafNode
       * disNextVar and disNextDeg is set to zero
       */
      JanetLeafNodeImpl(const nodeData newTriple)
          : JanetNodeBase(0, 0),
            myTriple(newTriple)  //, triple(newTriple)
      {
      }

      /**
       * This function clones the JanetNode. You need this function for
       * the class JanetHandle.
       */
      inline JanetLeafNodeImpl* myClone() const {
        return new JanetLeafNodeImpl(*this);
      }

      /**
       *This function returns a pointer on the JanetTriple.
       */
      inline JanetTriple* myGetTriplePtr() {
        return &(*myTriple);
      }

      inline JanetTriple& myGetTriple() {
        return *myTriple;
      }

      inline nodeData myGetTripleListPointer() {
        return myTriple;
      }

      /**
       * This function returns always true. It checks if the node is a
       * leafNode
       */
      inline bool IamLeafNode() {
        return true;
      }

      /**
       * only dummy function which we don't need in this class. We
       * return an empty arm
       */
      inline std::list<JanetHandle>* myNextVarArm() {
        CoCoA_THROW_ERROR(
            "You don't have the permisson to deal with arms in JanetLeafNodes, maybe there is a mistake in the datastructure",
            "myNextVarArm");  //??maybe remove this in the final version??
        return nullptr;
      }

      /**
       * only dummy function which we don't need in this class. If it
       * isn't a leaf node then it set's a new arm
       * !!WARNING because we don't need nextArm!!
       */
      inline void mySetNextVarArm(const std::list<JanetHandle>& /*NextArm*/) {
        CoCoA_THROW_ERROR(
            "You don't have the permisson to deal with arms JanetLeafNodes, maybe there is a mistake in the datastructure",
            "mySetNextVarArm");
      }

    private:
      /**
       *Stores the data for the leaf nodes. It is an iterator of a list
       *where we store the janet triples.
       */
      nodeData myTriple;
    };

    /***********************************************************************************************************************************/
    /** \brief A child-class of JanetNode. In the JanetTree this class
     * represent the internal nodes
     *
     * The class represents the internal nodes in the JanetTree. The
     * class contains the data of the janetsubtree in var direction
     */
    class JanetInternalNodeImpl : public JanetNodeBase {
    public:

      /**
       * the constructor of the class JanetInternalNode. It initializes
       * an empty arm of the JanetTree. The arm represents the
       * degree-direction in the next variable
       */
      JanetInternalNodeImpl(/*const SparsePolyRing& polyRing*/)
          : JanetNodeBase(0, 0),
            myArm() {
      }

      /**
       * This function clones the JanetInternalNode. You need this function for
       * the class JanetHandle.
       */
      inline JanetInternalNodeImpl* myClone() const {
        return new JanetInternalNodeImpl(*this);
      }

      /**
       * This function returns a pointer to the degree-arm in the next
       * variable
       */
      inline std::list<JanetHandle>* myNextVarArm() {
        return &myArm;
      }

      /**
       * This function checks if there is a node which is a leafNode and returns
       * false because it is an internal node
       */
      inline bool IamLeafNode() {
        return false;
      }

      /**
       * This function replaces the next degree-arm with a new one
       */
      inline void mySetNextVarArm(const std::list<JanetHandle>& NextArm) {
        myArm = NextArm;
      }

    private:

      /**
       * The list where we store the next nodes in var direction
       */
      std::list<JanetHandle> myArm;

    };


    /*************************************************************************************************************************************/
    /**
     * \brief The JanetTree.
     *
     * This class stores only the beginning of the tree. The rest of the
     * tree is stored in the JanetInternalNodes
     */

    class JanetTree {
    public:
      /**
       * The constructor. The constructor initializes a list with one Element. The element is a JanetInternalNode.
       * @param deg is the beginning degree of the tree
       * @param var is the beginning variable of the tree
       */
      JanetTree(const SparsePolyRing& PolyRing, const long& deg, const long& var)
          : myArm(),
            myBeginDeg(deg),
            myBeginVar(var),
            myPolyRing(PolyRing)
      {
        myArm.push_back(JanetHandle(new JanetInternalNodeImpl()));
      }

      JanetTree(const SparsePolyRing& PolyRing)
          : myArm(),
            myBeginDeg(0),
            myBeginVar(0),
            myPolyRing(PolyRing)
      {
        myArm.push_back(JanetHandle(new JanetInternalNodeImpl()));
      }

      /**
       * Returns the whole rootArm as std::list<JanetHandle>
       * JanetIterator needs this function
       */
      inline std::list<JanetHandle>* myGetRootArm() {
        return &myArm;
      }

      /**
       * Same as above, but readonly
       */
      inline std::list<JanetHandle>* myGetRootArm() const {
        return const_cast<std::list<JanetHandle> *>(&myArm);
      }

      /**
       * Returns the beginning degree of the tree
       */
      inline long myGetBeginDeg() const {
        return myBeginDeg;
      }

      /**
       *Returns the beginning variable of the tree
       */
      inline long myGetBeginVar() const {
        return myBeginVar;
      }

      /**
       *Returns the polyRing of the tree
       */
      inline const SparsePolyRing& myGetPolyRing() const {
        return myPolyRing;
      }

      /**
       *Deletes the content of the tree, but not the tree itself.
       *Attention: the polyRing remains!
       */
      inline void myDelete() {
        myArm.clear();
        myArm.push_back(JanetHandle(new JanetInternalNodeImpl()));
      }

      /**
       * myJInsertWithoutProlong insert the JanetTriple ps into the JanetTree JTree but doesn't perform any prolongations
       * @param JTree the JanetTree where we insert the triple
       * @param ps    the Itertor to the triple which we want to insert
       */
      void myInsert(const nodeData& ps);


      void myJTailNormalForm(RingElem& elem);

      void myJHeadNormalForm(RingElem& elem);

      void myJFullNormalForm(RingElem& elem);

      // TODO: Buggy!!!!!!!
      // TODO: Remove me, it is only used three times in easier usecases!
      /**
       * Delete the current JanetTree and add the new JanetTree
       * @param tree the new JanetTree
       * Only works if current Tree is empty!!!
       */
      void myAddAtBegin(JanetTree& tree);

      JanetTriple* myJDivisor(ConstRefPPMonoidElem w) const;

      /**
       * Check whether the PPMonidElem is at the tree
       */
      bool IamInternalNode(ConstRefPPMonoidElem w) const;

      /**
       * Insert a JanetTriple onto the tree i.e. a part of the tree will be deleted
       * It returns the iterators of the deleted nodes
       */
      std::vector<nodeData> myInsertInTreeStructure(const nodeData& triple);

      /**
       * Check whether the JanetTree is empty
       * @return true if the JanetTree is empty
       */
      inline bool IamEmpty() const
      {
        return myArm.empty();
      }

      /**
       * myOneLeafTree create a JanetTree from w (with starting variable i) to ps
       * @param  ps is an iterator to a JanetTriple
       * @param  i  is the starting variable
       * @param  w  is the starting monomial
       * @return returns the JanetTree which is beginning in m and ending in ps
       */
      static JanetTree OneLeafTree(const SparsePolyRing& polyRing, const std::list<JanetTriple>::iterator& ps, long i,
                                   PPMonoidElem w);
    private:
      /**
       * The first degree arm of the JanetTree
       */
      std::list<JanetHandle> myArm;

      /**
       * The beginning degree of the JanetTree
       */
      long myBeginDeg;

      /**
       * The beginning variable of the JanetTree
       */
      long myBeginVar;

      /**
       * The polyRing of the Tree
       */
      // const SparsePolyRing myPolyRing;
      SparsePolyRing myPolyRing;

      /**
       * only declaration, there isn't a definition because we can't
       * copy the tree
       * Nobody shall use this operator, therefore private
       */
      // JanetTree& operator=(const JanetTree &);
    };

    /*************************************************************************************************************************************/
    /**
     * \brief With an object of JanetIterator we can navigate through
     * the JanetTreẹ
     *
     * I try to programm this class similar to the STL (because of the
     * tree-structure I had to add some extra functionalities)
     *
     * After thinking about it: I guess, the only similarity is the name ;-)
     */
    class JanetIterator {
    public:

      /**
       * The constructor of the JanetIterator
       * @param tree is the JanetTree where the Iterator navigates through
       */
      JanetIterator(const JanetTree &tree)
          : myCurTree(&tree),
            myCurArm(tree.myGetRootArm()),
            myMonomial(NumIndets(tree.myGetPolyRing())),
            myCurVar(tree.myGetBeginVar()) {
        // myCurTree = &tree;
        myCurIter = myCurArm->begin();
        myMonomial[myCurVar] = tree.myGetBeginDeg();
      }

      /**
       * Changes the type of node
       * @param triple is the new Triple in the LeafNode
       *
       * Changes any kind of node in a leaf node.  This function deletes
       * all the nodes behind the current node
       */
      void myChangeToLeafNode(const nodeData triple);

      /**
       * This function sets the iterator on the next degree-node and
       * returns the degree-distance.  If it returns 0 then there is no
       * next node
       */
      long myNextDeg();

      /**
       * This function set the iterator on the next variable-node and
       * returns the variable-distance.  If it returns 0 then there is no
       * next node
       */
      long myNextVar();

      /**
       * This function returns the current degree of the current
       * variable on which the iterator points
       */
      inline long myCurrentDeg() const {
        return myMonomial[myCurVar];
      }

      /**
       * This function returns the current variable on which the iterator points
       */
      inline long myCurrentVar() const {
        return myCurVar;
      }

      /**
       * This function sets a new variable-node. This node is an
       * internal-node
       *
       * If there is another node in this direction, this node lies
       * behind the new node.  If the variable-distance of
       * the next node is lower or equal to dis then
       * this node will be deleted
       *
       * If the current node is a leaf Node we change this node
       * into an internal node
       *
       * @param dis is the variable distance between the current
       * and the new node
       */
      void mySetNextVar(long dis);

      /**
       * This function sets a new degree-node. This node is an
       * internal-node
       *
       * If there is another node in this direction then the new  node lies
       * between the current node and this node.  If the degree-distance of
       * the next node is lower or equal to dis
       * this node will be deleted and the new node will be inserted.
       *
       * ?? There is an assertion where we check if it isn't a leaf Node ??
       *
       * @param dis is the degree distance between the current node and
       * the new node
       */
      void mySetNextDeg(long dis);

      /**
       * This function connects two JanetTrees
       *
       * tree is appended on the current postion of the
       * iterator in deg-direction. The tree will be copied after this position in degree
       * direction.
       *
       * ?? WARNING: NO TESTS IF THE CONNECTION MAKES SENSE ??
       * Shall I warn?
       *
       * @param tree is the JanetTree which shall be connect with the current tree
       */
      void myConnectJanetTreeDeg(JanetTree& tree);

      /**
       * Checks if current Node is a leaf node
       */
      inline bool IamLeafNode() const {
        return (*myCurIter)->IamLeafNode();
      }

      /**
       * This function sets a nonMultiplicative variable
       *
       * This functions sets the varibale only if we are in a leaf
       * node. Then it returns true. If we are in a internal node
       * we do nothing and return false
       *
       * @param index is the index of the variable which we want to set
       */
      bool IamSetNM(const long &index);

      /**
       * This function returns the current polynomial of the current
       * JanetTriple as RingElem
       */
      inline const RingElem& myGetPol() const {
        CoCoA_ASSERT(IamLeafNode());
        JanetLeafNodeImpl* node(dynamic_cast<JanetLeafNodeImpl*>(myCurIter->get()));
        CoCoA_ASSERT(node != nullptr);

        return node->myGetTriple().myGetPol();
      }

      /**
       * This function returns the current ancestor of the current
       * JanetTriple as PPMonoidElem
       */
      inline PPMonoidElem myGetAnc() const {
        CoCoA_ASSERT(IamLeafNode());
        JanetLeafNodeImpl* node(dynamic_cast<JanetLeafNodeImpl*>(myCurIter->get()));
        CoCoA_ASSERT(node != nullptr);

        return node->myGetTriple().myGetAnc();
      }

      /**
       * This function checks if a variable is already prolonged
       *
       * If the variable is prolonged we return true
       * otherwise false
       *
       * @param index: The variable with the index which we test
       */
      inline bool IamProlonged(const long &index) const {
        CoCoA_ASSERT(IamLeafNode());
        JanetLeafNodeImpl* node(dynamic_cast<JanetLeafNodeImpl*>(myCurIter->get()));
        CoCoA_ASSERT(node != nullptr);

        return node->myGetTriple().IamProlonged(index);
      }

      /**
       * The iterator returns to the beginning of the JanetTree
       */
      void myReturnToBegin();

      /**
       * This function checks if current node is a leaf node
       */
      inline bool IamLeafNode(){
        return (*myCurIter)->IamLeafNode();
      }

      /**
       * This function returns the distance to the next degNode. If
       * there is no deg node the function returns zero
       */
      inline long myDisNextDeg() const {
        return (*myCurIter)->myGetDisNextDeg();
      }

      /**
       * This function returns the distance to the next varNode. If
       * there is no deg node the function returns zero
       */
      inline long myDisNextVar() const {
        return (*myCurIter)->myGetDisNextVar();
      }

      /**
       * The iterator jumps to the previous degNode. If there is no
       * previous degNode this function returns the current node
       */
      long myPrevDeg();

      /**
       *Returns the current monomial
       */
      inline PPMonoidElem myGetMonomial() const {
        return PPMonoidElem(PPM(myCurTree->myGetPolyRing()), myMonomial);
      }

      /**
       * Returns the highest Iterator starting from current node
       *
       */
      JanetIterator myGotoHighestNode();

      std::vector<nodeData> myAllNodesAboveMeIncludingMe() const;
    private:

      /**
       * A pointer to the tree on which the iterator works
       */
      const JanetTree* myCurTree;

      /**
       * A pointer on the arm on which the iterator is working at the
       * moment
       */
      std::list<JanetHandle>* myCurArm;

      /**
       * The iterator which points at the current node
       */
      std::list<JanetHandle>::iterator myCurIter;

      /**
       * The current monomial which we have at the current position
       */
      std::vector<long> myMonomial;

      /**
       * The current variable
       */
      long myCurVar;

    };

    /*
     * This class contains the Janet Tree and a list of JanetTriples
     * The Janet Tree contains iterators to this list
     * The class couples the Janet Tree to a corresponding list of JanetTriples
     */
    class JanetContainer {
    public:
      /*
       * Constructor: It takes the polynomial ring of the Janet Tree and initializes
       * a new tree with the same polynomial ring. After that it inserts all triples
       * from list. Beside the polynomial ring we do not use the JanetTree
       */
      JanetContainer(const JanetTree& tree, const std::list<JanetTriple>& list)
          : myTree(tree.myGetPolyRing())
      {
        myList.insert(myList.begin(), list.begin(), list.end());
        myInitializeTreeFromList();
      }

      /*
       * Constructor: It initializes a new tree with the polynomial ring. After
       * that it inserts all triples from list.
       */
      JanetContainer(const SparsePolyRing& ring, const std::list<JanetTriple>& list)
          : myTree(ring)
      {
        myList.insert(myList.begin(), list.begin(), list.end());
        myInitializeTreeFromList();
      }

      /*
       * Constructor: It initializes an empty tree with the given polynomial ring.
       */
      JanetContainer(const SparsePolyRing& ring)
          : myTree(ring)
      {
      }

      /*
       * Copy Constructor
       */
      JanetContainer(const JanetContainer& orig)
          : myTree(orig.myGetPolyRing())
      {
        myList.insert(myList.begin(), orig.myList.begin(), orig.myList.end());
        myInitializeTreeFromList();
      }

      /*
       * Returns an iterator to the begin of the internal list
       */
      std::list<JanetTriple>::iterator myListBegin()
      {
        return myList.begin();
      }

      /*
       * Returns a const_iterator to the begin of the internal list
       */
      std::list<JanetTriple>::const_iterator myListBegin() const
      {
        return myList.begin();
      }

      /*
       * Returns an iterator to the end of the internal list
       */
      std::list<JanetTriple>::iterator myListEnd()
      {
        return myList.end();
      }

      /*
       * Returns a const_iterator to the end of the internal list
      */
      std::list<JanetTriple>::const_iterator myListEnd() const
      {
        return myList.end();
      }

      /*
       * Returns a const reference to the internal JanetTree
       */
      const JanetTree& myGetTree() const
      {
        return myTree;
      }

      /*
       * It changes the basis to the list given by the iterators
       */
      void myChangeBasis(std::list<JanetTriple>::const_iterator newListBegin,
                         std::list<JanetTriple>::const_iterator newListEnd);

      /*
       * Returns the polynomial ring of the janet tree
       */
      const SparsePolyRing& myGetPolyRing() const
      {
        return myTree.myGetPolyRing();
      }
      /*
       * Checks if the container represents a constant ideal
       */
      inline bool IamConstant() const
      {
        return IsConstant(myList.front().myGetPol());
      }

    private:

      std::list<JanetTriple> myList;
      JanetTree myTree;

      /*
       * This method inserts all elements from myList into myTree.
       * Attention: It does not clear myTree before inserting.
       */
      void myInitializeTreeFromList();

    };
  } // end of namespace Involutive
}  // end of namespace CoCoA
#endif
