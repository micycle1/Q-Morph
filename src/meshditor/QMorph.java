package meshditor;

import java.util.ArrayList;
import java.util.List;

/**
 * This is the main class, implementing the triangle to quad conversion process.
 * The algorithm was invented by Steven J. Owen, Matthew L. Staten, Scott A.
 * Cannan, and Sunil Saigal and described in their paper "Advancing Front
 * Quadrilateral Meshing Using Triangle Transformations" (1998).
 * 
 * @see <a href="http://www.andrew.cmu.edu/user/sowen/abstracts/Ow509.html"
 *      TARGET="_top">the abstract of the paper</a>
 * @author Karl Erik Levik, karll@ifi.uio.no
 *
 */
public class QMorph extends GeomBasics {
	
	public QMorph() {
	}

	private List<Edge> frontList;
	private boolean finished = false;
	private int level = 0;
	private int nrOfFronts = 0;
	private boolean evenInitNrOfFronts = false;

	/** Initialize the class */
	public void init() {

		if (leftmost == null || lowermost == null || rightmost == null || uppermost == null) {
			findExtremeVertices();
		}

		if (doCleanUp) {
			topoCleanup = new TopoCleanup();
		}

		if (doSmooth) {
			globalSmooth = new GlobalSmooth();
		}

		if (doTri2QuadConversion) {
			setCurMethod(this);
			elementList = new ArrayList<>();
			finished = false;
			level = 0;

			repairZeroAreaTriangles();

			if (!verifyTriangleMesh(triangleList)) {
				Msg.error("This is not an all-triangle mesh.");
			}

			Edge.clearStateList();
			frontList = defineInitFronts(edgeList);
			printEdgeList(frontList);
			classifyStateOfAllFronts(frontList);

			nrOfFronts = countNOFrontsAtCurLowestLevel(frontList);
			if ((nrOfFronts & 1) == 0) {
				evenInitNrOfFronts = true;
			} else {
				evenInitNrOfFronts = false;
			}
		} else if (doCleanUp) {
			setCurMethod(topoCleanup);
			topoCleanup.init();
		} else if (doSmooth) {
			setCurMethod(globalSmooth);
			globalSmooth.init();
		}
	}

	/** Run the implementation on the given triangle mesh */
	public void run() {
		if (doTri2QuadConversion) {
			if (!step) {
				

				// The program's main loop from where all the real action originates
				while (!finished) {
					step();
				}

			}
		} else if (!step) {
			// Post-processing methods
			if (doCleanUp) {
				topoCleanup.init();
				topoCleanup.run();
			}

			if (doSmooth) {
				globalSmooth.init();
				globalSmooth.run();
			}

			printElements(elementList);
		}
	}

	/** Step through the morphing process one front edge at the time. */
	@Override
	public void step() {
		Quad q;
		Edge e;
		int i, oldBaseState;

		e = Edge.getNextFront(/* frontList, */);
		if (e != null) {
			oldBaseState = e.getState();
			if (nrOfFronts <= 0) {
				level++;
				nrOfFronts = countNOFrontsAtCurLowestLevel(frontList);
			}

			Edge.printStateLists();
			q = makeQuad(e);
			if (q == null) {
				while (q == null && e != null) {
					e.selectable = false;
					e = Edge.getNextFront();
					if (e == null) {
						break;
					}
					oldBaseState = e.getState();
					q = makeQuad(e);
				}
				Edge.markAllSelectable();
				Msg.warning("Main loop: makeQuad(..) returned null, so I chose another edge. It's alright now.");
			}
			if (q.firstVertex != null) {
				elementList.add(q); // is not really a quad, but an area that needs
									// local smoothing and update.
			}

			printEdgeList(frontList);
			Edge.printStateLists();

			if (q.firstVertex != null) {
				localSmooth(q, frontList);
			}
			i = localUpdateFronts(q, level, frontList);
			printEdgeList(frontList);
			Edge.printStateLists();
			nrOfFronts -= i;
		} else if (!finished) {
			// Post-processing methods
			if (doCleanUp) {
				topoCleanup.init();
				if (!step) {
					topoCleanup.run();
				} else {
					setCurMethod(topoCleanup);
					topoCleanup.step();
					return;
				}
			}
			if (doSmooth) {
				globalSmooth.init();
				globalSmooth.run();
			}

			printElements(elementList);
			finished = true;
		}
	}

	/** @return the frontList, that is, the list of front edges */
	public List<Edge> getFrontList() {
		return frontList;
	}

	/**
	 * Count the number of front edges in the new loops created if we were to create
	 * a new quad with edge b as base edge and one or both of edges l and r were
	 * promoted to front edges and one or both of their top vertices were located on
	 * the front. Return a byte signifying which, if any, of these loops contains an
	 * odd number of edges.
	 *
	 * @param b     the base edge of the quad we want to make. Should be a front
	 *              edge.
	 * @param l     the left edge of the quad we want to make
	 * @param r     the right edge of the quad we want to make
	 * @param lLoop boolean indicating whether to count the edges in the loop at l
	 * @param rLoop boolean indicating whether to count the edges in the loop at r
	 * @return a bit pattern indicating which, if any, of the resulting front loops
	 *         that have an odd number of edges. Also indicate if the two loops are
	 *         actually the same.
	 */
	private byte oddNOFEdgesInLoopsWithFEdge(Edge b, Edge l, Edge r, boolean lLoop, boolean rLoop) {
		byte ret = 0;
		int ln = 0, rn = 0;
		bothSidesInLoop = false;

		if (lLoop) {
			ln = countFrontsInNewLoopAt(b, l, r);
		}
		if (!bothSidesInLoop && rLoop) {
			rn = countFrontsInNewLoopAt(b, r, l);
		}

		if ((ln & 1) == 1) {
			ret += 1;
		}
		if ((rn & 1) == 1) {
			ret += 2;
		}
		if (bothSidesInLoop) {
			ret += 4;
		}

		return ret;
	}

	private boolean bothSidesInLoop = false;

	/**
	 * Supposing that the edges side and otherSide are promoted to front edges. The
	 * method parses a new loop involving the edges side, otherSide and possibly
	 * edge b. If two separate loops result from promoting the edges side and
	 * otherSide to front edges, then return the number of edges in the loop
	 * involving edges b and side. Otherwise, return the total number of edges in
	 * the loop involving edges side and otherSide.
	 *
	 * @param b         an edge that is considered for being promoted to a front
	 *                  edge
	 * @param side      one of b's front neighbors
	 * @param otherSide the other of b's front neighbors
	 * @return as described above.
	 */
	private int countFrontsInNewLoopAt(Edge b, Edge side, Edge otherSide) {
		Edge cur = null, prev, tmp;
		Vertex n1 = side.commonVertex(b), n2 = otherSide.commonVertex(b), n3 = side.otherVertex(n1), n4 = otherSide.otherVertex(n2);
		int count1stLoop = 1, count2ndLoop = 1, count3rdLoop = 1, n3n4Edges = 0, count = 0;
		boolean n4Inn3Loop = false;

		// First find the front edge where we want to start
		prev = b;
		cur = prev.frontNeighborAt(n1);
		

		// Parse the current loop of fronts:
		do {
			tmp = cur.nextFrontNeighbor(prev);
			prev = cur;
			cur = tmp;
			count1stLoop++;
			if (cur.hasVertex(n1)) { // Case 2,3,4,5 or 6
				bothSidesInLoop = true;
				prev = n3.anotherFrontEdge(null);
				cur = prev.frontNeighborAt(n3);

				double alpha = cur.computeCCWAngle(side);
				double beta = prev.computeCCWAngle(side);

				// Should we go the other way?
				if (beta < alpha) {
					tmp = cur;
					cur = prev;
					prev = tmp;
				}
				if (cur.hasVertex(n4)) {
					n4Inn3Loop = true;
					n3n4Edges = 1;
				}
				

				// Parse the n3 loop and count number of edges here
				do {
					tmp = cur.nextFrontNeighbor(prev);
					prev = cur;
					cur = tmp;
					count2ndLoop++;
					if (!n4Inn3Loop && cur.hasVertex(n4)) {
						n4Inn3Loop = true;
						n3n4Edges = count2ndLoop;
					}
				} while (!cur.hasVertex(n3) /* && count2ndLoop < 300 */);

				if (!n4Inn3Loop && !otherSide.isFrontEdge() && n4.frontVertex()) { // Case 5 only
					prev = n4.anotherFrontEdge(null);
					cur = prev.frontNeighborAt(n4);
					

					// Parse the n4 loop and count number of edges here
					do {
						tmp = cur.nextFrontNeighbor(prev);
						prev = cur;
						cur = tmp;
						count3rdLoop++;
					} while (!cur.hasVertex(n4) /* && count3rdLoop < 300 */);
				}

				if (n4Inn3Loop && !otherSide.isFrontEdge()) {// Case 2,3
					count = count1stLoop + count2ndLoop + 1 - n3n4Edges;
				} else if (!n4Inn3Loop && otherSide.isFrontEdge()) {// Case 4
					count = count1stLoop + count2ndLoop;
				} else if (!n4Inn3Loop && !otherSide.isFrontEdge() && !n4.frontVertex()) {// Case 6
					count = count1stLoop + count2ndLoop + 2;
				} else if (!n4Inn3Loop && !otherSide.isFrontEdge() && n4.frontVertex()) {// Case 5
					count = count1stLoop + count2ndLoop + count3rdLoop + 2;
				}
				break;
			}
		} while (!cur.hasVertex(n1) && !cur.hasVertex(n3) /* && count1stLoop < 300 */);
		/*
		 * if (count1stLoop>= 300) Msg.
		 * error("QMorph.countFrontsInNewLoopAt(..): Suspiciously many fronts in loop!"
		 * );
		 */

		if (!bothSidesInLoop) {
			count = count1stLoop + 1; // Add the new edge between n1 and n3
		}

		return count;
	}

	/** Returns the number of front edges at the currently lowest level loop(s). */
	private int countNOFrontsAtCurLowestLevel(List<Edge> frontList) {
		Edge cur;
		int lowestLevel, count = 0;

		// Get nr of fronts at the lowest level:
		if (frontList.size() > 0) {
			cur = (Edge) frontList.get(0);
			lowestLevel = cur.level;
			count = 1;
			for (int i = 1; i < frontList.size(); i++) {
				cur = (Edge) frontList.get(i);
				if (cur.level < lowestLevel) {
					lowestLevel = cur.level;
					count = 1;
				} else if (cur.level == lowestLevel) {
					count++;
				}
			}
		}
		return count;
	}

	/** Make sure the triangle mesh consists exclusively of triangles */
	private boolean verifyTriangleMesh(List<Triangle> triangleList) {
		Object o;
		for (Object element : triangleList) {
			o = element;
			if (!(o instanceof Triangle)) {
				return false;
			}
		}
		return true;
	}

	/**
	 * Construct a new quad
	 * 
	 * @param e the base edge of this new quad
	 * @return the new quad
	 */
	private Quad makeQuad(Edge e) {
		Quad q;
		Triangle t;
		Edge top;
		int index;
		boolean initL = false, initR = false;
		Edge[] sideEdges = new Edge[2];
		Edge leftSide, rightSide;
		ArrayList<Triangle> tris;

		if (e.leftSide) {
			initL = true;
		}
		if (e.rightSide) {
			initR = true;
		}

		q = handleSpecialCases(e);
		if (q != null) {
			return q;
		}
		sideEdges = makeSideEdges(e);
		leftSide = sideEdges[0];
		rightSide = sideEdges[1];
		if (leftSide.commonVertex(rightSide) != null) {
			Msg.warning("Cannot create quad, so we settle with a triangle:");

			t = new Triangle(e, leftSide, rightSide);
			index = triangleList.indexOf(t);
			if (index != -1) {
				t = (Triangle) triangleList.get(index);
			}

			q = new Quad(t);
			clearQuad(q, t.edgeList[0].getTriangleElement());

			return q;
		} else {
			Vertex nC = leftSide.otherVertex(e.leftVertex);
			Vertex nD = rightSide.otherVertex(e.rightVertex);
			top = recoverEdge(nC, nD);

			// We need some way to handle errors from recoverEdge(..). This is solved by
			// a recursive approach, where makeQuad(e) is called with new parameter
			// values for e until it returns a valid Edge.
			if (top == null) {
				e.removeFromStateList();
				e.classifyStateOfFrontEdge();

				if (leftSide.frontEdge) {
					leftSide.removeFromStateList();
					leftSide.classifyStateOfFrontEdge();
				}
				if (rightSide.frontEdge) {
					rightSide.removeFromStateList();
					rightSide.classifyStateOfFrontEdge();
				}

				return null;
			}
			q = new Quad(e, leftSide, rightSide, top);

			// Quad.trianglesContained(..) and clearQuad(..) needs one of q's interior
			// triangles as parameter. The base edge borders to only one triangle, and
			// this lies inside of q, so it is safe to choose this:

			tris = q.trianglesContained(e.getTriangleElement());
			if (q.containsHole(tris)) {
				e.removeFromStateList();
				e.classifyStateOfFrontEdge();

				if (leftSide.frontEdge) {
					leftSide.removeFromStateList();
					leftSide.classifyStateOfFrontEdge();
				}
				if (rightSide.frontEdge) {
					rightSide.removeFromStateList();
					rightSide.classifyStateOfFrontEdge();
				}

				return null;
			}

			clearQuad(q, tris/* e.getTriangleElement() */);
			// preSmoothUpdateFronts(q, frontList);
			return q;
		}
	}

	/**
	 * @param nK     front Vertex to be smoothed
	 * @param nJ     Vertex behind the front, connected to nK
	 * @param myQ    an arbitrary selected quad connected to Vertex nK
	 * @param front1 a front edge adjacent nK
	 * @param front2 another front edge adjacent nK and part of the same loop as
	 *               front1
	 * @return a new Vertex with the smoothed position
	 */
	private Vertex smoothFrontVertex(Vertex nK, Vertex nJ, Quad myQ, Edge front1, Edge front2) {
		List<Element> adjQuads = nK.adjQuads();
		double tr, ld = 0;
		Quad q;
		Vertex newVertex;
		int n = 0, seqQuads = myQ.nrOfQuadsSharingAnEdgeAt(nK) + 1;

		if (front1 == null) {
			Msg.error("front1 is null");
		}
		if (front2 == null) {
			Msg.error("front2 is null");
		}

		Edge eD = new Edge(nK, nJ);
		if (nK.edgeList.contains(eD)) {
			eD = (Edge) nK.edgeList.get(nK.edgeList.indexOf(eD));
		}

		if (seqQuads == 2) { // || nK.nrOfFrontEdges()> 2) {
			if (front1.length() > front2.length()) {
				tr = front1.length() / front2.length();
			} else {
				tr = front2.length() / front1.length();
			}

			if (tr > 20) {
				// Msg.debug("******************* tr > 20");
				ld = nK.meanNeighborEdgeLength();
				newVertex = eD.otherVertexGivenNewLength(ld, nJ);
			} else if (tr <= 20 && tr > 2.5) {
				
				// The mean of (some of) the other edges in all the adjacent elements
				ld = 0;
				// First add the lengths of edges from Triangles ahead of the front
				for (Edge e : nK.edgeList) {
					if (e != eD && e != front1 && e != front2) {
						ld += e.length();
						n++;
					}
				}
				// Then add the lengths of edges from the two Quads behind the front
				q = front1.getQuadElement();
				ld += q.edgeList[base].length();
				if (q.edgeList[left] != eD) {
					ld += q.edgeList[left].length();
				} else {
					ld += q.edgeList[right].length();
				}

				q = front2.getQuadElement();
				ld += q.edgeList[base].length();
				if (q.edgeList[left] != eD) {
					ld += q.edgeList[left].length();
				} else {
					ld += q.edgeList[right].length();
				}

				ld = ld / (4.0 + n);
				newVertex = eD.otherVertexGivenNewLength(ld, nJ);
			} else {
				// Msg.debug("******************* tr<= 2.5");
				newVertex = nK.blackerSmooth(nJ, front1, front2, eD.length());
			}
		} else {
			newVertex = nK.blackerSmooth(nJ, front1, front2, eD.length());
		}

		return newVertex;
	}

	/**
	 * Calculate smoothed position of Vertex. (Called from localSmooth)
	 * 
	 * @see #smoothFrontVertex(Vertex, Vertex, Quad, Edge, Edge)
	 * @see Vertex#modifiedLWLaplacianSmooth()
	 * @param n Vertex to be smoothed
	 * @param q Quad to which n belongs
	 */
	Vertex getSmoothedPos(Vertex n, Quad q) {
		Vertex newN, behind;
		Edge front1, front2, e;
		Quad q2, qn;
		Element neighbor;

		if (q.hasFrontEdgeAt(n) && !n.boundaryVertex()) {
			if (n == q.edgeList[left].otherVertex(q.edgeList[base].leftVertex)) {
				if (q.edgeList[top].isFrontEdge()) {
					front1 = q.edgeList[top];
				} else {
					front1 = q.edgeList[left];
				}
				front2 = front1.trueFrontNeighborAt(n);
				q2 = front2.getQuadElement();
				e = q.commonEdgeAt(n, q2);
				if (e == null) {
					if (q.edgeList[top] != front1 && q.edgeList[top] != front2) {
						e = q.edgeList[top];
					} else {
						e = q.edgeList[left];
					}
				}
			} else if (n == q.edgeList[right].otherVertex(q.edgeList[base].rightVertex)) {
				if (q.edgeList[top].isFrontEdge()) {
					front1 = q.edgeList[top];
				} else {
					front1 = q.edgeList[right];
				}
				front2 = front1.trueFrontNeighborAt(n);
				q2 = front2.getQuadElement();
				e = q.commonEdgeAt(n, q2);
				if (e == null) {
					if (q.edgeList[top] != front1 && q.edgeList[top] != front2) {
						e = q.edgeList[top];
					} else {
						e = q.edgeList[right];
					}
				}
			} else if (n == q.edgeList[base].leftVertex) {
				if (q.edgeList[left].isFrontEdge()) {
					front1 = q.edgeList[left];
				} else {
					front1 = q.edgeList[base];
				}
				front2 = front1.trueFrontNeighborAt(n);
				q2 = front2.getQuadElement();
				e = q.commonEdgeAt(n, q2);
				if (e == null) {
					if (q.edgeList[left] != front1 && q.edgeList[left] != front2) {
						e = q.edgeList[left];
					} else {
						e = q.edgeList[base];
					}
				}
			} else if (n == q.edgeList[base].rightVertex) {
				if (q.edgeList[right].isFrontEdge()) {
					front1 = q.edgeList[right];
				} else {
					front1 = q.edgeList[base];
				}
				front2 = front1.trueFrontNeighborAt(n);
				q2 = front2.getQuadElement();
				e = q.commonEdgeAt(n, q2);
				if (e == null) {
					if (q.edgeList[right] != front1 && q.edgeList[right] != front2) {
						e = q.edgeList[right];
					} else {
						e = q.edgeList[base];
					}
				}
			} else {
				Msg.error("getSmoothedPos(..): Vertex n is not in Quad q");
				newN = null;
				front1 = null;
				front2 = null;
				e = null;
			}
			behind = e.otherVertex(n);
			newN = smoothFrontVertex(n, behind, q, front1, front2);
		} else if (!n.boundaryVertex()) {
			newN = n.modifiedLWLaplacianSmooth();
		} else {
			newN = n;
		}
		return newN;
	}

	/**
	 * Smooth as explained in Owen's paper Each Vertex in the newly formed quad is
	 * smoothed. So is every Vertex directly connected to these.
	 */
	private void localSmooth(Quad q, List<Edge> frontList) {
		Quad tempQ1, tempQ2;
		Vertex n, nNew, nOld;
		Edge e, fe1, fe2;
		Vertex bottomLeft, bottomRight, bottomLeftNew, bottomRightNew, bottomLeftOld, bottomRightOld;
		List<Vertex> adjVertices;
		List<Vertex> adjVerticesNew;

		if (q.isFake) {

			Vertex top = q.edgeList[left].otherVertex(q.edgeList[base].leftVertex);
			bottomLeft = q.edgeList[base].leftVertex;
			bottomRight = q.edgeList[base].rightVertex;

			Vertex topNew;

			Vertex topOld = top.copyXY();
			bottomLeftOld = bottomLeft.copyXY();
			bottomRightOld = bottomRight.copyXY();

			adjVertices = q.getAdjVertices();
			adjVerticesNew = new ArrayList<>(adjVertices.size());

			// Calculate smoothed pos for each element Vertex and those vertices connected to
			// the element. If the element has become inverted, then repair it.

			topNew = getSmoothedPos(top, q);
			bottomLeftNew = getSmoothedPos(bottomLeft, q);
			bottomRightNew = getSmoothedPos(bottomRight, q);

			for (int i = 0; i < adjVertices.size(); i++) {
				n = (Vertex) adjVertices.get(i);
				if (n.frontVertex() && !n.boundaryVertex()) {
					fe1 = n.anotherFrontEdge(null);
					fe2 = n.anotherFrontEdge(fe1);
					tempQ1 = fe1.getQuadElement();
					tempQ2 = fe2.getQuadElement();
					e = tempQ1.commonEdgeAt(n, tempQ2);
					if (e == null) {
						e = tempQ1.neighborEdge(n, fe1);
					}
					adjVerticesNew.add(smoothFrontVertex(n, e.otherVertex(n), tempQ1, fe1, fe2));
				} else if (!n.boundaryVertex()) {
					adjVerticesNew.add(n.modifiedLWLaplacianSmooth());
				} else {
					adjVerticesNew.add(n);
				}
			}

			if (!top.equals(topNew)) {
				top.moveTo(topNew);
				inversionCheckAndRepair(top, topOld);
				top.update();
			}
			// else
			// top.updateAngles();

			if (!bottomLeft.equals(bottomLeftNew)) {
				bottomLeft.moveTo(bottomLeftNew);
				inversionCheckAndRepair(bottomLeft, bottomLeftOld);
				bottomLeft.update();
			}
			// else
			// bottomLeft.updateAngles();

			if (!bottomRight.equals(bottomRightNew)) {
				bottomRight.moveTo(bottomRightNew);
				inversionCheckAndRepair(bottomRight, bottomRightOld);
				bottomRight.update();
			}
			// else
			// bottomRight.updateAngles();

			for (int i = 0; i < adjVertices.size(); i++) {
				n = (Vertex) adjVertices.get(i);
				nOld = n.copyXY();
				nNew = (Vertex) adjVerticesNew.get(i);
				if (!n.equals(nNew)) {
					n.moveTo(nNew);
					inversionCheckAndRepair(n, nOld);
					n.update();
				}
				// else
				// n.updateAngles();
			}

			return;
		} else {
			Vertex topLeft = q.edgeList[left].otherVertex(q.edgeList[base].leftVertex);
			Vertex topRight = q.edgeList[right].otherVertex(q.edgeList[base].rightVertex);
			bottomLeft = q.edgeList[base].leftVertex;
			bottomRight = q.edgeList[base].rightVertex;

			Vertex topLeftNew, topRightNew;

			Vertex topLeftOld = topLeft.copyXY();
			Vertex topRightOld = topRight.copyXY();
			bottomLeftOld = bottomLeft.copyXY();
			bottomRightOld = bottomRight.copyXY();

			adjVertices = q.getAdjVertices();
			adjVerticesNew = new ArrayList<Vertex>(adjVertices.size());

			// Calculate smoothed pos for each element Vertex and those vertices connected to
			// the element. If the element has become inverted, then repair it.

			topLeftNew = getSmoothedPos(topLeft, q);
			topRightNew = getSmoothedPos(topRight, q);
			bottomLeftNew = getSmoothedPos(bottomLeft, q);
			bottomRightNew = getSmoothedPos(bottomRight, q);

			for (int i = 0; i < adjVertices.size(); i++) {
				n = (Vertex) adjVertices.get(i);
				if (n.frontVertex() && !n.boundaryVertex()) {
					fe1 = n.anotherFrontEdge(null);
					fe2 = n.anotherFrontEdge(fe1);
					tempQ1 = fe1.getQuadElement();
					tempQ2 = fe2.getQuadElement();
					e = tempQ1.commonEdgeAt(n, tempQ2);
					if (e == null) {
						e = tempQ1.neighborEdge(n, fe1);
					}
					adjVerticesNew.add(smoothFrontVertex(n, e.otherVertex(n), tempQ1, fe1, fe2));
				} else if (!n.boundaryVertex()) {
					adjVerticesNew.add(n.modifiedLWLaplacianSmooth());
				} else {
					adjVerticesNew.add(n);
				}
			}

			if (!topLeft.equals(topLeftNew)) {
				topLeft.moveTo(topLeftNew);
				inversionCheckAndRepair(topLeft, topLeftOld);
				topLeft.update();
			}
			// else
			// topLeft.updateAngles();

			if (!topRight.equals(topRightNew)) {
				topRight.moveTo(topRightNew);
				inversionCheckAndRepair(topRight, topRightOld);
				topRight.update();
			}
			// else
			// topRight.updateAngles();

			if (!bottomLeft.equals(bottomLeftNew)) {
				bottomLeft.moveTo(bottomLeftNew);
				inversionCheckAndRepair(bottomLeft, bottomLeftOld);
				bottomLeft.update();
			}
			// else
			// bottomLeft.updateAngles();

			if (!bottomRight.equals(bottomRightNew)) {
				bottomRight.moveTo(bottomRightNew);
				inversionCheckAndRepair(bottomRight, bottomRightOld);
				bottomRight.update();
			}
			// else
			// bottomRight.updateAngles();

			for (int i = 0; i < adjVertices.size(); i++) {
				n = (Vertex) adjVertices.get(i);
				nOld = n.copyXY();
				nNew = (Vertex) adjVerticesNew.get(i);
				if (!n.equals(nNew)) {
					n.moveTo(nNew);
					inversionCheckAndRepair(n, nOld);
					n.update();
				}
				// else
				// n.updateAngles();
			}
		}
	}

	/**
	 * Delete all interior triangles within the edges of this quad
	 * 
	 * @param q    the quad
	 * @param tris the list of triangles to be deleted
	 */
	private void clearQuad(Quad q, ArrayList<Triangle> tris) {
		int VertexInd, edgeInd, triInd;
		Vertex vertex;
		Edge e;

		for (int j = 0; j < tris.size(); j++) {
			Triangle t = tris.get(j);
			for (int i = 0; i < 3; i++) {
				e = t.edgeList[i];
				if (!q.hasEdge(e)) {
					edgeInd = edgeList.indexOf(e);
					if (edgeInd != -1) {
						edgeList.remove(edgeInd);
					}

					e.tryToDisconnectVertices();

					// Remove the leftVertex if not a vertex of q:
					if (!q.hasVertex(e.leftVertex)) {
						VertexInd = vertexList.indexOf(e.leftVertex);
						if (VertexInd != -1) {
							vertexList.remove(VertexInd);
						}
					}

					// Remove the rightVertex if not a vertex of q:
					if (!q.hasVertex(e.rightVertex)) {
						VertexInd = vertexList.indexOf(e.rightVertex);
						if (VertexInd != -1) {
							vertexList.remove(VertexInd);
						}
					}

					e.disconnectFromElement(t);
				} else {
					e.disconnectFromElement(t);
				}
			}

			triInd = triangleList.indexOf(t);
			if (triInd != -1) {
				triangleList.remove(triInd);
			}
		}

		q.connectEdges();
	}

	/**
	 * "Virus" that removes all triangles and their edges and vertices inside of this
	 * quad. Assumes that only triangles are present, not quads, inside of q
	 */
	private void clearQuad(Quad q, Triangle first) {
		Element neighbor;
		Element cur; // triangle
		ArrayList<Element> n = new ArrayList<>();
		Edge e;
		int VertexInd, edgeInd, triInd;

		n.add(first);
		for (int j = 0; j < n.size(); j++) {
			cur = n.get(j);
			for (int i = 0; i < 3; i++) {
				e = cur.edgeList[i];
				if (!q.hasEdge(e)) {

					neighbor = cur.neighbor(e);
					if (neighbor != null && !n.contains(neighbor)) {
						n.add(neighbor);
					}

					edgeInd = edgeList.indexOf(e);
					if (edgeInd != -1) {
						edgeList.remove(edgeInd);
					}

					e.tryToDisconnectVertices();

					// Remove the leftVertex if not a vertex of q:
					if (!q.hasVertex(e.leftVertex)) {
						VertexInd = vertexList.indexOf(e.leftVertex);
						if (VertexInd != -1) {
							vertexList.remove(VertexInd);
						}
					}
					// Remove the rightVertex if not a vertex of q:
					if (!q.hasVertex(e.rightVertex)) {
						VertexInd = vertexList.indexOf(e.rightVertex);
						if (VertexInd != -1) {
							vertexList.remove(VertexInd);
						}
					}

					e.disconnectFromElement(cur);
				} else {
					e.disconnectFromElement(cur);
				}
			}

			triInd = triangleList.indexOf(cur);
			if (triInd != -1) {
				triangleList.remove(triInd);
			}

		}
		q.connectEdges();
	}

	/** Updates fronts in fake quads (which are triangles, really) */
	private int localFakeUpdateFronts(Quad q, int lowestLevel, List<Edge> frontList) {
		int curLevelEdgesRemoved = 0;
		Edge e;
		Element neighbor;
		List<Edge> needsNewFN = new ArrayList<Edge>();
		List<Edge> lostFNList = new ArrayList<Edge>();
		List<Edge> needsReclassification = new ArrayList<Edge>();

		Edge.printStateLists();

		// Decide whether the base edge belongs in the frontlist, and act accordingly
		e = q.edgeList[base];
		if (!e.frontEdge && e.isFrontEdge()) {
			e.promoteToFront(q.edgeList[base].level + 1, frontList);
			needsNewFN.add(e);
		} else if (e.frontEdge && !e.isFrontEdge()) {
			if (e.removeFromFront(frontList) && e.level == lowestLevel) {
				curLevelEdgesRemoved++;
			}
			e.removeFromStateList();
			if (e.leftFrontNeighbor != null && e.leftFrontNeighbor.isFrontEdge()) {
				lostFNList.add(e.leftFrontNeighbor);
			}
			if (e.rightFrontNeighbor != null && e.rightFrontNeighbor.isFrontEdge()) {
				lostFNList.add(e.rightFrontNeighbor);
			}
		}

		// Decide whether the left edge belongs in the frontlist, and act accordingly
		e = q.edgeList[left];
		if (!e.frontEdge && e.isFrontEdge()) {
			e.promoteToFront(q.edgeList[base].level + 1, frontList);
			needsNewFN.add(e);
		} else if (e.frontEdge && !e.isFrontEdge()) {
			if (e.removeFromFront(frontList) && e.level == lowestLevel) {
				curLevelEdgesRemoved++;
			}
			e.removeFromStateList();
			if (e.leftFrontNeighbor != null && e.leftFrontNeighbor.isFrontEdge()) {
				lostFNList.add(e.leftFrontNeighbor);
			}
			if (e.rightFrontNeighbor != null && e.rightFrontNeighbor.isFrontEdge()) {
				lostFNList.add(e.rightFrontNeighbor);
			}
		}

		Edge.printStateLists();

		// Decide whether the right edge belongs in the frontlist, and act accordingly
		e = q.edgeList[right];
		if (!e.frontEdge && e.isFrontEdge()) {
			e.promoteToFront(q.edgeList[base].level + 1, frontList);
			needsNewFN.add(e);
		} else if (e.frontEdge && !e.isFrontEdge()) {
			if (e.removeFromFront(frontList) && e.level == lowestLevel) {
				curLevelEdgesRemoved++;
			}
			Edge.printStateLists();
			if (e.leftFrontNeighbor != null && e.leftFrontNeighbor.isFrontEdge()) {
				lostFNList.add(e.leftFrontNeighbor);
			}
			if (e.rightFrontNeighbor != null && e.rightFrontNeighbor.isFrontEdge()) {
				lostFNList.add(e.rightFrontNeighbor);
			}
		}

		// Set new front neighbors when necessary:
		for (int i = 0; i < needsNewFN.size(); i++) {
			e = (Edge) needsNewFN.get(i);
			if (e.frontEdge) {
				if (e.setFrontNeighbors(frontList)) {
					if (e.hasFalseFrontNeighbor()) {
						needsNewFN.add(e);
					} else {
						needsReclassification.add(e);
					}
				}
			} else {
				Msg.error("ERROR");
			}
		}

		// Run through the list of Edges that have lost their old frontNeighbor and see
		// if they need to have their frontNeighbor set again.
		Edge l, r;
		for (int i = 0; i < lostFNList.size(); i++) {
			e = (Edge) lostFNList.get(i);

			if (e.isFrontEdge() && e.hasFalseFrontNeighbor()) {
				if (e.setFrontNeighbors(frontList)) {
					if (e.hasFalseFrontNeighbor()) {
						lostFNList.add(e);
					} else {
						needsReclassification.add(e);
					}
				}
			}
		}

		// Reclassify front edges when necessary:
		for (Object element : needsReclassification) {
			e = (Edge) element;
			e.removeFromStateList();
			e.classifyStateOfFrontEdge();
		}

		return curLevelEdgesRemoved;
	}

	// Do some neccessary updating of the fronts before localSmooth(..) is run
	private void preSmoothUpdateFronts(Quad q, List<Edge> frontList) {
		q.edgeList[top].setFrontNeighbors(frontList);
	}

	/**
	 * Define new fronts, remove old ones. Set new frontNeighbors. Reclassify front
	 * edges.
	 * 
	 * @return nr of edges removed that belonged to the currently lowest level.
	 */
	private int localUpdateFronts(Quad q, int lowestLevel, List<Edge> frontList) {
		if (q.isFake) {
			return localFakeUpdateFronts(q, lowestLevel, frontList);
		} else {
			int curLevelEdgesRemoved = 0;
			Edge e;
			List<Edge> lostFNList = new ArrayList<>();
			List<Edge> needsNewFN = new ArrayList<>();
			List<Edge> needsReclassification = new ArrayList<>();

			printEdgeList(q.edgeList[top].rightVertex.frontEdgeList());

			// Remove the base edge from frontList and from one of stateList[0..2]
			e = q.edgeList[base];
			if (e.frontEdge && !e.isFrontEdge()) {
				if (e.removeFromFront(frontList) && e.level == lowestLevel) {
					curLevelEdgesRemoved++;
				}
				e.removeFromStateList();
			}
			if (e.leftFrontNeighbor != null && e.leftFrontNeighbor.isFrontEdge()) {
				lostFNList.add(e.leftFrontNeighbor);
			}
			if (e.rightFrontNeighbor != null && e.rightFrontNeighbor.isFrontEdge()) {
				lostFNList.add(e.rightFrontNeighbor);
			}

			// Decide whether the top edge belongs in the frontlist, and act accordingly
			e = q.edgeList[top];
			if (e.getTriangleElement() != null) {
			}
			if (!e.frontEdge && e.isFrontEdge()) {
				e.promoteToFront(q.edgeList[base].level + 1, frontList);
				needsNewFN.add(e);
			} else if (e.frontEdge && !e.isFrontEdge()) {
				if (e.removeFromFront(frontList) && e.level == lowestLevel) {
					curLevelEdgesRemoved++;
				}
				e.removeFromStateList();
				if (e.leftFrontNeighbor != null && e.leftFrontNeighbor.isFrontEdge()) {
					lostFNList.add(e.leftFrontNeighbor);
				}
				if (e.rightFrontNeighbor != null && e.rightFrontNeighbor.isFrontEdge()) {
					lostFNList.add(e.rightFrontNeighbor);
				}
			}

			// Decide whether the left edge belongs in the frontlist, and act accordingly
			e = q.edgeList[left];
			if (e.getTriangleElement() != null) {
			}
			if (!e.frontEdge && e.isFrontEdge()) {
				e.promoteToFront(q.edgeList[base].level + 1, frontList);
				needsNewFN.add(e);
			} else if (e.frontEdge && !e.isFrontEdge()) {
				if (e.removeFromFront(frontList) && e.level == lowestLevel) {
					curLevelEdgesRemoved++;
				}
				e.removeFromStateList();
				if (e.leftFrontNeighbor != null && e.leftFrontNeighbor.isFrontEdge()) {
					lostFNList.add(e.leftFrontNeighbor);
				}
				if (e.rightFrontNeighbor != null && e.rightFrontNeighbor.isFrontEdge()) {
					lostFNList.add(e.rightFrontNeighbor);
				}
			}

			// Decide whether the right edge belongs in the frontlist, and act accordingly
			e = q.edgeList[right];
			if (e.getTriangleElement() != null) {
			}

			if (!e.frontEdge && e.isFrontEdge()) {
				e.promoteToFront(q.edgeList[base].level + 1, frontList);
				needsNewFN.add(e);
			} else if (e.frontEdge && !e.isFrontEdge()) {
				if (e.removeFromFront(frontList) && e.level == lowestLevel) {
					curLevelEdgesRemoved++;
				}
				e.removeFromStateList();
				if (e.leftFrontNeighbor != null && e.leftFrontNeighbor.isFrontEdge()) {
					lostFNList.add(e.leftFrontNeighbor);
				}
				if (e.rightFrontNeighbor != null && e.rightFrontNeighbor.isFrontEdge()) {
					lostFNList.add(e.rightFrontNeighbor);
				}
			}

			// Set new front neighbors when necessary:
			for (int i = 0; i < needsNewFN.size(); i++) {
				e = (Edge) needsNewFN.get(i);
				if (e.setFrontNeighbors(frontList)) {
					if (e.hasFalseFrontNeighbor()) {
						needsNewFN.add(e);
					} else {
						needsReclassification.add(e);
					}
				}
			}

			// Parse the list of Edges that have lost their old frontNeighbor and see
			// if they need to have their frontNeighbor set again.
			Edge l, r;
			for (int i = 0; i < lostFNList.size(); i++) {
				e = (Edge) lostFNList.get(i);

				if (e.isFrontEdge() && e.hasFalseFrontNeighbor()) {
					if (e.setFrontNeighbors(frontList)) {
						if (e.hasFalseFrontNeighbor()) {
							lostFNList.add(e);
						} else {
							needsReclassification.add(e);
						}
					}
				}
			}

			// Reclassify front edges when necessary:
			for (int i = 0; i < needsReclassification.size(); i++) {
				e = (Edge) needsReclassification.get(i);

				if (needsNewFN.contains(e)) {
				}
				if (lostFNList.contains(e)) {
				}

				e.removeFromStateList();

				e.classifyStateOfFrontEdge();
			}

			// NOTE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			// When Quad::smooth are implemented, I need also
			// to update *ALL* edges affected by that smoothing

			return curLevelEdgesRemoved;
		}
	}

	private List<Edge> defineInitFronts(List<Edge> edgeList) {
		List<Edge> frontList = new ArrayList<>();
		Edge leftEdge, rightEdge;
		for (Edge e : edgeList) {
			if (e.hasElement(null)) {
				e.promoteToFront(0, frontList);
				e.swappable = false;
			}
		}

		for (Edge e : frontList) {
			e.setFrontNeighbors(frontList);
		}

		// for safety...
		for (Edge e : frontList) {
			if (e.leftFrontNeighbor == null) {
				Msg.warning("e.leftFrontNeighbor is null.");
			} else if (e.leftFrontNeighbor == e) {
				Msg.warning("e.leftFrontNeighbor points to e.");
			}
			if (e.rightFrontNeighbor == null) {
				Msg.warning("e.rightFrontNeighbor is null.");
			} else if (e.rightFrontNeighbor == e) {
				Msg.warning("e.rightFrontNeighbor points to e.");
			}
		}

		return frontList;
	}

	private void classifyStateOfAllFronts(List<Edge> frontList) {
		for (Edge e : frontList) {
			e.classifyStateOfFrontEdge();
		}
	}

	/** Performs seaming operation as described in Owen's paper */
	private Quad doSeam(Edge e1, Edge e2, Vertex nK) {
		Msg.warning("Entering doSeam(..)...");
		Quad e1Quad = e1.getQuadElement();
		Vertex safePos = null;

		// We must avoid that nKp1 becomes a boundaryVertex because nKp1 is the Vertex that
		// is to
		// be relocated.
		Edge tmp;
		if (e1.otherVertex(nK).boundaryVertex()) {
			tmp = e1;
			e1 = e2;
			e2 = tmp;
		}

		// Clear everything inside triangle defined by e0, e1, and e2:
		Vertex nKm1 = e2.otherVertex(nK);
		Vertex nKp1 = e1.otherVertex(nK);
		Edge e0 = recoverEdge(nKm1, nKp1);

		if (e0 == null || e0.element1 instanceof Quad || e0.element2 instanceof Quad) {
			return null;
		}
		Triangle ta = (Triangle) e0.element1, tb = (Triangle) e0.element2;
		Edge et1, et2;
		Vertex nta = ta.oppositeOfEdge(e0), ntb = tb.oppositeOfEdge(e0), nT;
		if (nta.inHalfplane(e0, nK) == 1) {
			nT = ntb;
			et1 = tb.neighborEdge(e0.leftVertex, e0);
			et2 = tb.neighborEdge(e0.rightVertex, e0);
		} else {
			nT = nta;
			et1 = ta.neighborEdge(e0.leftVertex, e0);
			et2 = ta.neighborEdge(e0.rightVertex, e0);
		}

		Edge l, r, t;
		if (e1.leftVertex == nK) {
			l = e2;
			if (et1.hasVertex(e1.rightVertex)) {
				r = et1;
				t = et2;
			} else {
				r = et2;
				t = et1;
			}
		} else {
			r = e2;
			if (et1.hasVertex(e1.leftVertex)) {
				l = et1;
				t = et2;
			} else {
				l = et2;
				t = et1;
			}
		}

		Quad q = new Quad(e1, l, r, t);

		if (!nKm1.boundaryVertex()) {
			safePos = safeNewPosWhenCollapsingQuad(q, nKp1, nKm1); // thus: merged Vertex can be moved
		}

		if (safePos == null) {
			if (q.anyInvertedElementsWhenCollapsed(nKm1, nKp1, nKm1, nKp1.adjElements(), nKm1.adjElements())) {
				return null;
			} else {
				safePos = nKm1.copyXY();
			}
		}

		q.firstVertex = null; // indicate that this is not really quad
		clearQuad(q, e1.getTriangleElement());
		q.disconnectEdges();

		// Merge nKp1 and nKm1 at a location midway between their initial positions:
		q.closeQuad(e1, e2); // (Copy the edges belonging to nKm1 into nKp1.)

		nKp1.moveTo(safePos);

		/*
		 * if (!nKp1.boundaryVertex()) { if (safePos!= null) nKp1.moveTo(safePos); //
		 * e0.midPoint() }
		 */
		// Update: Remove the two edges (that have collapsed) of the quad from edgeList
		if (!l.hasVertex(nKp1)) {
			edgeList.remove(edgeList.indexOf(l));
		}
		if (!r.hasVertex(nKp1)) {
			edgeList.remove(edgeList.indexOf(r));
		}
		if (!t.hasVertex(nKp1)) {
			edgeList.remove(edgeList.indexOf(t));
		}

		vertexList.remove(vertexList.indexOf(nKm1));

		// A little trick so that we will be able to use localUpdateFront(..) later:
		e2.removeFromFront(frontList);
		e2.removeFromStateList();
		if (e2.level == level) {
			nrOfFronts--;
		}

		if (t.frontEdge) {
			t.removeFromFront(frontList);
			t.removeFromStateList();
			if (t.level == level) {
				nrOfFronts--;
			}
		}

		q.replaceEdge(e2, e1);
		q.replaceEdge(t, e1);

		// Do a local smooth (perhaps the other vertices involved could need some, too?)
		Vertex old;
		Vertex nKNew, nKp1New;

		if (!nK.boundaryVertex()) {
			nKNew = nK.laplacianSmooth(); // modifiedLWLaplacianSmooth();
		} else {
			nKNew = nK;
		}

		if (!nKp1.boundaryVertex()) {
			nKp1New = nKp1.laplacianSmooth(); // modifiedLWLaplacianSmooth();
		} else {
			nKp1New = nKp1;
		}

		if (!nK.equals(nKNew)) {
			old = nK.copyXY();
			nK.moveTo(nKNew);
			inversionCheckAndRepair(nK, old);
			nK.update();
		}

		if (!nKp1.equals(nKp1New)) {
			old = nKp1.copyXY();
			nKp1.moveTo(nKp1New);
			inversionCheckAndRepair(nKp1, old);
			nKp1.update();
		}

		Msg.warning("Leaving doSeam(..)...");
		return q;
	}

	/** Performs the transition seam operation as described in Owen's paper. */
	private Quad doTransitionSeam(Edge e1, Edge e2, Vertex nK) {
		Edge longer, shorter;
		if (e1.len > e2.len) {
			longer = e1;
			shorter = e2;
		} else {
			longer = e2;
			shorter = e1;
		}

		// Get all the edges, vertices, triangles, etc. and create 4 new edges
		Vertex mid = longer.midPoint();
		mid.color = java.awt.Color.blue; // Blue color indicates it was created by split
		Vertex nKm1 = longer.otherVertex(nK), nKp1 = shorter.otherVertex(nK);
		Edge frontBeyondKm1 = longer.frontNeighborAt(nKm1);
		Edge frontBeyondKp1 = shorter.frontNeighborAt(nKp1);
		Quad q = longer.getQuadElement();
		Triangle tL = longer.getTriangleElement();
		Edge eTLKm1 = tL.neighborEdge(nKm1, longer), eTLK = tL.neighborEdge(nK, longer);

		Edge eF, eFL, eMidKm1, eKMid, eMidTT;
		Vertex nF;

		if (q.largestAngle() > DEG_179) {// Then we can't use the approach in the Owen et al paper!
			eF = q.neighborEdge(nKm1, longer);
			eMidKm1 = new Edge(mid, nKm1);
			eKMid = new Edge(nK, mid);
			eMidTT = new Edge(mid, tL.oppositeOfEdge(longer));
			nF = eF.otherVertex(nKm1);
			eFL = new Edge(mid, nF);

			eFL.connectVertices();
			eMidKm1.connectVertices();
			eKMid.connectVertices();
			eMidTT.connectVertices();

			// Update edgeList
			edgeList.add(eFL);
			edgeList.add(eKMid);
			edgeList.add(eMidKm1);
			edgeList.add(eMidTT);

			Triangle t1New = new Triangle(eFL, eF, eMidKm1);
			Triangle t2New = new Triangle(eKMid, eTLK, eMidTT);
			Triangle t3New = new Triangle(eMidKm1, eMidTT, eTLKm1);

			Edge b = q.oppositeEdge(longer), l = q.neighborEdge(b.leftVertex, b), r = q.neighborEdge(b.rightVertex, b);
			if (l == eF) {
				l = eFL;
			} else {
				r = eFL;
			}
			Quad q1New = new Quad(b, l, r, eKMid);

			q.disconnectEdges();
			tL.disconnectEdges();

			q1New.connectEdges();

			t1New.connectEdges();
			t2New.connectEdges();
			t3New.connectEdges();

			// Update triangleList and elementList
			elementList.remove(elementList.indexOf(q));
			triangleList.remove(triangleList.indexOf(tL));

			triangleList.add(t1New);
			triangleList.add(t2New);
			triangleList.add(t3New);
			// elementList.add(q1New); // Will be added later...

			// Update edgeList
			edgeList.remove(edgeList.indexOf(longer));

			// Update vertexList
			longer.disconnectVertices();
			vertexList.add(mid);

			// Update front
			longer.removeFromStateList();
			longer.removeFromFront(frontList);

			eKMid.promoteToFront(longer.level, frontList);
			eFL.promoteToFront(longer.level, frontList);

			if (frontBeyondKm1 != eF) {
				eF.promoteToFront(longer.level, frontList);

				frontBeyondKm1.setFrontNeighbor(eF);

				eF.setFrontNeighbor(frontBeyondKm1);
				eF.setFrontNeighbor(eFL);

				eFL.setFrontNeighbor(eF);
				eFL.setFrontNeighbor(eKMid);

				eKMid.setFrontNeighbor(eFL);
				eKMid.setFrontNeighbor(shorter);

				shorter.setFrontNeighbor(eKMid);

				eF.classifyStateOfFrontEdge();

				nrOfFronts += 2;
			} else {
				Edge beyond = eF.frontNeighborAt(nF);
				eF.removeFromStateList();
				eF.removeFromFront(frontList);

				beyond.setFrontNeighbor(eFL);

				eFL.setFrontNeighbor(beyond);
				eFL.setFrontNeighbor(eKMid);

				eKMid.setFrontNeighbor(eFL);
				eKMid.setFrontNeighbor(shorter);

				shorter.setFrontNeighbor(eKMid);

				nrOfFronts += 0;
			}

			eKMid.classifyStateOfFrontEdge();
			eFL.classifyStateOfFrontEdge();

			localSmooth(q1New, frontList);
			nrOfFronts -= localUpdateFronts(q1New, level, frontList);

			return q1New;
		}

		eF = q.neighborEdge(nK, longer);
		nF = eF.otherVertex(nK);

		eFL = new Edge(mid, nF);
		eMidKm1 = new Edge(mid, nKm1);
		eKMid = new Edge(nK, mid);
		eMidTT = new Edge(mid, tL.oppositeOfEdge(longer));

		eFL.connectVertices();
		eMidKm1.connectVertices();
		eKMid.connectVertices();
		eMidTT.connectVertices();

		// Update edgeList
		edgeList.add(eFL);
		edgeList.add(eMidKm1);
		edgeList.add(eMidTT);
		// edgeList.add(eKMid);

		Edge b = q.oppositeEdge(longer), l, r;
		if (eF.hasVertex(b.leftVertex)) {
			l = eFL;
			r = q.neighborEdge(nKm1, longer);
		} else {
			l = q.neighborEdge(nKm1, longer);
			r = eFL;
		}

		q.disconnectEdges();
		tL.disconnectEdges();

		// Create 1 quad and 3 triangles
		Quad q1New = new Quad(b, l, r, eMidKm1);
		q1New.connectEdges();

		Triangle t1New = new Triangle(eF, eFL, eKMid), t2New = new Triangle(eMidKm1, eMidTT, eTLKm1),
				t3New = new Triangle(eKMid, eMidTT, eTLK);

		t1New.connectEdges();
		t2New.connectEdges();
		t3New.connectEdges();

		// Update triangleList and edgeList
		// triangleList.add(t1New);
		triangleList.add(t2New);
		triangleList.add(t3New);

		// Create second quad
		Edge top = recoverEdge(mid, nKp1);

		if (shorter.hasVertex(eF.leftVertex)) {
			l = shorter;
			r = eFL;
		} else {
			l = eFL;
			r = shorter;
		}
		Quad q2New = new Quad(eF, l, r, top);
		clearQuad(q2New, t1New);

		// Update, update, update....
		longer.disconnectVertices();
		edgeList.remove(edgeList.indexOf(longer));

		// if (longer.level== level)
		// nrOfFronts+= 2; // Longer and shorter is split into two front edges

		longer.removeFromStateList();
		longer.removeFromFront(frontList);

		shorter.removeFromStateList();
		shorter.removeFromFront(frontList);

		eMidKm1.promoteToFront(longer.level, frontList);
		top.promoteToFront(longer.level + 1, frontList);

		frontBeyondKm1.setFrontNeighbor(eMidKm1);

		eMidKm1.setFrontNeighbor(frontBeyondKm1);
		eMidKm1.setFrontNeighbor(top);

		top.setFrontNeighbor(eMidKm1);
		top.setFrontNeighbor(frontBeyondKp1);

		frontBeyondKp1.setFrontNeighbor(top);

		eMidKm1.classifyStateOfFrontEdge();
		top.classifyStateOfFrontEdge();

		// Update vertexList
		vertexList.add(mid);

		// Update triangleList and elementList
		elementList.remove(elementList.indexOf(q));
		elementList.add(q1New);

		triangleList.remove(triangleList.indexOf(tL));

		localSmooth(q1New, frontList);
		nrOfFronts -= localUpdateFronts(q1New, level, frontList);

		return q2New; // This quad will be added to elementList shortly.
	}

	/** Performs the transition split operation as described in Owen's paper. */
	private Quad doTransitionSplit(Edge e1, Edge e2, Vertex nK) {
		Edge longer, shorter;
		if (e1.len > e2.len) {
			longer = e1;
			shorter = e2;
		} else {
			longer = e2;
			shorter = e1;
		}

		// Get all edges, vertices, etc. and also create some new edges.
		Quad q1 = longer.getQuadElement();
		Triangle t1 = longer.getTriangleElement();

		Vertex nKm1 = longer.otherVertex(nK), nKp1 = shorter.otherVertex(nK), opposite = q1.oppositeVertex(nK);

		Edge frontBeyondKp1 = shorter.frontNeighborAt(nKp1);
		Edge frontBeyondKm1 = longer.frontNeighborAt(nKm1);

		Edge q1Top = q1.oppositeEdge(longer), q1nK = q1.neighborEdge(nK, longer), eKm1Opp = q1.oppositeEdge(q1nK);
		Edge eT1K = t1.neighborEdge(nK, longer), eT1Km1 = t1.neighborEdge(nKm1, longer);
		Vertex c = q1.centroid(), mid = longer.midPoint();
		Edge eF = new Edge(nK, c), eFL = new Edge(c, mid), eCOpp = new Edge(c, opposite), eMidTT = new Edge(mid, t1.oppositeOfEdge(longer)),
				eKMid = new Edge(nK, mid), eMidKm1 = new Edge(mid, nKm1);

		longer.disconnectVertices();

		eF.connectVertices();
		eFL.connectVertices();
		eCOpp.connectVertices();
		eMidTT.connectVertices();
		eKMid.connectVertices();
		eMidKm1.connectVertices();

		q1.disconnectEdges();
		t1.disconnectEdges();

		// Longer is split into three front edges
		if (longer.level == level) {
			nrOfFronts += 2;
		}

		// Remove longer from the front
		longer.removeFromStateList();
		longer.removeFromFront(frontList);

		Edge l, r;
		if (q1nK.hasVertex(q1Top.leftVertex)) {
			l = q1nK;
			r = eCOpp;
		} else {
			l = eCOpp;
			r = q1nK;
		}
		Quad q11New = new Quad(q1Top, l, r, eF);
		q11New.connectEdges();

		if (eFL.hasVertex(eCOpp.leftVertex)) {
			l = eFL;
			r = eKm1Opp;
		} else {
			l = eKm1Opp;
			r = eFL;
		}
		Quad q12New = new Quad(eCOpp, l, r, eMidKm1);
		q12New.connectEdges();

		Triangle t1New = new Triangle(eF, eKMid, eFL);
		Triangle t2New = new Triangle(eKMid, eMidTT, eT1K);
		Triangle t3New = new Triangle(eMidKm1, eMidTT, eT1Km1);
		t1New.connectEdges();
		t2New.connectEdges();
		t3New.connectEdges();

		// Update stuff....
		vertexList.add(c);
		vertexList.add(mid);
		c.color = java.awt.Color.blue;
		mid.color = java.awt.Color.blue;

		triangleList.add(t1New);
		triangleList.add(t2New);
		triangleList.add(t3New);

		triangleList.remove(t1);

		int i = elementList.indexOf(q1);
		if (i != -1) {
			elementList.remove(i);
		}
		elementList.add(q11New);

		// Remember to remove longer from stateList and frontList...
		// ...set new front edges ... ...and set new frontNeighbors....

		eMidKm1.promoteToFront(q1.edgeList[base].level + 1, frontList);
		eFL.promoteToFront(q1.edgeList[base].level + 1, frontList);
		eF.promoteToFront(q1.edgeList[base].level + 1, frontList);

		frontBeyondKm1.setFrontNeighbor(eMidKm1);

		eMidKm1.setFrontNeighbor(frontBeyondKm1);
		eMidKm1.setFrontNeighbor(eFL);

		eFL.setFrontNeighbor(eMidKm1);
		eFL.setFrontNeighbor(eF);

		eF.setFrontNeighbor(eFL);
		eF.setFrontNeighbor(shorter);

		shorter.setFrontNeighbor(eF);

		eF.leftSide = true; // Edge eF is classified as a state 1-1 edge, regardless..
		eF.rightSide = true;
		Edge.stateList[2].add(eF);

		eMidKm1.classifyStateOfFrontEdge();
		eFL.classifyStateOfFrontEdge();

		edgeList.add(eF);
		edgeList.add(eFL);
		edgeList.add(eKMid);
		edgeList.add(eCOpp);
		edgeList.add(eMidTT);
		edgeList.add(eMidKm1);

		edgeList.remove(edgeList.indexOf(longer));

		localSmooth(q11New, frontList);
		localUpdateFronts(q11New, level, frontList); // Param 4 is dummy

		if (q11New.inverted()) {
			Msg.warning("...q11New was inverted initially!!");
		} else {
		}

		if (q12New.inverted()) {
			Msg.warning("...q12New was inverted initially!!");
		} else {
		}

		return q12New; // ...so that q12New will be smoothed and updated as well.
	}

	private boolean canSeam(Vertex n, Edge e1, Edge e2) {
		Vertex n1 = e1.otherVertex(n);
		Vertex n2 = e1.otherVertex(n);

		if (!n1.boundaryVertex() || !n2.boundaryVertex()) {
			return true;
		} else {
			return false;
		}
	}

	/**
	 * Check wether a seam operation is needed.
	 * 
	 * @param e1 an edge
	 * @param e2 an edge that might need to be merged with e1
	 * @param n  the common Vertex of e1 and e2
	 * @param nQ number of Quads adjacent to Vertex n
	 */
	private boolean needsSeam(Edge e1, Edge e2, Vertex n, int nQ) {
		Element elem;
		double ang;

		elem = e1.getTriangleElement();
		if (elem == null) {
			return false;
		}
		ang = e1.sumAngle(elem, n, e2);

		if (nQ >= 5) {
			if (ang < EPSILON1) {
				return true;
			} else {
				return false;
			}
		} else if (ang < EPSILON2) {
			return true;
		} else {
			return false;
		}
	}

	/**
	 * @return not null if a special case was encountered. If so, this means that
	 *         that the base edge has been destroyed, and that a new one has to be
	 *         chosen.
	 */
	private Quad handleSpecialCases(Edge e) {
		Quad q = null;
		Triangle eTri = e.getTriangleElement();
		if (!e.boundaryEdge() && !e.leftFrontNeighbor.boundaryEdge() && !e.rightFrontNeighbor.boundaryEdge()) {

			// Check for special cases:
			int nQLeft = e.leftVertex.nrOfAdjQuads();
			if (needsSeam(e, e.leftFrontNeighbor, e.leftVertex, nQLeft) && canSeam(e.leftVertex, e, e.leftFrontNeighbor)) {

				if (e.isLargeTransition(e.leftFrontNeighbor)) {
					q = doTransitionSeam(e, e.leftFrontNeighbor, e.leftVertex);
				} else {
					q = doSeam(e, e.leftFrontNeighbor, e.leftVertex);
				}
			} else if (e.getState() != 2 && e.isLargeTransition(e.leftFrontNeighbor)
					&& e.sumAngle(eTri, e.leftVertex, e.leftFrontNeighbor) < Math.PI) {
				q = doTransitionSplit(e, e.leftFrontNeighbor, e.leftVertex);
			}

			int nQRight = e.rightVertex.nrOfAdjQuads();
			if (needsSeam(e, e.rightFrontNeighbor, e.rightVertex, nQRight) && canSeam(e.rightVertex, e, e.rightFrontNeighbor)) {

				if (e.isLargeTransition(e.rightFrontNeighbor)) {
					q = doTransitionSeam(e, e.rightFrontNeighbor, e.rightVertex);
				} else {
					q = doSeam(e, e.rightFrontNeighbor, e.rightVertex);
				}
			} else if (e.getState() != 2 && e.isLargeTransition(e.rightFrontNeighbor)
					&& e.sumAngle(eTri, e.rightVertex, e.rightFrontNeighbor) < Math.PI) {
				q = doTransitionSplit(e, e.rightFrontNeighbor, e.rightVertex);
			}
		}

		return q;
	}

	// jobber her...
	/**
	 * If not already set, this method sets the sideEdges, based on state. Called
	 * from makeQuad(Edge)
	 */
	private Edge[] makeSideEdges(Edge e) {
		Edge[] sideEdges = new Edge[2];
		List<Edge> altLSE;
		List<Edge> altRSE;

		int lState = 0, rState = 0;
		if (e.leftSide) {
			lState = 1;
		}
		if (e.rightSide) {
			rState = 1;
		}

		if (lState == 1 && rState == 1) { // Then side edges are already set.
			sideEdges[0] = e.leftFrontNeighbor;
			sideEdges[1] = e.rightFrontNeighbor;

			/*
			 * if (altLSE.size()==1 && altRSE.size()==1) { // Make sure that we use legal
			 * edges e.leftSide= (Edge)altLSE.get(0); e.rightSide= (Edge)altRSE.get(0); }
			 */
		} else if (lState == 1 && rState == 0) {
			sideEdges[0] = e.leftFrontNeighbor;
			sideEdges[0].swappable = false;
			altRSE = getPotSideEdges(e, e.rightVertex);
			if (altRSE.size() > 1) {
				sideEdges[1] = defineSideEdge(e, e.rightVertex, sideEdges[0], null, altRSE);
			} else {
				sideEdges[1] = (Edge) altRSE.get(0);
			}
			sideEdges[0].swappable = true;
		} else if (lState == 0 && rState == 1) {
			sideEdges[1] = e.rightFrontNeighbor;
			sideEdges[1].swappable = false;
			altLSE = getPotSideEdges(e, e.leftVertex);
			if (altLSE.size() > 1) {
				sideEdges[0] = defineSideEdge(e, e.leftVertex, null, sideEdges[1], altLSE);
			} else {
				sideEdges[0] = (Edge) altLSE.get(0);
			}
			sideEdges[1].swappable = true;
		} else if (lState == 0 && rState == 0) {
			altLSE = getPotSideEdges(e, e.leftVertex);
			if (altLSE.size() > 1) {
				sideEdges[0] = defineSideEdge(e, e.leftVertex, null, null, altLSE);
			} else {
				sideEdges[0] = (Edge) altLSE.get(0);
			}

			sideEdges[0].swappable = false;
			altRSE = getPotSideEdges(e, e.rightVertex);
			if (altRSE.size() > 1) {
				sideEdges[1] = defineSideEdge(e, e.rightVertex, sideEdges[0], null, altRSE);
			} else {
				sideEdges[1] = (Edge) altRSE.get(0);
			}
			sideEdges[0].swappable = true;
		}

		// Now, if we're about to close a front, then we have to check whether we get
		// loops with an even or odd number of edges. If this number is odd, then the
		// side edge must be split, so that a subsequent side edge selection may create
		// loops with an even number of edges.
		Edge lSide = sideEdges[0], rSide = sideEdges[1];
		byte lrc;
		boolean lLoop = false, rLoop = false;

		if (((lSide.otherVertex(e.leftVertex).frontVertex() && !lSide.isFrontEdge())
				|| (rSide.otherVertex(e.rightVertex).frontVertex() && !rSide.isFrontEdge()))) {

			/*
			 * if (evenInitNrOfFronts && ( (lSide.otherVertex(e.leftVertex).frontVertex() &&
			 * !lSide.frontEdge) || (rSide.otherVertex(e.rightVertex).frontVertex() &&
			 * !rSide.frontEdge)) ) {
			 */

			

			// nM (the "otherVertex") lies on an opposing front, and if the side Edge
			// isn't a front edge, it may only be used if the number of edges on
			// each resulting front loop is even.
			if (lSide.otherVertex(e.leftVertex).frontVertex() && !lSide.isFrontEdge()) {
			} else {
			}

			if (/* lSide.getQuadElement()== null && */ !lSide.boundaryEdge() && lSide.otherVertex(e.leftVertex).frontVertex()
					&& !lSide.isFrontEdge()) {
				lLoop = true;
			}

			if (/* rSide.getQuadElement()== null && */ !rSide.boundaryEdge() && rSide.otherVertex(e.rightVertex).frontVertex()
					&& !rSide.isFrontEdge()) {
				rLoop = true;
			}

			lrc = oddNOFEdgesInLoopsWithFEdge(e, lSide, rSide, lLoop, rLoop);
			if ((lrc & 1) == 1) {
				sideEdges[0] = lSide.splitTrianglesAtMyMidPoint(triangleList, edgeList, vertexList, e);
			}
			if (((lrc & 2) == 2) && ((lrc & 4) == 0)) {
				sideEdges[1] = rSide.splitTrianglesAtMyMidPoint(triangleList, edgeList, vertexList, e);
			}

		}

		/*
		 * // There are more than one alt side edge for both left and right if
		 * (altLSE.size()> 1 && altRSE.size()> 1) { e.leftSide= defineSideEdge(e,
		 * e.leftVertex, PIdiv6, altLSE); e.leftSide.swappable= false; e.rightSide=
		 * defineSideEdge(e, e.rightVertex, PIdiv6, altRSE); e.leftSide.swappable= true; }
		 * // Only right side has more than one alt side edge else if (altLSE.size()== 1
		 * && altRSE.size()> 1) { e.leftSide= (Edge)altLSE.get(0); e.leftSide.swappable=
		 * false; e.rightSide= defineSideEdge(e, e.rightVertex, PIdiv6, altRSE);
		 * e.leftSide.swappable= true; } // Only left side has more than one alt side
		 * edge else if (altLSE.size()>1 && altRSE.size()== 1) { e.rightSide=
		 * (Edge)altRSE.get(0); e.rightSide.swappable= false; e.leftSide=
		 * defineSideEdge(e, e.leftVertex, PIdiv6, altLSE); e.rightSide.swappable= true; }
		 * // Neither left or right has more than one alt side edge else { e.rightSide=
		 * (Edge)altRSE.get(0); e.leftSide= (Edge)altLSE.get(0); }
		 */

		return sideEdges;
	}

	/** @return a list of potential side edges for this base edge at Vertex n */
	private List<Edge> getPotSideEdges(Edge baseEdge, Vertex n) {
		Edge cur;
		List<Edge> list = new ArrayList<>();

		for (int i = 0; i < n.edgeList.size(); i++) {
			cur = (Edge) n.edgeList.get(i);
			if (cur != baseEdge && cur.bordersToTriangle()) {
				list.add(cur);
			}
		}
		printEdgeList(list);
		if (list.size() > 0) {
			return list;
		} else {
			Msg.error("Cannot find a side edge for this base edge.");
			return null;
		}
	}

	/**
	 * @param eF1  the base edge of the quad to be created
	 * @param nK   the Vertex at eF1 that the side edge will be connected to
	 * @param list list of candidate edges from which we might select a side edge
	 * @return A side edge: a reused edge OR one created in a swap/split operation
	 */
	private Edge defineSideEdge(Edge eF1, Vertex nK, Edge leftSide, Edge rightSide, List<Edge> list) {
		Edge selected = null, closest, eF2;
		Vertex noVertex = null;
		double curAng, selAng, closestAng;
		Element curElement, curElement2;

		
		// Set eF2 and noVertex (a Vertex that must not connect with the edge
		// which we are trying to find), or if the state bits permits, return the
		// already found side edge.
		if (nK == eF1.leftVertex) {
			eF2 = eF1.leftFrontNeighbor;
			if (leftSide != null) { // the left bit set
				return leftSide;
			}
			if (rightSide != null) {
				noVertex = rightSide.otherVertex(eF1.rightVertex);
			}
		} else {
			eF2 = eF1.rightFrontNeighbor;
			if (rightSide != null) { // the right bit set
				return rightSide;
			}
			if (leftSide != null) {
				noVertex = leftSide.otherVertex(eF1.leftVertex);
			}
		}

		if (noVertex != null) {
		}

		// Select the neighboring elements ahead of eF1 and eF2
		curElement = eF1.getTriangleElement();
		if (curElement == null) {
			Msg.error("...oisann!");
		}
		curElement2 = eF2.getTriangleElement();
		if (curElement2 == null) {
			if (eF2.element1 != null) {
				curElement2 = eF2.element1;
			} else {
				Msg.error("defineSideEdge: Cannot find element ahead of " + eF2.descr());
				return null;
			}
		}

		double bisected = eF1.sumAngle(curElement, nK, eF2) / 2.0;
		

		// *TRY* to select an initial edge that is not EF1 or EF2. If that's not
		// possible, then select EF2. Compute angle between the selected edge and
		// vector Vk (see figure...uhm). However, the selected edge should not
		// contain the noVertex.

		if (!eF2.hasVertex(noVertex)) {
			selected = eF2;
		} else {
			for (Edge current : list) {
				if (current != selected && current != eF2 && !current.hasVertex(noVertex)) {
					selected = current;
					break;
				}
			}
		}

		if (selected == null) { // If no edges are found that don't contain the noVertex
			selected = eF2; // then we have to select eF2
			selAng = Math.abs(bisected - eF1.sumAngle(curElement, nK, selected));
			closest = selected;
			closestAng = selAng;
		} else {
			selAng = Math.abs(bisected - eF1.sumAngle(curElement, nK, selected));
			closest = selected;
			closestAng = selAng;

			for (Edge current : list) {
				if (current != selected && current != eF2) {
					curAng = Math.abs(bisected - eF1.sumAngle(curElement, nK, current));

					if (curAng < closestAng) {
						closestAng = curAng;
						closest = current;
					}

					if (curAng < selAng && !current.hasVertex(noVertex)) {
						selAng = curAng;
						selected = current;
					}

				}
			}
		}

		if (selAng < EPSILON) {
			return selected;
		}

		// According to the algorithm, when nM (the opposite Vertex of nK on the selected
		// Edge) lies on an opposing front, I should allow for a larger angle
		// (EPSILONLARGER). If that fails, well, the program fails. So EPSILONLARGER
		// should really be set to a laaaaarge value, or, I could accept all values.
		// Yeah, that sounds like a good idea:
		if (selected.otherVertex(nK).frontVertex()) {
			// if (selAng < EPSILONLARGER) {
			return selected;
			/*
			 * } else
			 * Msg.error("Cannot close front because value EPSILONLARGER is too small."+
			 * "Needed value: larger than "+Math.toDegrees(selAng)+
			 * ", EPSILONLARGER== "Math.toDegrees(EPSILONLARGER));
			 */
		}

		

		// First find the triangle that vK passes through.
		// This is the plan:
		// The edge named "closest" is the edge that lies closest to vK. Closest belongs
		// to (one or) two elements. We check these two elements whether the bisected
		// angle passes through. (That is, bisected is larger than the first angle,
		// but smaller than the other...) When we have correctly identified
		// the element, e0 is that *triangle's* third edge.

		Element elem1 = closest.element1, elem2 = closest.element2;
		Triangle bisectTriangle;
		Edge otherEdge = elem1.neighborEdge(nK, closest), e0;

		double a1 = eF1.sumAngle(curElement, nK, closest);
		double a2 = eF1.sumAngle(curElement, nK, otherEdge);

		double ang1 = Math.min(a1, a2);
		double ang2 = Math.max(a1, a2);

		if ((bisected >= ang1 && bisected <= ang2) || elem2 == null) { // Does vK go through elem1? or is elem1 the only element connected
																		// to selected?
			if (elem1 instanceof Quad) {
				return selected; // have to give up
			}
			// selInd= elem1.indexOf(selected);
			// otherInd= elem1.indexOf(otherEdge);
			bisectTriangle = (Triangle) elem1;
		} else if (elem2 != null) { // Nope, it seems vK goes through elem2, not elem1.
			if (elem2 instanceof Quad) {
				return selected; // have to give up
			}
			// selInd= elem2.indexOf(selected);
			otherEdge = elem2.neighborEdge(nK, closest);
			// otherInd= elem2.indexOf(otherEdge);

			bisectTriangle = (Triangle) elem2;
		} else {
			Msg.error("defineSideEdge(..): Cannot find which element vK runs through.");
			return null;
		}

		// Now get e0:
		e0 = bisectTriangle.otherEdge(closest, otherEdge);

		// e0 cannot be a front edge, that is already ruled out in the reuse part.
		// If e0 belongs to a quad or lies at the boundary,
		// I cannot swap, so I just return:
		Element neighborTriangle = bisectTriangle.neighbor(e0);
		if (e0.frontEdge || neighborTriangle == null || neighborTriangle instanceof Quad) {
			return selected;
		}
		// Then, use e0 to get to nM.
		Vertex nM = e0.oppositeVertex(nK);
		Edge eK = new Edge(nK, nM); // new edge Ek

		// The angle between eF1 and eK has to be less than 180 degrees. The same goes
		// for the angle between eF2 and selected. Therefore, we can use
		// computePosAngle here.
		double beta = Math.abs(bisected - Math.min(eF1.computePosAngle(eK, nK), eF2.computePosAngle(eK, nK)));

		if (beta < EPSILON && eK.len < (eF1.len + eF2.len) * sqrt3div2) {
			
			// ok, swap edges: remove e0 and introduce eK.... update stuff...

			// Remove old triangles from list...
			triangleList.remove(triangleList.indexOf(e0.element1));
			triangleList.remove(triangleList.indexOf(e0.element2));

			// ... swap ...
			e0.swapToAndSetElementsFor(eK);

			// ... and replace with new ones:
			triangleList.add((Triangle) eK.element1);
			triangleList.add((Triangle) eK.element2);

			// Update "global" edge list
			edgeList.remove(edgeList.indexOf(e0));
			edgeList.add(eK);
			return eK;
		} else {
			// Then the only way this might work is splitting bisectTriangle
			// and it's neighbor at edge e0.
			MyVector vEF1 = new MyVector(nK, eF1.otherVertex(nK));
			MyVector vEF2 = new MyVector(nK, eF2.otherVertex(nK));

			// Take special care to construct Ray rK correctly:
			Ray rK;
			if (vEF2.isCWto(vEF1)) {
				if (bisected < PIdiv2) {
					rK = new Ray(nK, eF2, bisected);
				} else {
					rK = new Ray(nK, eF1, bisected);
				}
			} else if (bisected < PIdiv2) {
				rK = new Ray(nK, eF1, bisected);
			} else {
				rK = new Ray(nK, eF2, bisected);
			}

			// MyVector vK= new MyVector(eF1, nK,bisected, Double.POSITIVE_INFINITY);
			// (eF1, nK,bisected, Double.POSITIVE_INFINITY);

			MyVector v0 = e0.getVector();
			Vertex nN = rK.pointIntersectsAt(v0);
			if (nN == null) {
				Msg.error("defineSideEdge(..): Cannot split edge e0==" + e0.descr());
			}
			if (!vertexList.contains(nN)) {
				vertexList.add(nN);
			}

			nN.color = java.awt.Color.blue;

			// Pick some more necessary vertices and edges:
			Vertex nO = e0.commonVertex(closest);
			Vertex nP = e0.commonVertex(otherEdge);
			Edge eO = neighborTriangle.neighborEdge(nO, e0);
			Edge eP = neighborTriangle.neighborEdge(nP, e0);

			// Create 4 new edges:
			eK = new Edge(nK, nN);
			Edge eM = new Edge(nN, nM);
			Edge e1 = new Edge(nN, nO);
			Edge e2 = new Edge(nN, nP);

			// Update "local" edgeLists... at each affected Vertex:
			eK.connectVertices();
			eM.connectVertices();
			e1.connectVertices();
			e2.connectVertices();

			e0.disconnectVertices();

			// Update "global" edgeList:
			edgeList.add(eK);
			edgeList.add(eM);
			edgeList.add(e1);
			edgeList.add(e2);
			edgeList.remove(edgeList.indexOf(e0));

			// Update "global" triangleList, disconnect edges, etc:
			triangleList.remove(triangleList.indexOf(neighborTriangle));
			triangleList.remove(triangleList.indexOf(bisectTriangle));

			bisectTriangle.disconnectEdges();

			Triangle ta = new Triangle(eK, e1, closest);
			triangleList.add(ta);
			Triangle tb = new Triangle(eK, e2, otherEdge);
			triangleList.add(tb);

			neighborTriangle.disconnectEdges();

			Triangle tc = new Triangle(eM, e1, eO);
			triangleList.add(tc);
			Triangle td = new Triangle(eM, e2, eP);
			triangleList.add(td);

			ta.connectEdges();
			tb.connectEdges();
			tc.connectEdges();
			td.connectEdges();

			Msg.warning("Splitting is not yet thouroughly tested.");
			return eK;
		}
	}

	/**
	 * Create an edge from nC to nD, or if it already exists, return that edge.
	 * Remove all edges intersected by the new edge. This is accomplished through a
	 * swapping procedure. All local and "global" updating is taken care of.
	 * 
	 * @param nC the start Vertex
	 * @param nD the end Vertex
	 * @return null if the Edge could not be recovered. This might happen when the
	 *         line segment from nC to nD intersects a Quad or a boundary Edge.
	 */
	private Edge recoverEdge(Vertex nC, Vertex nD) {
		Edge S = new Edge(nD, nC);
		printEdgeList(nC.edgeList);

		if (nC.edgeList.contains(S)) {
			Edge edge = nC.edgeList.get(nC.edgeList.indexOf(S));
			return edge;
		}
		/* ---- First find the edges connecting vertices nC and nD: ---- */
		/* (Implementation of algorithm 2 in Owen) */
		List<Edge> intersectedEdges = new ArrayList<>();
		Edge eK, eKp1, eI = null, eJ = null, eN, eNp1;
		Element tK = null;
		Triangle tI = null, tIp1;
		MyVector vK, vKp1 = null, vI;
		MyVector vS = new MyVector(nC, nD);
		Vertex nI;

		List<MyVector> V = nC.ccwSortedVectorList();
		V.add(V.get(0)); // First add first edge to end of list to avoid crash in loops
		printVectors(V);

		// Aided by V, fill T with elements adjacent nC, in
		// ccw order (should work even for vertices with only two edges in their list):
		List<Element> T = new ArrayList<>();

		for (int k = 0; k < V.size() - 1; k++) {
			vK = (MyVector) V.get(k);
			vKp1 = (MyVector) V.get(k + 1);
			eK = vK.edge;
			eKp1 = vKp1.edge;

			if (eK == null || eKp1 == null) {
			}
			if (!T.contains(eK.element1) && (eK.element1 == eKp1.element1 || eK.element1 == eKp1.element2)) {
				T.add(eK.element1);
			}
			if (eK.element2 != null && !T.contains(eK.element2) && (eK.element2 == eKp1.element2 || eK.element2 == eKp1.element1)) {
				T.add(eK.element2);
			}
		}

		// Now, get the element attached to nC that contains a part of S, tI:
		for (int k = 0; k < T.size(); k++) {
			vK = (MyVector) V.get(k);
			vKp1 = (MyVector) V.get(k + 1);

			tK = (Element) T.get(k);

			// We could optimize by using isCWto(..) directly instead of dot(..)
			// And I can't get the dot(...) to work properly anyway, sooo...

			
			// Msg.debug("vS.dot(vK)=="+vS.dot(vK)+" and vS.dot(vKp1)=="+vS.dot(vKp1));

			if (!vS.isCWto(vK) && vS.isCWto(vKp1)) {
				// tI= tK;
				// eI= vKp1.edge; // just something I tried.. no good
				break; // got it, escape from loop
			}
		}
		if (tK == null) {
			Msg.error("Oida.. valgt tK er null");
		}

		Element elemI, elemIp1;
		elemI = tK;
		if (elemI instanceof Quad) {
			Msg.warning("Leaving recoverEdge(..): intersecting quad, returning null.");
			return null;
		} else { // elemI is a Triangle...
			tI = (Triangle) elemI;
			eI = tI.oppositeOfVertex(nC);
		}
		if (!eI.frontEdge) {
			intersectedEdges.add(eI);
		} else {
			Msg.warning("Leaving recoverEdge: eI=" + eI.descr() + " is part of the front.");
			return null;
		}

		// Quad qI;
		while (true) {
			elemIp1 = elemI.neighbor(eI);
			// tIp1= (Triangle)elem;
			if (elemIp1.hasVertex(nD)) {
				break;
			}
			elemI = elemIp1;
			if (elemI instanceof Triangle) {
				tI = (Triangle) elemI;
				nI = tI.oppositeOfEdge(eI);
				vI = new MyVector(nC, nI);
				eN = tI.nextCCWEdge(eI);
				eNp1 = tI.nextCWEdge(eI);

				
				// if (vS.dot(vI)<0) // Not convinced that dot(..) works properly
				if (vS.isCWto(vI)) {
					eI = eN;
				} else {
					eI = eNp1;
				}

				if (!eI.frontEdge) {
					intersectedEdges.add(eI);
				} else {
					Msg.warning("Leaving recoverEdge: eI=" + eI.descr() + " is part of the front.");
					return null;
				}
			} else { // elemI is instanceof Quad
				Msg.warning("Leaving recoverEdge: intersecting quad.");
				return null;
			}
		}

		/* ---- ---- ---- ---- ---- ---- */

		/* ---- Secondly, do the actual recovering of the top edge ---- */
		if (intersectedEdges.size() == 0) {
			eJ = S;
			Msg.warning("recoverEdge: intersectedEdges.size()== 0");
		} else if (intersectedEdges.size() == 1) {
			// According to alg. 2 in Owen, intersectedEdges always contains at least
			// one edge.
			eI = (Edge) intersectedEdges.get(0);
			if (eI.equals(S)) {
				return eI;
			}
		}

		

		// When this loop is done, the edge should be recovered
		Element oldEIElement1, oldEIElement2;
		ArrayList<Element> removeList = new ArrayList<Element>();
		Triangle t;
		Quad q;
		int index;
		int n = intersectedEdges.size();
		// for (int i=0; i< n; i++) {
		Element old1, old2;
		Vertex na, nb, nc, nd;
		Edge qe1, qe2;

		while (intersectedEdges.size() > 0) {
			eI = (Edge) intersectedEdges.get(0);
			

			// We must avoid creating inverted or degenerate triangles.
			q = new Quad(eI);
			old1 = eI.element1;
			old2 = eI.element2;

			na = ((Triangle) old1).oppositeOfEdge(eI);
			nb = ((Triangle) old2).oppositeOfEdge(eI);
			nc = eI.leftVertex;
			nd = eI.rightVertex;

			double cross1 = cross(na, nc, nb, nc); // The cross product nanc x nbnc
			double cross2 = cross(na, nd, nb, nd); // The cross product nand x nbnd

			// For a stricly convex quad, the vertices must be located on different sides of
			// eI:
			if ((cross1 > 0 && cross2 < 0) || (cross1 < 0 && cross2 > 0) /* q.isStrictlyConvex() */) {
				// ... but this seems ok
				eI.swappable = true;
				eJ = eI.getSwappedEdge();

				eI.swapToAndSetElementsFor(eJ);

				if (eJ.element1.inverted()) {
					Msg.error("eJ.element1 is inverted: " + eJ.element1.descr());
				}

				if (eJ.element2.inverted()) {
					Msg.error("eJ.element2 is inverted: " + eJ.element2.descr());
				}

				// Add old triangles to removeList...
				if (!removeList.contains(old1)) {
					removeList.add(old1);
				}

				if (!removeList.contains(old2)) {
					removeList.add(old2);
				}

				// ... and replace with new ones:
				triangleList.add((Triangle) eJ.element1);
				triangleList.add((Triangle) eJ.element2);
				

				// Update "global" edge list
				edgeList.remove(edgeList.indexOf(eI));
				edgeList.add(eJ);

				intersectedEdges.remove(0);
				MyVector vEj = new MyVector(eJ.leftVertex, eJ.rightVertex);

				if (!eJ.hasVertex(nC) && !eJ.hasVertex(nD) && vEj.innerpointIntersects(vS)) {
					intersectedEdges.add(eJ);
				}

			} else if (intersectedEdges.size() > 1 && eI.swappable) {
				intersectedEdges.remove(0); // put eI last in the list
				intersectedEdges.add(eI);
				eI.swappable = false;
			} else {
				// We have to give up
				Msg.warning("Leaving recoverEdge: Cannot swap edge - non-convex quad.");

				for (Object intersectedEdge : intersectedEdges) {
					eI = (Edge) intersectedEdge;
					eI.swappable = true;
				}
				// Remove all triangles found in the removeList from triangleList:
				for (int i = 0; i < removeList.size(); i++) {
					t = (Triangle) removeList.get(i);
					index = triangleList.indexOf(t);
					if (index != -1) {
						triangleList.remove(index);
					}
				}
				return null;
			}
		}

		// Remove all triangles found in the removeList from triangleList:
		for (int i = 0; i < removeList.size(); i++) {
			t = (Triangle) removeList.get(i);
			index = triangleList.indexOf(t);
			if (index != -1) {
				triangleList.remove(index);
			}
		}

		/* ---- ---- ---- ---- ---- ---- */

		// eJ should be the recovered edge now, according to fig.6(d) in Owen.
		// From the swapToAndSetElementsFor(..) method, eJ has already got its
		// adjacent triangles. So everything should be juuuuust fine by now...
		return eJ;
	}

} // End of class QMorph
