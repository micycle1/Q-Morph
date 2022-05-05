package meshditor;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.List;

/**
 * This class holds information for vertices, and has methods for the management
 * of issues regarding vertices.
 */
public class Vertex extends Constants {

	/** Boolean indicating whether the Vertex has been moved by the OBS */
	public boolean movedByOBS = false; // Used by the smoother
	/** The coordinates */
	public double x, y;
	/** A valence pattern for this Vertex */
	public byte[] pattern;
	// byte state= 0; // For front Vertices only
	List<Edge> edgeList;
	java.awt.Color color = java.awt.Color.cyan;

	/** Create new Vertex with position (x,y). */
	public Vertex(double x, double y) {
		this.x = x;
		this.y = y;
		edgeList = new ArrayList<>();
	}

	@Override
	public boolean equals(Object elem) {
		if (!(elem instanceof Vertex)) {
			return false;
		}
		Vertex vertex = (Vertex) elem;
		if (x == vertex.x && y == vertex.y) {
			return true;
		} else {
			return false;
		}
	}

	@Override
	public int hashCode() {
		int result = 17;
		result = 37 * result + hashCode(x);
		result = 37 * result + hashCode(y);
		return result;
	}

	/**
	 * Computes a hash code for a double value, using the algorithm from Joshua
	 * Bloch's book <i>Effective Java"</i>
	 * 
	 * @param x the value to compute for
	 * @return a hashcode for x
	 */
	private static int hashCode(double x) {
		long f = Double.doubleToLongBits(x);
		return (int) (f ^ (f >>> 32));
	}

	/** @return a "real" copy of this Vertex with a shallow copy of its edgeList. */
	public Vertex copy() {
		Vertex n = new Vertex(x, y);
		n.edgeList = new ArrayList<>(edgeList);
		return n;
	}

	/** @return a new Vertex with the same position as this. */
	public Vertex copyXY() {
		return new Vertex(x, y);
	}

	/** Relocate this Vertex to the same position as n. */
	public void setXY(Vertex n) {
		x = n.x;
		y = n.y;
	}

	/** Relocate this Vertex to position (x,y) */
	public void setXY(double x, double y) {
		this.x = x;
		this.y = y;
	}

	void update() {
		updateLRinEdgeList();
		updateEdgeLengths();
		updateAngles();
	}

	public void updateLRinEdgeList() {
		for (Edge e : edgeList) {
			if ((e.leftVertex.x > e.rightVertex.x) || (e.leftVertex.x == e.rightVertex.x && e.leftVertex.y < e.rightVertex.y)) {
				Vertex vertex = e.leftVertex;
				e.leftVertex = e.rightVertex;
				e.rightVertex = vertex;

				if (e.frontEdge) {
					Edge temp = e.leftFrontNeighbor;
					e.leftFrontNeighbor = e.rightFrontNeighbor;
					e.rightFrontNeighbor = temp;
					boolean btemp = e.leftSide;
					e.leftSide = e.rightSide;
					e.rightSide = btemp;
				}

				if (e.element1 instanceof Quad) {
					Quad q = (Quad) e.element1;
					if (e == q.edgeList[base]) {
						q.updateLR();
					}
				}

				if (e.element2 instanceof Quad) {
					Quad q = (Quad) e.element2;
					if (e == q.edgeList[base]) {
						q.updateLR();
					}
				}
			}
		}
	}

	/** Change the position of the Vertex to position (x,y) */
	public void moveToPos(double x, double y) {
		this.x = x;
		this.y = y;

		updateLRinEdgeList();
	}

	/** Change the position of the Vertex to the position of the specified Vertex */
	public void moveTo(Vertex n) {
		x = n.x;
		y = n.y;

		updateLRinEdgeList();
	}

	/** Update all lengths of edges around this Vertex */
	public void updateEdgeLengths() {
		for (Edge e : edgeList) {
			e.len = e.computeLength();
		}
	}

	/**
	 * Update (almost) all angles in all elements adjacent to this Vertex. NOTE
	 * still experimental, not tested thoroughly.
	 */
	public void updateAngles() {
		Edge e, ne;
		Vertex other1, other2;
		Quad q;
		ArrayList<Element> list = new ArrayList<>();

		for (Edge element : edgeList) {
			e = element;
			if (!list.contains(e.element1)) {
				list.add(e.element1);
				ne = e.element1.neighborEdge(this, e);
				other1 = e.otherVertex(this);
				other2 = ne.otherVertex(this);

				if (e.element1 instanceof Triangle) {
					e.element1.updateAngles();
				} else {
					q = (Quad) e.element1;
					q.updateAnglesExcept(q.oppositeVertex(this));
				}
			}
			if (e.element2 != null && !list.contains(e.element2)) {
				list.add(e.element2);
				ne = e.element2.neighborEdge(this, e);
				other1 = e.otherVertex(this);
				other2 = ne.otherVertex(this);

				if (e.element2 instanceof Triangle) {
					e.element2.updateAngles();
				} else {
					q = (Quad) e.element2;
					q.updateAnglesExcept(q.oppositeVertex(this));
				}
			}
		}
	}

	public double cross(Vertex n) {
		return x * n.y - n.x * y;
	}

	public void connectToEdge(Edge edge) {
		if (!edgeList.contains(edge)) {
			edgeList.add(edge);
		}
	}

	// Rewrite of ccwSortedEdgeList().
	// We use vector representations instead of the edges directly.
	public List<MyVector> ccwSortedVectorList() {
		Element elem, start;
		MyVector v, v0, v1;
		Edge e;
		ArrayList<MyVector> boundaryVectors = new ArrayList<>();
		ArrayList<MyVector> vectors = new ArrayList<>();
		double ang;
		for (Edge element : edgeList) {
			e = element;
			v = e.getVector(this);
			v.edge = e;

			if (e.boundaryEdge()) {
				boundaryVectors.add(v);
			} else {
				vectors.add(v);
			}
		}

		// If the edgeList contains boundary edges, then select the vector of most CW
		// of these.
		// Else select the vector of an arbitrary edge.
		// The selected vector is put into v0.
		// Sets elem to the element that is ccw to v0 around this Vertex

		if (!boundaryVectors.isEmpty()) { // this size is always 0 or 2
			v0 = boundaryVectors.get(0);
			v1 = boundaryVectors.get(1);

			elem = v0.edge.element1;
			e = elem.neighborEdge(this, v0.edge);
			v = e.getVector(this);
			v.edge = e;

			if (v0.isCWto(v)) {
				if (elem.concavityAt(this)) {
					v0 = v1;
					elem = v1.edge.element1;
				}
			} else if (!elem.concavityAt(this)) {
				v0 = v1;
				elem = v1.edge.element1;
			}
		} else {
			v0 = vectors.get(0);
			elem = v0.edge.element1;
			e = elem.neighborEdge(this, v0.edge);
			v1 = e.getVector(this);
			v1.edge = e;

			if (v0.isCWto(v1)) {
				if (elem.concavityAt(this)) {
					v0 = v1;
				}
			} else if (!elem.concavityAt(this)) {
				v0 = v1;
			}
		}

		// Sort vectors in ccw order starting with v0.
		// Uses the fact that elem initially is the element ccw to v0 around this
		// Vertex.
		ArrayList<MyVector> VS = new ArrayList<>();
		e = v0.edge;

		start = elem;
		do {
			v = e.getVector(this);
			v.edge = e;
			VS.add(v);

			e = elem.neighborEdge(this, e);
			elem = elem.neighbor(e);
		} while (elem != start && elem != null);

		if (elem == null) {
			v = e.getVector(this);
			v.edge = e;
			VS.add(v);
		}

		return VS;
	}
	/*
	 * // b1: First boundary edge // b2: Second boundary edge public ArrayList
	 * ccwSortedVectorList(Edge b0, Edge b1) {
	 * Msg.debug("Entering Vertex.ccwSortedVectorList(Edge b0, Edge b1)");
	 * Msg.debug("b0: "+b0.descr()); Msg.debug("b1: "+b1.descr()); Element elem,
	 * start; MyVector v, v0, v1; Edge e; ArrayList vectors= new ArrayList();
	 * 
	 * for (int i= 0; i< edgeList.size(); i++) { e= edgeList.get(i); if (e!= b0 &&
	 * e!= b1) { v= e.getVector(this); v.edge= e; vectors.add(v); } }
	 * 
	 * v0= b0.getVector(this); v0.edge= b0; v1= b1.getVector(this); v1.edge= b1;
	 * 
	 * elem= b0.element1; if (v0.isCWto(v1)) { if (elem.concavityAt(this)) { v0= v1;
	 * elem= b1.element1; } } else if (!elem.concavityAt(this)) { v0= v1; elem=
	 * b1.element1; }
	 * Msg.debug("Vertex.ccwSortedVectorList(Edge, Edge): 0: "+v0.edge.descr());
	 * 
	 * // Sort vectors in ccw order starting with v0. // Uses the fact that elem
	 * initially is the element ccw to v0 around this Vertex. ArrayList VS= new
	 * ArrayList(); e= v0.edge;
	 * 
	 * start= elem; do { v= e.getVector(this); v.edge= e;
	 * Msg.debug("... VS.add("+v.descr()+")"); VS.add(v);
	 * 
	 * e= elem.neighborEdge(this, e); elem= elem.neighbor(e); } while (elem!= start
	 * && elem!= null);
	 * 
	 * return VS; }
	 */

	/**
	 * Assumes: b0 and b1 form a convex boundary, but not neccessarily *strictly*
	 * convex
	 * 
	 * @param b1 First boundary edge
	 * @param b2 Second boundary edge
	 */
	public List<Edge> calcCCWSortedEdgeList(Edge b0, Edge b1) {
		MyVector v, v0, v1;
		Edge e;
		ArrayList<MyVector> vectors = new ArrayList<>();

		for (Edge element : edgeList) {
			e = element;
			if (e != b0 && e != b1) {
				v = e.getVector(this);
				v.edge = e;
				vectors.add(v);
			}
		}

		// Initially put the two vectors of b0 and b1 in list.
		// Select the most CW boundary edge to be first in list.

		ArrayList<MyVector> VS = new ArrayList<>();
		v0 = b0.getVector(this);
		v0.edge = b0;
		v1 = b1.getVector(this);
		v1.edge = b1;

		if (vectors.size() > 0) {
			v = vectors.get(0);
		} else {
			v = v1;
		}

		if (v0.isCWto(v)) {
			VS.add(v0);
			VS.add(v1);
		} else {
			VS.add(v1);
			VS.add(v0);
		}

		// Sort vectors in ccw order. I will not move the vector that lies first in VS.
		for (int i = 0; i < vectors.size(); i++) {
			v = vectors.get(i);

			for (int j = 0; j < VS.size(); j++) {
				v0 = VS.get(j);
				if (j + 1 == VS.size()) {
					v1 = VS.get(0);
				} else {
					v1 = VS.get(j + 1);
				}

				if (!v.isCWto(v0) && v.isCWto(v1)) {
					VS.add(j + 1, v);
					break;
				}
			}
		}

		List<Edge> edges = new ArrayList<>(VS.size());
		for (MyVector element : VS) {
			v = element;
			edges.add(v.edge);
		}
		return edges;
	}

	/**
	 * NOTE *ALL* vertices in a neighboring quad is regarded as neighbors, not only
	 * those that are directly connected to this Vertex by edges.
	 * 
	 * @return a ccw sorted list of the neighboring vertices to this, but returns
	 *         null if this Vertex is part of any triangle.
	 */
	public Vertex[] ccwSortedNeighbors() {
		Element elem;
		MyVector v, v0, v1;
		Edge e;

		// First try to find two boundary edges
		int j = 0;
		MyVector[] b = new MyVector[2];
		for (Edge edge : edgeList) {
			e = edge;
			if (e.boundaryEdge()) {
				b[j] = e.getVector(this);
				b[j++].edge = e;
				if (j == 2) {
					break;
				}
			}
		}

		// If these are found, then v0 is the vector of the most cw edge.
		if (j == 2) {
			elem = b[0].edge.element1;
			e = elem.neighborEdge(this, b[0].edge);
			v1 = e.getVector(this);
			v1.edge = e;

			if (b[0].isCWto(v1)) {
				v0 = b[0];
			} else {
				v0 = b[1]; // that is, the other boundary vector
				elem = b[1].edge.element1;
			}
		} else {
			/*
			 * Failing to find any boundary edges, we select the vector of an arbitrary edge
			 * to be v0. Sets elem to the element that is ccw to v0 around this Vertex.
			 */
			e = edgeList.get(0);
			v0 = e.getVector(this);
			v0.edge = e;
			elem = e.element1;
			e = elem.neighborEdge(this, e);
			v1 = e.getVector(this);
			v1.edge = e;

			if (v0.isCWto(v1)) {
				v0 = v0;
			} else {
				v0 = v1;
			}
		}

		/*
		 * Sort vertices in ccw order starting with otherVertex of v0 edge. Uses the
		 * fact that elem initially is the element ccw to v0 around this Vertex.
		 */
		Vertex[] ccwvertexList = new Vertex[edgeList.size() * 2];
		Element start = elem;
		Quad q;
		e = v0.edge;
		int i = 0;
		do {
			ccwvertexList[i++] = e.otherVertex(this);
			if (!(elem instanceof Quad)) {
				return null;
			}
			q = (Quad) elem;
			ccwvertexList[i++] = q.oppositeVertex(this);
			e = elem.neighborEdge(this, e);
			elem = elem.neighbor(e);
		} while (elem != start && elem != null);

		if (elem == null) {
			ccwvertexList[i++] = e.otherVertex(this);
		}

		return ccwvertexList;
	}

	public double meanNeighborEdgeLength() {
		double sumLengths = 0.0, len, j = 0;

		for (Edge e : edgeList) {
			len = e.length();
			if (len != 0) {
				j++;
				sumLengths += len;
			}
		}
		return sumLengths / j;
	}

	public int nrOfAdjElements() {
		List<Element> list = adjElements();
		return list.size();
	}

	public List<Element> adjElements() {
		Edge e;
		List<Element> list = new ArrayList<>();

		for (int i = 0; i < edgeList.size(); i++) {
			e = edgeList.get(i);
			if (!list.contains(e.element1)) {
				list.add(e.element1);
			}
			if (e.element2 != null && !list.contains(e.element2)) {
				list.add(e.element2);
			}
		}
		return list;
	}

	public int nrOfAdjQuads() {
		List<Element> list = adjQuads();
		return list.size();
	}

	public List<Element> adjQuads() {
		Edge e;
		ArrayList<Element> list = new ArrayList<>();

		for (Edge element : edgeList) {
			e = element;
			if (e.element1 instanceof Quad && !list.contains(e.element1)) {
				list.add(e.element1);
			} else if (e.element2 != null && e.element2 instanceof Quad && !list.contains(e.element2)) {
				list.add(e.element2);
			}
		}
		return list;
	}

	public int nrOfAdjTriangles() {
		List<Element> list = adjTriangles();
		return list.size();
	}

	// Hmm. Should I include fake quads as well?
	public List<Element> adjTriangles() {
		Edge e;
		ArrayList<Element> list = new ArrayList<>();

		for (Edge element : edgeList) {
			e = element;
			if (e.element1 instanceof Triangle && !list.contains(e.element1)) {
				list.add(e.element1);
			} else if (e.element2 != null && e.element2 instanceof Triangle && !list.contains(e.element2)) {
				list.add(e.element2);
			}
		}
		return list;
	}

	/**
	 * Classic Laplacian smooth. Of course, to be run on internal vertices only.
	 * 
	 * @return the vector from the old to the new position.
	 */
	public MyVector laplacianMoveVector() {
		MyVector c, cJSum = new MyVector(origin, origin);
		Edge e;
		Vertex nJ;

		int n = edgeList.size();
		for (int i = 0; i < n; i++) {
			e = edgeList.get(i);
			nJ = e.otherVertex(this);
			c = new MyVector(this, nJ);
			cJSum = cJSum.plus(c);
		}
		cJSum = cJSum.div(n);
		return cJSum;
	}

	/**
	 * Classic Laplacian smooth. Of course, to be run on internal vertices only.
	 * 
	 * @return the new position of Vertex
	 */
	public Vertex laplacianSmooth() {
		MyVector c, cJSum = new MyVector(origin, origin);
		Edge e;
		Vertex nJ;

		int n = edgeList.size();
		for (int i = 0; i < n; i++) {
			e = edgeList.get(i);
			nJ = e.otherVertex(this);
			c = new MyVector(this, nJ);
			cJSum = cJSum.plus(c);
		}
		cJSum = cJSum.div(n);
		return new Vertex(x + cJSum.x, y + cJSum.y);
	}

	/**
	 * Classic Laplacian smooth, but exclude the given neighbor Vertex from the
	 * calculation. Of course, to be run on internal vertices only.
	 * 
	 * @param vertex the Vertex to be excluded
	 * @return the new position of Vertex
	 */
	public Vertex laplacianSmoothExclude(Vertex vertex) {
		MyVector c, cJSum = new MyVector(origin, origin);
		Edge e;
		Vertex nJ;

		int n = edgeList.size();
		for (int i = 0; i < n; i++) {
			e = edgeList.get(i);
			nJ = e.otherVertex(this);
			if (nJ != vertex) {
				c = new MyVector(this, nJ);
				cJSum = cJSum.plus(c);
			}
		}
		cJSum = cJSum.div(n - 1); // -1 because Vertex is excluded
		return new Vertex(x + cJSum.x, y + cJSum.y);
	}

	/**
	 * Run this on internal vertices (not part of the boundary or front) Does a
	 * modified length weighted Laplacian smooth.
	 * 
	 * @return a new Vertex with the smoothed position.
	 */
	public Vertex modifiedLWLaplacianSmooth() {
		Vertex nJ;
		double cJLengthSum = 0, len;
		Edge e, bEdge1, bEdge2;
		MyVector c, cJLengthMulcJSum = new MyVector(origin, origin), deltaCj, deltaI;
		int n = edgeList.size();
		if (n == 0) {
			Msg.error("...edgeList.size()== 0");
		}
		for (int i = 0; i < n; i++) {
			e = edgeList.get(i);
			nJ = e.otherVertex(this);
			c = new MyVector(this, nJ);
			if (nJ.boundaryVertex()) {
				bEdge1 = nJ.anotherBoundaryEdge(null);
				bEdge2 = nJ.anotherBoundaryEdge(bEdge1);

				// This should be correct:
				deltaCj = nJ.angularSmoothnessAdjustment(this, bEdge1, bEdge2, e.length());
				c = c.plus(deltaCj);
			}

			len = c.length();
			c = c.mul(len);
			cJLengthMulcJSum = cJLengthMulcJSum.plus(c);
			cJLengthSum += len;
		}
		deltaI = cJLengthMulcJSum.div(cJLengthSum);

		return new Vertex(x + deltaI.x, y + deltaI.y);
	}

	public int nrOfFrontEdges() {
		int fronts = 0;
		for (Edge e : edgeList) {
			if (e.frontEdge) {
				fronts++;
			}
		}
		return fronts;
	}

	/**
	 * An implementation of an algorithm described in a paper by Blacker &
	 * Stephenson.
	 * 
	 * @param nJ     the other Vertex that lies behind this Vertex (not on the
	 *               front/boundary)
	 * @param ld     length from this to nJ
	 * @param front1 front/boundary neighbor edge to this
	 * @param front2 front/boundary neighbor edge to this
	 * @return a new Vertex (with a smoothed positing) that can replace this Vertex.
	 */
	public Vertex blackerSmooth(Vertex nJ, Edge front1, Edge front2, double ld) {
		Vertex nI = this;
		Vertex origin = new Vertex(0, 0);
		Vertex n1, n2, n3, n4;
		Quad q;
		List<Element> adjQuads = adjQuads();

		// Step 1, the isoparametric smooth:
		MyVector vI = new MyVector(origin, nI);
		MyVector vMXsum = new MyVector(origin, origin);
		MyVector vMJ;
		MyVector vMK;
		MyVector vML;

		for (Element e : adjQuads) {
			q = (Quad) e;

			n1 = q.edgeList[base].leftVertex;
			n2 = q.edgeList[base].rightVertex;
			n3 = q.edgeList[left].otherVertex(q.edgeList[base].leftVertex);
			n4 = q.edgeList[right].otherVertex(q.edgeList[base].rightVertex);

			// Sorting vMJ, vMK, and vML in ccw order:
			if (nI == n1) {
				vMJ = new MyVector(origin, n2);
				vMK = new MyVector(origin, n4);
				vML = new MyVector(origin, n3);
			} else if (nI == n2) {
				vMJ = new MyVector(origin, n4);
				vMK = new MyVector(origin, n3);
				vML = new MyVector(origin, n1);
			} else if (nI == n3) {
				vMJ = new MyVector(origin, n1);
				vMK = new MyVector(origin, n2);
				vML = new MyVector(origin, n4);
			} else { // if (nI==n4) {
				vMJ = new MyVector(origin, n3);
				vMK = new MyVector(origin, n1);
				vML = new MyVector(origin, n2);
			}

			vMXsum = vMXsum.plus(vMJ);
			vMXsum = vMXsum.plus(vML);
			vMXsum = vMXsum.minus(vMK);
		}

		MyVector vImarked = vMXsum.div(adjQuads.size());
		MyVector deltaA = vImarked.minus(vI);

		if (adjQuads.size() != 2 || nrOfFrontEdges() > 2) {
			return new Vertex(x + deltaA.x, y + deltaA.y);
		}
		// Step 2, length adjustment:
		else {
			MyVector vJ = new MyVector(origin, nJ);
			MyVector vIJ = new MyVector(nJ, vImarked.x, vImarked.y);
			double la = vIJ.length();

			MyVector deltaB = deltaA.plus(vI);
			deltaB = deltaB.minus(vJ);
			deltaB = deltaB.mul(ld / la);
			deltaB = deltaB.plus(vJ);
			deltaB = deltaB.minus(vI);

			// Step 3, angular smoothness:
			MyVector deltaC = angularSmoothnessAdjustment(nJ, front1, front2, ld);
			MyVector deltaI = deltaB.plus(deltaC);
			deltaI = deltaI.mul(0.5);
			return new Vertex(x + deltaI.x, y + deltaI.y);
		}
	}

	/**
	 * Performs an angular smoothness adjustment as described in the paper by
	 * Blacker and Stephenson. Assumes that this is a Vertex that lies on the
	 * front/boundary.
	 * 
	 * @param nJ the Vertex connected to this, that lies behind the front/boundary
	 * @param f1 front/boundary neighbor edge to this
	 * @param f2 front/boundary neighbor edge to this
	 * @return a vector that should replace the edge between this and nJ
	 */
	public MyVector angularSmoothnessAdjustment(Vertex nJ, Edge f1, Edge f2, double ld) {
		Vertex nI = this;
		if (Double.isNaN(ld)) {
			Msg.error("ld is NaN!!!");
		}

		if (f2.length() == 0) {
			Msg.error("f2.length()== 0");
		}

		Vertex nIm1 = f1.otherVertex(nI);
		Vertex nIp1 = f2.otherVertex(nI);

		if (nIm1.equals(nI)) {
			Msg.error("nIm1.equals(nI)");
		}

		if (nIp1.equals(nI)) {
			Msg.error("nIp1.equals(nI)");
		}

		// if (nIm1.equals(nIp1))
		// this should be okay, in fact...
		// Msg.error("nIm1.equals(nIp1)");

		MyVector pI1 = new MyVector(nJ, nIm1);
		MyVector pI = new MyVector(nJ, nI);
		MyVector pI2 = new MyVector(nJ, nIp1);

		double pI1Angle = pI1.posAngle();
		double pI2Angle = pI2.posAngle();
		double pIp1Angle = Math.max(pI1Angle, pI2Angle);
		double pIm1Angle = Math.min(pI1Angle, pI2Angle);
		double pIAngle = pI.posAngle();
		double pIm1p1Angle = 0;
		if (pIAngle < pIm1Angle || pIAngle > pIp1Angle) {
			pIm1p1Angle = TWO_PI - pIp1Angle + pIm1Angle;
		} else {
			pIm1p1Angle = pIp1Angle - pIm1Angle;
		}

		// Check if the sum of angles between pIp1 and pI and the angle between pIm1 and
		// PI is greater or equal to 180 degrees. If so, I choose ld as the length of
		// pB2.
		if (pIm1p1Angle > Math.PI) {
			double pB1Angle = pIm1p1Angle * 0.5 + pIp1Angle;
			if (pB1Angle >= TWO_PI) {
				pB1Angle = Math.IEEEremainder(pB1Angle, TWO_PI);
			}
			double pB1pIMax = Math.max(pB1Angle, pIAngle);
			double pB1pIMin = Math.min(pB1Angle, pIAngle);
			double pB2Angle = pB1pIMin + 0.5 * (pB1pIMax - pB1pIMin);
			if (pB1pIMax - pB1pIMin > Math.PI) {
				pB2Angle += Math.PI;
			}

			MyVector pB2 = new MyVector(pB2Angle, ld, nJ);
			MyVector deltaC = pB2.minus(pI);
			return deltaC;
		}

		MyVector line = new MyVector(nIp1, nIm1);

		// pB1 should be the halved angle between pIp1 and pIm1, in the direction of pI:
		double pB1Angle = pIm1Angle + 0.5 * (pIp1Angle - pIm1Angle);
		if (pIp1Angle - pIm1Angle > Math.PI) {
			pB1Angle += Math.PI;
		}

		if (Double.isNaN(pB1Angle)) {
			Msg.error("pB1Angle is NaN!!!");
		}
		double pB1pIMax = Math.max(pB1Angle, pIAngle);
		double pB1pIMin = Math.min(pB1Angle, pIAngle);
		double pB2Angle = pB1pIMin + 0.5 * (pB1pIMax - pB1pIMin);
		if (pB1pIMax - pB1pIMin > Math.PI) {
			pB2Angle += Math.PI;
		}

		if (Double.isNaN(pB2Angle)) {
			Msg.error("pB2Angle is NaN!!!");
		}
		Ray pB2Ray = new Ray(nJ, pB2Angle);
		// MyVector pB2= new MyVector(pB2Angle, 100.0, nJ);

		Vertex q = pB2Ray.pointIntersectsAt(line);
		double lq = q.length(nJ);
		if (Double.isNaN(lq)) {
			Msg.error("lq is NaN!!!");
		}

		MyVector pB2;
		if (ld > lq) {
			pB2 = new MyVector(pB2Ray, (lq + ld) * 0.5);
			// pB2.setLengthAndAngle((lq+ld)*0.5, pB2Angle);
		} else {
			pB2 = new MyVector(pB2Ray, ld);
			// pB2.setLengthAndAngle(ld, pB2Angle);
		}

		MyVector deltaC = pB2.minus(pI);
		return deltaC;
	}

	/**
	 * Test whether any of the adjacent elements has become inverted or their areas
	 * are zero.
	 * 
	 * @param elements the list of elements to parse
	 * @return true if the movement of a Vertex has caused any of it's adjacent
	 *         elements to become inverted or get an area of size zero.
	 */
	public boolean invertedOrZeroAreaElements(List<Element> elements) {
		for (Element elem : elements) {
			if (elem.invertedOrZeroArea()) {
				return true;
			}
		}
		return false;
	}

	/**
	 * Incrementally adjust the location of the Vertex (along a vector) until none
	 * of it's neighboring elements are inverted. Use increments of size vector
	 * component divided by 50 in each direction, unless ONE of these increments is
	 * less than a given lower limit. If so, the increments in the direction of the
	 * shortest component should be equal to that limit, while the other direction
	 * is dictated by the first, of course.
	 * 
	 * @return true on success else false.
	 */
	public boolean incrAdjustUntilNotInvertedOrZeroArea(Vertex old, List<Element> elements) {
		MyVector back = new MyVector(this, old);
		double startX = x, startY = y;
		double xstep = back.x / 50.0, ystep = back.y / 50.0;
		double xinc, yinc;
		int steps, i;

		if (Math.abs(xstep) < COINCTOL || Math.abs(ystep) < COINCTOL) {

			if (COINCTOL < Math.abs(back.x) && COINCTOL < Math.abs(back.y)) {// && or || ?
				if (Math.abs(back.x) < Math.abs(back.y)) {
					if (back.x < 0) {
						xinc = -COINCTOL;
					} else {
						xinc = COINCTOL;
					}

					yinc = Math.abs(old.y) * COINCTOL / Math.abs(old.x);
					if (back.y < 0) {
						yinc = -yinc;
					}

					steps = (int) (back.x / xinc);
				} else {
					if (back.y < 0) {
						yinc = -COINCTOL;
					} else {
						yinc = COINCTOL;
					}

					xinc = Math.abs(old.x) * COINCTOL / Math.abs(old.y);
					if (back.x < 0) {
						xinc = -xinc;
					}

					steps = (int) (back.y / yinc);
				}

				for (i = 1; i <= steps; i++) {
					x = startX + xinc * i;
					y = startY + yinc * i;

					if (!invertedOrZeroAreaElements(elements)) {
						return true;
					}
				}
			}
		} else {
			for (i = 1; i <= 49; i++) {
				x = startX + back.x * i / 50.0;
				y = startY + back.y * i / 50.0;

				if (!invertedOrZeroAreaElements(elements)) {
					return true;
				}
			}
		}

		x = old.x;
		y = old.y;
		if (!invertedOrZeroAreaElements(elements)) {
			return true;
		}

		return false;
	}

	/** @return true if the Vertex is part of the boundary of the mesh. */
	public boolean boundaryVertex() {
		for (Edge e : edgeList) {
			if (e.boundaryEdge()) {
				return true;
			}
		}
		return false;
	}

	/**
	 * @return true if the Vertex is part of the boundary of the mesh or a triangle.
	 */
	public boolean boundaryOrTriangleVertex() {
		for (Edge e : edgeList) {
			if (e.boundaryOrTriangleEdge()) {
				return true;
			}
		}
		return false;
	}

	/** @return true if the Vertex is truely a part of the front. */
	public boolean frontVertex() {
		for (Edge e : edgeList) {
			if (e.isFrontEdge()) {
				return true;
			}
		}
		return false;
	}

	/**
	 * @param known is the front edge that is already known. (Use null if no such
	 *              edge is known.)
	 * @return a front edge found in this Vertex's edgelist.
	 */
	public Edge anotherFrontEdge(Edge known) {
		for (Edge e : edgeList) {
			if (e != known && e.isFrontEdge()) {
				return e;
			}
		}
		return null;
	}

	/**
	 * @param known is the boundary edge that is already known. (Use null if no such
	 *              edge is known.)
	 * @return a boundary edge found in this Vertex's edgelist.
	 */
	public Edge anotherBoundaryEdge(Edge known) {
		for (Edge e : edgeList) {
			if (e != known && e.boundaryEdge()) {
				return e;
			}
		}
		return null;
	}

	public double length(Vertex n) {
		double xDiff = x - n.x;
		double yDiff = y - n.y;
		return Math.sqrt(xDiff * xDiff + yDiff * yDiff);
	}

	public double length(double x, double y) {
		double xDiff = this.x - x;
		double yDiff = this.y - y;
		return Math.sqrt(xDiff * xDiff + yDiff * yDiff);
	}

	/**
	 * Determine if a Vertex is on the line (of infinite length) that e is a part
	 * of.
	 */
	public boolean onLine(Edge e) {
		BigDecimal x1 = new BigDecimal(e.leftVertex.x);
		BigDecimal y1 = new BigDecimal(e.leftVertex.y);
		BigDecimal x2 = new BigDecimal(e.rightVertex.x);
		BigDecimal y2 = new BigDecimal(e.rightVertex.y);
		BigDecimal x3 = new BigDecimal(x);
		BigDecimal y3 = new BigDecimal(y);

		BigDecimal zero = new BigDecimal(0.0);
		BigDecimal l_cross_r = (x1.multiply(y2)).subtract(x2.multiply(y1));
		BigDecimal xdiff = x1.subtract(x2);
		BigDecimal ydiff = y1.subtract(y2);
		BigDecimal det1 = l_cross_r.subtract(xdiff.multiply(y3)).add(ydiff.multiply(x3));

		int eval1 = det1.compareTo(zero);
		if (eval1 == 0) {
			return true;
		} else {
			return false;
		}
	}

	/**
	 * Determine if a Vertex is in a given halfplane. The method is based on the
	 * determinant as described in Schewchuk's paper.
	 * 
	 * @return 1 if this Vertex is in the halfplane defined by Triangle t and Edge
	 *         e, 0 if the Vertex is on Edge e, -1 if the Vertex is not in the
	 *         halfplane defined by Triangle t and Edge e.
	 */
	public int inHalfplane(Triangle t, Edge e) {
		return inHalfplane(e.leftVertex, e.rightVertex, t.oppositeOfEdge(e));
	}

	// @return 1 if this Vertex is on the same side of Edge e as Vertex n is,
	// 0 if this Vertex is on the line that extends Edge e, and
	// -1 if this Vertex is on the other side of Edge e than Vertex n is.
	public int inHalfplane(Edge e, Vertex n) {
		return inHalfplane(e.leftVertex, e.rightVertex, n);
	}

	// @return 1 if this Vertex is on the same side of the line (l1, l2) as Vertex n
	// is,
	// 0 if this Vertex is on the line that extends line (l1, l2), and
	// -1 if this Vertex is on the other side of line (l1, l2) than Vertex n is.
	public int inHalfplane(Vertex l1, Vertex l2, Vertex n) {
		BigDecimal x1 = new BigDecimal(l1.x);
		BigDecimal y1 = new BigDecimal(l1.y);
		BigDecimal x2 = new BigDecimal(l2.x);
		BigDecimal y2 = new BigDecimal(l2.y);
		BigDecimal x3 = new BigDecimal(x);
		BigDecimal y3 = new BigDecimal(y);
		BigDecimal x4 = new BigDecimal(n.x);
		BigDecimal y4 = new BigDecimal(n.y);

		BigDecimal zero = new BigDecimal(0.0);
		BigDecimal l_cross_r = (x1.multiply(y2)).subtract(x2.multiply(y1));
		BigDecimal xdiff = x1.subtract(x2);
		BigDecimal ydiff = y1.subtract(y2);
		BigDecimal det1 = l_cross_r.subtract(xdiff.multiply(y3)).add(ydiff.multiply(x3));

		int eval1 = det1.compareTo(zero);
		if (eval1 == 0) {
			return 0;
		}

		BigDecimal det2 = l_cross_r.subtract(xdiff.multiply(y4)).add(ydiff.multiply(x4));
		int eval2 = det2.compareTo(zero);
		if ((eval1 < 0 && eval2 < 0) || (eval1 > 0 && eval2 > 0)) {
			return 1;
		} else {
			return -1;
		}
	}

	/**
	 * Test to see if this Vertex lies in the plane bounded by the two parallel
	 * lines intersecting the Vertices of Edge e that are normal to Edge e.
	 */
	public boolean inBoundedPlane(Edge e) {
		Edge normal1 = e.unitNormalAt(e.leftVertex);
		Edge normal2 = e.unitNormalAt(e.rightVertex);

		int a = inHalfplane(normal1, e.rightVertex);
		int b = inHalfplane(normal2, e.leftVertex);

		if ((a == 1 || a == 0) && (b == 1 || b == 0)) {
			return true;
		} else {
			return false;
		}
	}

	/**
	 * Return true if the circle intersecting the Vertices p1, p2, and p3 contains
	 * this Vertex in its interior. p1, p2, p3, and p4 are ccw sorted. Note that
	 * testing for convexity of the quad should not be necessary.
	 */
	public boolean inCircle(Vertex p1, Vertex p2, Vertex p3) {

		/*
		 * a * b = a_1*b_1 + a_2*b_2 Scalar product a x b = (a_1*b_2 - b_1*a_2) Vector
		 * product
		 * 
		 * e_3 is the unit vector (0,0,1). All the points are defined in 3D space, with
		 * z-values equal to 0. The v_i vectors are unit vectors. The points are ccw
		 * ordered.
		 * 
		 * sin(alpha) = (v_1 x v_2) * e_3= x3y1 -x3y2 -x2y1 -y3x1 +y3x2 +y2x1 sin(beta)
		 * = (v_3 x v_4) * e_3= y3x1 -x1y4 -x4y3 -x3y1 +y1x4 +y4x3 cos(alpha) = v_1 *
		 * v_2 = (x3 -x2)(x1 -x2) +(y3 -y2)(y1 -y2) cos(beta) = v_3 * v_4 = (x1 -x4)(x3
		 * -x4) +(y1 -y4)(y3 -y4)
		 * 
		 */

		double cosAlpha = (p3.x - p2.x) * (p1.x - p2.x) + (p3.y - p2.y) * (p1.y - p2.y);
		double cosBeta = (p1.x - x) * (p3.x - x) + (p1.y - y) * (p3.y - y);

		if (cosAlpha < 0 && cosBeta < 0) { // (if both angles > than 90 degrees)
			return true;
		} else if (cosAlpha > 0 && cosBeta > 0) { // (if both angles < than 90 degrees)
			return false;
		} else {
			double sinAlpha = p3.x * p1.y - p3.x * p2.y - p2.x * p1.y - p3.y * p1.x + p3.y * p2.x + p2.y * p1.x;
			double sinBeta = p3.y * p1.x - p1.x * y - x * p3.y - p3.x * p1.y + p1.y * x + y * p3.x;
			if (cosAlpha * sinBeta + sinAlpha * cosBeta < 0) {
				return true;
			} else {
				return false;
			}
		}
	}

	/**
	 * Pretending this and n has the same location, copy the edges in n's edgelist
	 * that this Vertex doesn't already have, and put them into this Vertex's
	 * edgeList. If this and n have any common edges, these must be removed.
	 */
	public void merge(Vertex n) {
		Vertex oldN = n.copyXY();
		Edge e;
		int ind;
		n.setXY(this);
		for (int i = 0; i < n.edgeList.size(); i++) {
			e = n.edgeList.get(i);
			ind = edgeList.indexOf(e);
			if (ind == -1) {
				e.replaceVertex(n, this);
				edgeList.add(e);
			} else if (e.leftVertex == e.rightVertex) {
				edgeList.remove(ind);
			}
		}
		n.setXY(oldN);
	}

	public List<Edge> frontEdgeList() {
		ArrayList<Edge> list = new ArrayList<>();
		Edge e;
		for (Edge element : edgeList) {
			e = element;
			if (e.frontEdge) {
				list.add(e);
			}
		}
		return list;
	}

	/**
	 * Parse the edgeList to look for Edge e.
	 * 
	 * @return true if found, else false
	 */
	public boolean hasEdge(Edge e) {
		for (Edge curEdge : edgeList) {
			if (e == curEdge) {
				return true;
			}
		}
		return false;
	}

	/** Determine the Vertex valence. */
	public byte valence() {
		byte temp = (byte) edgeList.size();
		if (!boundaryVertex()) {
			return temp;
		} else {
			Edge b1 = anotherBoundaryEdge(null);
			Edge b2 = anotherBoundaryEdge(b1);
			double ang = b1.sumAngle(b1.element1, this, b2);

			// Determine which kind of boundary Vertex we're dealing with
			if (ang <= PIx3div4) {
				return (byte) (temp + 2);
			} else if (ang < PIx5div4) {
				return (byte) (temp + 1);
			} else {
				return (temp);
			}

			/*
			 * if (Math.abs(ang- PIdiv2) < Math.abs(ang- Math.PI)) return (byte)(temp+2);
			 * else if (Math.abs(ang- Math.PI) < Math.abs(ang- TWO_PI)) return
			 * (byte)(temp+1); else return temp;
			 */
		}
	}

	/** Calculate the valence pattern for this Vertex and its neighbors. */
	public void createValencePattern(Vertex[] ccwVertices) {
		int j = edgeList.size() * 2;
		if (j >= 128) {
			Msg.error("Number of edges adjacent Vertex " + descr() + " was greater than expected (" + edgeList.size() + "-2 >= 64)");
		}
		byte ccwVerticesSize = (byte) j;
		pattern = new byte[ccwVerticesSize + 2]; // +2 for size and c.valence()
		pattern[0] = (byte) (ccwVerticesSize + 2);
		pattern[1] = valence();

		for (int i = 0; i < ccwVerticesSize; i++) {
			pattern[i + 2] = ccwVertices[i].valence();
		}
	}

	/** Calculate the valence pattern for this Vertex and its neighbors. */
	public void createValencePattern(byte ccwVerticesSize, Vertex[] ccwVertices) {
		pattern = new byte[ccwVerticesSize + 2]; // +2 for size and c.valence()
		pattern[0] = (byte) (ccwVerticesSize + 2);
		pattern[1] = valence();

		for (int i = 0; i < ccwVerticesSize; i++) {
			pattern[i + 2] = ccwVertices[i].valence();
		}
	}

	/**
	 * Return # of irregular vertices in the valence pattern (vertices whose
	 * valence!= 4) Note that calcMyValencePattern() must be called before calling
	 * this method.
	 */
	public int irregNeighborVertices() {
		int count = 0;
		for (int i = 1; i < pattern[0]; i++) {
			if (pattern[i] != 4) {
				count++;
			}
		}
		return count;
	}

	/**
	 * Compare the valence pattern of this Vertex to the special pattern in
	 * pattern2. In pattern2, the following codes apply:<br>
	 * <ul>
	 * 14 means 4- (4 or less)
	 * <li>24 means 4+ (4 or more)
	 * <li>5 means 5 or more
	 * <li>0 means any number
	 * </ul>
	 * Note that the length of the patterns must be an even number. Also note that
	 * the patterns have to be aligned in a 2-byte fashion. (This means that the
	 * index into the Vertex's pattern where they start matching have to be an even
	 * number.)
	 * 
	 * @param pattern2 A valence pattern
	 * @return If they match then return the index of the valence value in the
	 *         Vertex's pattern that corresponds to the first valence value in
	 *         pattern2, otherwise return -1. Note that calcMyValencePattern() must
	 *         be called before calling this method.
	 */
	public int patternMatch(byte[] pattern2) {
		if (pattern[0] != pattern2[0] || pattern[1] != pattern2[1]) {
			return -1; // Different length or different valence of central Vertex
		}

		int i, j, jstart = 2, matches = 0;

		while (jstart < pattern[0]) {
			// Find index of next valence in pattern2 that matches valence of pattern[2]
			for (j = jstart; j < pattern[0]; j += 2) {
				if (pattern[j] == 2 && (pattern2[2] == 2 || pattern2[2] == 14 || pattern2[2] == 0)) {
					matches = 1;
					jstart = j;
					break;
				} else if (pattern[j] == 3 && (pattern2[2] == 3 || pattern2[2] == 14 || pattern2[2] == 0)) {

					matches = 1;
					jstart = j;
					break;
				} else if (pattern[j] == 4 && (pattern2[2] == 4 || pattern2[2] == 14 || pattern2[2] == 24 || pattern2[2] == 0)) {

					matches = 1;
					jstart = j;
					break;
				} else if (pattern[j] >= 5 && (pattern2[2] == 5 || pattern2[2] == 24 || pattern2[2] == 0)) {

					matches = 1;
					jstart = j;
					break;
				}
			}

			if (matches == 0) {
				return -1; // Search completed, patterns don't match
			}

			if (jstart == pattern[0] - 1) {
				j = 1;
			} else {
				j = jstart + 1;
			}
			// Count nr of sequential matches starting at this index
			for (i = 3; matches < pattern2[0] - 2; i++, j++) {
				if (pattern[j] == 2 && (pattern2[i] == 2 || pattern2[i] == 14 || pattern2[i] == 0)) {
					matches++;
				} else if (pattern[j] == 3 && (pattern2[i] == 3 || pattern2[i] == 14 || pattern2[i] == 0)) {
					matches++;
				} else if (pattern[j] == 4 && (pattern2[i] == 4 || pattern2[i] == 14 || pattern2[i] == 24 || pattern2[i] == 0)) {
					matches++;
				} else if (pattern[j] >= 5 && (pattern2[i] == 5 || pattern2[i] == 24 || pattern2[i] == 0)) {
					matches++;
				} else {
					matches = 0;
					break;
				}

				/*
				 * if (pattern2[i]!= pattern[j] && pattern2[i]!= 0 && pattern[j]!= 0) { matches=
				 * 0; break; } else { matches++; Msg.debug("...pattern2["+i+"]: "+pattern2[i]+
				 * ", pattern["+j+"]: "+pattern[j]); }
				 */

				if (j == pattern[0] - 1) {
					j = 1;
				}
			}
			if (matches == pattern2[0] - 2) {
				return jstart; // Search completed, patterns match
			}
			jstart += 2;
		}
		return -1;
	}

	/**
	 * Compare the valence pattern of this Vertex to the special pattern in
	 * pattern2. Also make sure that the tagged vertices in vertexPat are vertices.
	 * (That is, the interior angles must be greater than any other interior angles
	 * around this Vertex.) In pattern2, the following codes apply:<br>
	 * <ul>
	 * 14 means 4- (4 or less)
	 * <li>24 means 4+ (4 or more)
	 * <li>5 means 5 or more
	 * <li>0 means any number
	 * </ul>
	 * Note that the length of the patterns must be an even number. Also note that
	 * the patterns have to be aligned in a 2-byte fashion. (This means that the
	 * index into the Vertex's pattern where they start matching have to be an even
	 * number.)
	 * 
	 * @param pattern2 A valence pattern
	 * @return If they match then return the index of the valence value in the
	 *         Vertex's pattern that corresponds to the first valence value in
	 *         pattern2, otherwise return -1. Note that calcMyValencePattern() must
	 *         be called before calling this method.
	 */
	public int patternMatch(byte[] pattern2, boolean[] vertexPat, double[] angles) {
		if (pattern[0] != pattern2[0] || pattern[1] != pattern2[1]) {
			return -1; // Different length or different valence of central Vertex
		}

		int i, j, jstart = 2, matches = 0;

		while (jstart < pattern[0]) {
			// Find index of next valence in pattern2 that matches valence of pattern[2]
			for (j = jstart; j < pattern[0]; j += 2) {
				if (pattern[j] == 2 && (pattern2[2] == 2 || pattern2[2] == 14 || pattern2[2] == 0)) {
					if (fitsVertexPat((byte) (j - 2), angles, vertexPat, pattern[0] - 2)) {
						matches = 1;
						jstart = j;
						break;
					}
				} else if (pattern[j] == 3 && (pattern2[2] == 3 || pattern2[2] == 14 || pattern2[2] == 0)) {
					if (fitsVertexPat((byte) (j - 2), angles, vertexPat, pattern[0] - 2)) {
						matches = 1;
						jstart = j;
						break;
					}
				} else if (pattern[j] == 4 && (pattern2[2] == 4 || pattern2[2] == 14 || pattern2[2] == 24 || pattern2[2] == 0)) {
					if (fitsVertexPat((byte) (j - 2), angles, vertexPat, pattern[0] - 2)) {
						matches = 1;
						jstart = j;
						break;
					}
				} else if (pattern[j] >= 5 && (pattern2[2] == 5 || pattern2[2] == 24 || pattern2[2] == 0)) {
					if (fitsVertexPat((byte) (j - 2), angles, vertexPat, pattern[0] - 2)) {
						matches = 1;
						jstart = j;
						break;
					}
				}
			}

			if (matches == 0) {
				return -1; // Search completed, patterns don't match
			}
			if (jstart == pattern[0] - 1) {
				j = 1; // Shouldn't it be 2???
			} else {
				j = jstart + 1;
			}
			// Count nr of sequential matches starting at this index
			for (i = 3; matches < pattern2[0] - 2; i++, j++) {
				if (pattern[j] == 2 && (pattern2[i] == 2 || pattern2[i] == 14 || pattern2[i] == 0)) {
					matches++;
				} else if (pattern[j] == 3 && (pattern2[i] == 3 || pattern2[i] == 14 || pattern2[i] == 0)) {
					matches++;
				} else if (pattern[j] == 4 && (pattern2[i] == 4 || pattern2[i] == 14 || pattern2[i] == 24 || pattern2[i] == 0)) {
					matches++;
				} else if (pattern[j] >= 5 && (pattern2[i] == 5 || pattern2[i] == 24 || pattern2[i] == 0)) {
					matches++;
				} else {
					matches = 0;
					break;
				}

				if (j == pattern[0] - 1) {
					j = 1;
				}
			}
			if (matches == pattern2[0] - 2) {
				return jstart; // Search completed, patterns match
			}
			jstart += 2;
		}
		return -1;
	}

	/**
	 * Confirm whether the vertices having the given interior angles have the
	 * correct vertex pattern.
	 * 
	 * @param start     start index for the ang array
	 * @param ang       an array of interior angles
	 * @param vertexPat a boolean array indicating which angles are at actual
	 *                  vertices
	 * @param len       the length of the two arrays
	 * @return True if the pattern matches. Otherwise false.
	 */
	public boolean fitsVertexPat(byte start, double[] ang, boolean[] vertexPat, int len) {
		int i, j = start, k = 0, l;
		do {
			// Check the corresponding boolean in vertexPat
			if (vertexPat[k]) {
				// Compare ang[j] to all other angles at non-vertex vertices
				i = j + 1;
				if (i == len) {
					i = 0;
				}
				l = k + 1;
				if (l == len) {
					l = 0;
				}

				while (i != j) {
					if (!vertexPat[l] && ang[i] < ang[j]) {
						return false;
					}
					i++;
					if (i == len) {
						i = 0;
					}
					l++;
					if (l == len) {
						l = 0;
					}
				}
			}

			j++;
			if (j == len) {
				j = 0;
			}
			k++;
		} while (j != start);

		return true;
	}

	/**
	 * Fill the angles array with the angles at the opposite vertices.
	 * 
	 * @param ccwNeighbors the surrounding vertices in ccw order
	 * @param len          the length of
	 * @return an array of doubles
	 */
	public double[] surroundingAngles(Vertex[] ccwNeighbors, int len) {
		Quad q, qa, qb;
		Edge e, ep, en;
		Vertex n, np, nn, no;
		double[] angles = new double[len];
		for (int i = 0; i < len; i++) {
			n = ccwNeighbors[i];
			e = commonEdge(n);
			if (e == null) {
				if (i - 1 >= 0) {
					np = ccwNeighbors[i - 1];
				} else {
					np = ccwNeighbors[len - 1];
				}

				if (i + 1 < len) {
					nn = ccwNeighbors[i + 1];
				} else {
					nn = ccwNeighbors[0];
				}

				ep = commonEdge(np);
				en = commonEdge(nn);
				q = (Quad) ep.commonElement(en);

				no = q.oppositeVertex(this);
				angles[i] = q.ang[q.angleIndex(no)];
			} else {
				no = e.otherVertex(this);
				qa = (Quad) e.element1;
				qb = (Quad) e.element2;

				angles[i] = qa.ang[qa.angleIndex(no)] + qb.ang[qb.angleIndex(no)];
			}
		}
		return angles;
	}

	/**
	 * Compare the valence pattern of this boundary Vertex to the special pattern in
	 * pattern2. In pattern2, the following codes apply:<br>
	 * <ul>
	 * 14 means 4- (4 or less)
	 * <li>24 means 4+ (4 or more)
	 * <li>5 means 5 or more
	 * <li>0 means any number
	 * </ul>
	 * Note that calcMyValencePattern() must be called before calling this method.
	 * 
	 * @param pattern2 A valence pattern
	 * @param bpat     a boolean pattern indicating which vertices are located on
	 *                 the boundary
	 * @return If they match then return the true, otherwise return false.
	 */
	public boolean boundaryPatternMatch(byte[] pattern2, boolean[] bpat, Vertex[] ccwNeighbors) {
		if (pattern[0] != pattern2[0] || pattern[1] != pattern2[1] || bpat[0] != boundaryVertex()) {
			return false;
		}
		int i;

		for (i = 2; i < pattern[0]; i++) {
			if (pattern[i] == 2 && (pattern2[i] == 2 || pattern2[i] == 14 || pattern2[i] == 0)) {
				if (bpat[i - 1] && !ccwNeighbors[i - 2].boundaryVertex()) {
					return false;
				}
			} else if (pattern[i] == 3 && (pattern2[i] == 3 || pattern2[i] == 14 || pattern2[i] == 0)) {
				if (bpat[i - 1] && !ccwNeighbors[i - 2].boundaryVertex()) {
					return false;
				}
			} else if (pattern[i] == 4 && (pattern2[i] == 4 || pattern2[i] == 14 || pattern2[i] == 24 || pattern2[i] == 0)) {
				if (bpat[i - 1] && !ccwNeighbors[i - 2].boundaryVertex()) {
					return false;
				}
			} else if (pattern[i] >= 5 && (pattern2[i] == 5 || pattern2[i] == 24 || pattern2[i] == 0)) {
				if (bpat[i - 1] && !ccwNeighbors[i - 2].boundaryVertex()) {
					return false;
				}
			} else {
				return false;
			}
		}
		return true;
	}

	/**
	 * Compare the valence pattern of this internal Vertex to the special pattern in
	 * pattern2. The boundary pattern must also fit. In pattern2, the following
	 * codes apply:<br>
	 * <ul>
	 * 14 means 4- (4 or less)
	 * <li>24 means 4+ (4 or more)
	 * <li>5 means 5 or more
	 * <li>0 means any number
	 * </ul>
	 * Note that calcMyValencePattern() must be called before calling this method.
	 * 
	 * @param pattern2     A valence pattern
	 * @param bpat         a boolean pattern indicating which vertices are located
	 *                     on the boundary
	 * @param ccwNeighbors the neighbor vertices in ccw order
	 * @return If they match then return the true, otherwise return false.
	 */

	public int boundaryPatternMatchSpecial(byte[] pattern2, boolean[] bpat, Vertex[] ccwNeighbors) {

		if (pattern[0] != pattern2[0] || pattern[1] != pattern2[1] || bpat[0] != boundaryVertex()) {
			return -1;
		}
		int i, j, k;
		boolean match;

		for (k = 2; k < pattern[0]; k++) {

			j = k;
			match = true;

			for (i = 2; i < pattern[0]; i++) {
				if (pattern[j] == 2 && (pattern2[i] == 2 || pattern2[i] == 14 || pattern2[i] == 0)) {
					if (bpat[i - 1] && !ccwNeighbors[j - 2].boundaryVertex()) {
						match = false;
						break;
					}
				} else if (pattern[j] == 3 && (pattern2[i] == 3 || pattern2[i] == 14 || pattern2[i] == 0)) {
					if (bpat[i - 1] && !ccwNeighbors[j - 2].boundaryVertex()) {
						match = false;
						break;
					}
				} else if (pattern[j] == 4 && (pattern2[i] == 4 || pattern2[i] == 14 || pattern2[i] == 24 || pattern2[i] == 0)) {
					if (bpat[i - 1] && !ccwNeighbors[j - 2].boundaryVertex()) {
						match = false;
						break;
					}
				} else if (pattern[j] >= 5 && (pattern2[i] == 5 || pattern2[i] == 24 || pattern2[i] == 0)) {
					if (bpat[i - 1]) {
					} else {
					}

					if (ccwNeighbors[j - 2].boundaryVertex()) {
					} else {
					}

					if (bpat[i - 1] && !ccwNeighbors[j - 2].boundaryVertex()) {
						match = false;
						break;
					}
				} else {
					match = false;
					break;
				}

				j++;
				if (j == pattern[0]) {
					j = 2;
				}
			}
			if (match) {
				return k;
			}
		}
		return -1;
	}

	public Edge commonEdge(Vertex n) {
		Vertex other;
		for (Edge e : edgeList) {
			other = e.otherVertex(this);
			if (other == n) {
				return e;
			}
		}
		return null;
	}

	public boolean replaceWithStdMesh() {
		return true;
	}

	/**
	 * Give a string representation of the Vertex.
	 * 
	 * @return a string representation of the Vertex.
	 */
	public String descr() {
		return "(" + x + ", " + y + ")";
	}

	public String valDescr() {
		String s = "" + pattern[1] + "-";
		for (int i = 2; i < pattern[0]; i++) {
			s = s + pattern[i];
		}

		return s;
	}

	/** Output a string representation of the Vertex. */
	public void printMe() {
		System.out.println(descr());
	}
}
