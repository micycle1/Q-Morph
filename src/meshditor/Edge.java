package meshditor;

import java.util.ArrayList;
import java.util.List;

/**
 * This class holds information for edges, and has methods for handling issues
 * involving edges.
 */
public class Edge extends Constants {

	public Vertex leftVertex, rightVertex; // This Edge has these two vertices
	Element element1 = null, element2 = null; // Belongs to these Elements (Quads/Triangles)
	Edge leftFrontNeighbor, rightFrontNeighbor;
	int level;

	boolean frontEdge = false;
	boolean swappable = true;
	boolean selectable = true;
	// Edge leftSide= null, rightSide= null; // Side edges when building a quad
	boolean leftSide = false, rightSide = false; // Indicates if frontNeighbor is
													// to be used as side edge in quad
	double len; // length of this edge
	java.awt.Color color = java.awt.Color.green;

	static ArrayList<Edge>[] stateList = new ArrayList[3];

	public Edge(Vertex Vertex1, Vertex Vertex2) {
		if ((Vertex1.x < Vertex2.x) || (Vertex1.x == Vertex2.x && Vertex1.y > Vertex2.y)) {
			leftVertex = Vertex1;
			rightVertex = Vertex2;
		} else {
			leftVertex = Vertex2;
			rightVertex = Vertex1;
		}

		len = computeLength();
	}

	// Create a clone of Edge e with all the important fields
	private Edge(Edge e) {
		leftVertex = e.leftVertex;
		rightVertex = e.rightVertex;
		len = e.len;
		e.element1 = element1;
		e.element2 = element2;
	}

	@Override
	public boolean equals(Object o) {
		if (o instanceof Edge) {
			Edge e = (Edge) o;
			if (leftVertex.equals(e.leftVertex) && rightVertex.equals(e.rightVertex)) {
				return true;
			}
		}
		return false;
	}

	@Override
	public int hashCode() {
		return leftVertex.hashCode() ^ rightVertex.hashCode(); // direction-agnostic
//		long bits0 = Double.doubleToLongBits(leftVertex.x);
//		bits0 ^= Double.doubleToLongBits(leftVertex.y) * 31;
//		int hash0 = (((int) bits0) ^ ((int) (bits0 >> 32)));
//
//		long bits1 = Double.doubleToLongBits(rightVertex.x);
//		bits1 ^= Double.doubleToLongBits(rightVertex.y) * 31;
//		int hash1 = (((int) bits1) ^ ((int) (bits1 >> 32)));
//
//		return hash0 ^ hash1;
	}

	// Return a copy of the edge
	public Edge copy() {
		return new Edge(this);
	}

	public static void clearStateList() {
		stateList[0] = new ArrayList<>();
		stateList[1] = new ArrayList<>();
		stateList[2] = new ArrayList<>();
	}

	// Removes an Edge from the stateLists
	// Returns true if the Edge was successfully removed, else false.
	public boolean removeFromStateList() {
		int i;
		int state = getState();
		i = stateList[state].indexOf(this);
		if (i == -1) {
			return false;
		}
		stateList[state].remove(i);
		return true;
	}

	// Removes an Edge from the stateLists
	// Returns true if the Edge was successfully removed, else false.
	public boolean removeFromStateList(int state) {
		int i;
		i = stateList[state].indexOf(this);
		if (i == -1) {
			return false;
		}
		stateList[state].remove(i);
		return true;
	}

	public int getState() {
		int ret = 0;
		if (leftSide) {
			ret++;
		}
		if (rightSide) {
			ret++;
		}
		return ret;
	}

	public int getTrueState() {
		return getState();
	}

	public boolean alterLeftState(boolean newLeftState) {
		int state = getState();
		int i = stateList[state].indexOf(this);
		if (i == -1) {
			return false;
		}
		leftSide = newLeftState;
		int newState = getState();
		if (state != newState) {
			stateList[state].remove(i);
			stateList[newState].add(this);
		}
		return true;
	}

	public boolean alterRightState(boolean newRightState) {
		int state = getState();
		int i = stateList[state].indexOf(this);
		if (i == -1) {
			return false;
		}
		rightSide = newRightState;
		int newState = getState();
		if (state != newState) {
			stateList[state].remove(i);
			stateList[newState].add(this);
		}
		return true;
	}

	// Determine whether the frontNeighbor is an appropriate side Edge for a future
	// Quad with Edge e as base Edge. Return the frontNeighbor if so, else null.
	// elem== null if n.boundaryVertex()== true
	public Edge evalPotSideEdge(Edge frontNeighbor, Vertex n) {
		Element tri = getTriangleElement(), quad = getQuadElement();
		double ang;

		if (tri != null) {
			ang = sumAngle(tri, n, frontNeighbor);
		} else {
			ang = TWO_PI - sumAngle(quad, n, frontNeighbor);
		}

		if (ang < PIx3div4) { // if (ang< PIdiv2+EPSILON) // Could this be better?
			return frontNeighbor;
		} else {
			return null;
		}
	}

	// Determine the state bit at both Vertices and set the left and right side
	// Edges.
	// If a state bit is set, then the corresponding front neighbor Edge must get
	// Edge this as a side Edge at that Vertex. If it is not set, then it must get a
	// null
	// value instead.

	// this: the edge to be classified (gets a state value, and is added to a
	// statelist)
	public void classifyStateOfFrontEdge() {
		final Edge lfn = leftFrontNeighbor, rfn = rightFrontNeighbor;

		Edge l, r;

		// Alter states and side Edges on left side:
		l = evalPotSideEdge(lfn, leftVertex);
		if (l != null) {
			leftSide = true;
			if (leftVertex == lfn.leftVertex) {
				lfn.alterLeftState(true);
			} else {
				lfn.alterRightState(true);
			}
		} else {
			leftSide = false;
			if (leftVertex == lfn.leftVertex) {
				lfn.alterLeftState(false);
			} else {
				lfn.alterRightState(false);
			}
		}

		// Alter states and side Edges on right side:
		r = evalPotSideEdge(rfn, rightVertex);
		if (r != null) {
			rightSide = true;
			if (rightVertex == rfn.leftVertex) {
				rfn.alterLeftState(true);
			} else {
				rfn.alterRightState(true);
			}
		} else {
			rightSide = false;
			if (rightVertex == rfn.leftVertex) {
				rfn.alterLeftState(false);
			} else {
				rfn.alterRightState(false);
			}
		}

		// Add this to a stateList:
		stateList[getState()].add(this);
	}

	public boolean isLargeTransition(Edge e) {
		double ratio;
		double e1Len = length();
		double e2Len = e.length();

		if (e1Len > e2Len) {
			ratio = e1Len / e2Len;
		} else {
			ratio = e2Len / e1Len;
		}

		if (ratio > 2.5) {
			return true;
		} else {
			return false;
		}
	}

	// Select the next front to be processed. The selection criteria is:
	// Primary: the edge state
	// Secondary: the edge level
	// If the candidate edge is part of a large transition on the front where
	// longest - shortest length ratio > 2.5, and the candidate edge is not in
	// state 1-1, then the shorter edge is selected.
	public static Edge getNextFront(/* ArrayList frontList, */) {
		Edge current, selected = null;
		int selState, curState = 2, i;

		// Select a front preferrably in stateList[2]

		while (curState >= 0 && selected == null) {
			for (i = 0; i < stateList[curState].size(); i++) {
				current = stateList[curState].get(i);
				if (current.selectable) {
					selected = current;
					break;
				}
			}
			curState--;
		}
		if (selected == null) {
			Msg.warning("getNextFront(): no selectable fronts found in stateLists.");
			return null;
		}

		selState = selected.getState();

		for (i = 0; i < stateList[selState].size(); i++) {
			current = stateList[selState].get(i);

			if (current.selectable
					&& (current.level < selected.level || (current.level == selected.level && current.length() < selected.length()))) {
				selected = current;
			}
		}

		if (selState != 2) {
			if (selected.isLargeTransition(selected.leftFrontNeighbor)) {
				if (selected.length() > selected.leftFrontNeighbor.length() && selected.leftFrontNeighbor.selectable) {
					return selected.leftFrontNeighbor;
				}
			}

			if (selected.isLargeTransition(selected.rightFrontNeighbor)) {
				if (selected.length() > selected.rightFrontNeighbor.length() && selected.rightFrontNeighbor.selectable) {
					return selected.rightFrontNeighbor;
				}
			}
		}
		return selected;
	}

	public static void markAllSelectable() {
		for (int i = 0; i < 3; i++) {
			for (Edge e : stateList[i]) {
				e.selectable = true;
			}
		}
	}

	public static void printStateLists() {
		if (Msg.debugMode) {
			System.out.println("frontsInState 1-1:");
			for (Edge edge : stateList[2]) {
				System.out.println("" + edge.descr() + ", (" + edge.getState() + ")");
			}
			System.out.println("frontsInState 0-1 and 1-0:");
			for (Edge edge : stateList[1]) {
				System.out.println("" + edge.descr() + ", (" + edge.getState() + ")");
			}
			System.out.println("frontsInState 0-0:");
			for (Edge edge : stateList[0]) {
				System.out.println("" + edge.descr() + ", (" + edge.getState() + ")");
			}
		}
	}

	// If e.leftVertex is leftmore than this.leftVertex, return true, else false
	public boolean leftTo(Edge e) {
		if ((leftVertex.x < e.leftVertex.x) || (leftVertex.x == e.leftVertex.x && leftVertex.y < e.leftVertex.y)) {
			return true;
		} else {
			return false;
		}
	}

	public boolean isFrontEdge() {
		if ((element1 instanceof Triangle && !(element2 instanceof Triangle))
				|| (element2 instanceof Triangle && !(element1 instanceof Triangle))) {
			return true;
		} else {
			return false;
		}

		/*
		 * if ((element1 instanceof Quad && element2 instanceof Quad) || (element1
		 * instanceof Triangle && element2 instanceof Triangle)) return false; else
		 * return true;
		 */
	}

	public String descr() {
		return "(" + leftVertex.x + ", " + leftVertex.y + "), (" + rightVertex.x + ", " + rightVertex.y + ")";
	}

	public void printMe() {
		System.out.println(descr());
	}

	public double length() {
		return len;
	}

	// Replace this edge's Vertex n1 with the Vertex n2:
	public boolean replaceVertex(Vertex n1, Vertex n2) {
		if (leftVertex.equals(n1)) {
			leftVertex = n2;
		} else if (rightVertex.equals(n1)) {
			rightVertex = n2;
		} else {
			return false;
		}
		len = computeLength();
		return true;
	}

	// Seam two edges together. The method assumes that they already have one common
	// Vertex. For each edge in the edgeList of the otherVertex (otherE) of edge e:
	// Provided that otherThis or the otherVertex of otherE is not found in any edge
	// in
	// the edgeList of the otherVertex of this edge (otherThis), then the edge is
	// added
	// to this' edgeList.

	// Assumes that the quad area defined by the three distinct vertices of the two
	// edges,
	// and the Vertex in the top of the triangle adjacent the otherVertices (of the
	// common
	// Vertex of the two edges), is empty.

	public void seamWith(Edge e) {
		Vertex nK = commonVertex(e);
		Vertex nKp1 = otherVertex(nK), nKm1 = e.otherVertex(nK), other;
		boolean found = false;
		Edge eI, eJ;

		for (int i = 0; i < nKm1.edgeList.size(); i++) {
			eI = nKm1.edgeList.get(i);
			other = eI.otherVertex(nKm1);

			if (other != nKp1) {
				for (int j = 0; j < nKp1.edgeList.size(); j++) {
					eJ = nKp1.edgeList.get(j);

					if (other == eJ.otherVertex(nKp1)) {
						found = true;

						other.edgeList.remove(other.edgeList.indexOf(eI));

						if (eI.element1.firstVertex == nKm1) {
							eI.element1.firstVertex = nKp1;
						}
						eI.element1.replaceEdge(eI, eJ);
						eJ.connectToElement(eI.element1);
						break;
					}
				}
				if (!found) {
					if (eI.element1.firstVertex == nKm1) {
						eI.element1.firstVertex = nKp1;
					}
					if (eI.element2.firstVertex == nKm1) {
						eI.element2.firstVertex = nKp1;
					}

					eI.replaceVertex(nKm1, nKp1);
					nKp1.edgeList.add(eI);
				} else {
					found = false;
				}
			} else {
				// Remove the edge between eKp1 and eKm1 (from the edgeList of eKp1)
				nKp1.edgeList.remove(nKp1.edgeList.indexOf(eI));
			}
		}
	}

	// Return the midpoint (represented by a new Vertex) of this edge:
	public Vertex midPoint() {
		double xDiff = rightVertex.x - leftVertex.x;
		double yDiff = rightVertex.y - leftVertex.y;

		return new Vertex(leftVertex.x + xDiff * 0.5, leftVertex.y + yDiff * 0.5);
	}

	public double computeLength() {
		double xdiff = rightVertex.x - leftVertex.x;
		double ydiff = rightVertex.y - leftVertex.y;
		return Math.sqrt(xdiff * xdiff + ydiff * ydiff);
	}

	public static double length(double x1, double y1, double x2, double y2) {
		double xdiff = x2 - x1;
		double ydiff = y2 - y1;
		return Math.sqrt(xdiff * xdiff + ydiff * ydiff);
	}

	public static double length(Vertex Vertex1, Vertex Vertex2) {
		double xdiff = Vertex2.x - Vertex1.x;
		double ydiff = Vertex2.y - Vertex1.y;
		return Math.sqrt(xdiff * xdiff + ydiff * ydiff);
	}

	// Returns angle relative to the x-axis (which is directed from the origin (0,0)
	// to
	// the right, btw) at Vertex n (which is leftVertex or rightVertex).
	// Range: <-180, 180>
	public double angleAt(Vertex n) {
		// acos returns values in the range 0.0 through pi, avoiding neg numbers.
		// If this is CW to x-axis, then return pos acos, else return neg acos
		// double x= leftVertex.x-rightVertex.x;
		// double y= leftVertex.y-rightVertex.y;
		Vertex other = otherVertex(n);

		double x = n.x - other.x;
		double y = n.y - other.y;

		if (x == 0 && y > 0) {
			return -PIdiv2;
		} else if (x == 0 && y < 0) {
			return PIdiv2;
		} else {
			double hyp = Math.sqrt(x * x + y * y);

			if (x > 0) {
				if (y > 0) {
					return Math.PI + acos(x / hyp);
				} else {
					return Math.PI - acos(x / hyp);
				}
			} else if (y > 0) {
				return Math.PI + Math.PI - acos(-x / hyp);
			} else {
				return acos(-x / hyp);
			}

			// double cLen= Math.sqrt(x*x + y*y);
			/*
			 * double aLen= Math.sqrt(x*x + y*y); if (y> 0) return -acos((aLen*aLen + x*x
			 * -y*y)/(2*aLen*Math.abs(x))); else return acos((aLen*aLen + x*x
			 * -y*y)/(2*aLen*Math.abs(x)));
			 */
		}
	}

	/**
	 * Returns the angle from this Edge to eEdge by summing the angles of the
	 * elements adjacent to vertex n.
	 */
	public double sumAngle(Element sElem, Vertex n, Edge eEdge) {
		Element curElem = sElem;
		Edge curEdge = this;
		double ang = 0, iang = 0;
		double d;

		while (curEdge != eEdge && curElem != null) {
			d = curElem.angle(curEdge, n);
			ang += d;
			curEdge = curElem.neighborEdge(n, curEdge);
			curElem = curElem.neighbor(curEdge);
		}

		// If n is a boundaryVertex, the situation gets more complicated:
		if (curEdge != eEdge) {
			// Sum the "internal" angles:
			curEdge = this;
			curElem = sElem.neighbor(curEdge);

			while (curEdge != eEdge && curElem != null) {
				d = curElem.angle(curEdge, n);
				iang += d;
				curEdge = curElem.neighborEdge(n, curEdge);
				curElem = curElem.neighbor(curEdge);
			}
			ang = TWO_PI - iang;
		}
		return ang;
	}

	// Return the first front Edge adjacent Vertex n starting check at
	// this Edge in Element sElem.
	public Edge firstFrontEdgeAt(Element sElem, Vertex n) {
		Element curElem = sElem;
		Edge curEdge = this;

		while (!curEdge.frontEdge) {
			curEdge = curElem.neighborEdge(n, curEdge);
			curElem = curElem.neighbor(curEdge);
		}
		return curEdge;
	}

	/**
	 * Compute the internal angle between this Edge and Edge edge at Vertex n.
	 * 
	 * @param edge
	 * @param n
	 * @return a positive value.
	 */
	public double computePosAngle(Edge edge, Vertex n) {
		Vertex middle = commonVertex(edge);
		Vertex a = leftVertex == middle ? rightVertex : leftVertex;
		Vertex b = edge.leftVertex == middle ? edge.rightVertex : edge.leftVertex;
//		return angleBetween(a, middle, b)+1;
		return 0;

//		double a, b, c;
//		if (edge == this) {
//			Msg.warning("Edge.computePosAngle(..): The parameter Edge is the same as this Edge.");
//			return 2 * Math.PI;
//		}
//
//		if (leftVertex.equals(n)) {
//			if (n.equals(edge.leftVertex)) {
//				c = length(rightVertex, edge.rightVertex);
//			} else if (n.equals(edge.rightVertex)) {
//				c = length(rightVertex, edge.leftVertex);
//			} else {
//				Msg.error("Edge::computePosAngle(..): These edges are not connected.");
//				return 0;
//			}
//		} else if (rightVertex.equals(n)) {
//			if (n.equals(edge.leftVertex)) {
//				c = length(leftVertex, edge.rightVertex);
//			} else if (n.equals(edge.rightVertex)) {
//				c = length(leftVertex, edge.leftVertex);
//			} else {
//				Msg.error("Edge::computePosAngle(..): These edges are not connected.");
//				return 0;
//			}
//		} else {
//			Msg.error("Edge::computePosAngle(..): These edges are not connected.");
//			return 0;
//		}
//		a = len;
//		b = edge.len;
//
//		// acos returns a value in the range [0, PI],
//		// and input *MUST BE STRICTLY* in the range [-1, 1] !!!!!!!!
//		// ^^^^^^^^
//		double itemp = (a * a + b * b - c * c) / (2 * a * b);
//		if (itemp > 1.0) {
//			return 0;
//		} else if (itemp < -1.0) {
//			return Math.PI;
//		} else {
//			return acos(itemp);
//		}
	}

	/**
	 * Computes the CCW-directed angle between this edge and other edge.
	 * 
	 * @param edge
	 * @return a positive value in range [0, 2*PI]
	 */
	public double computeCCWAngle(Edge edge) {
		Vertex middle = commonVertex(edge);
		Vertex a = leftVertex == middle ? rightVertex : leftVertex;
		Vertex b = edge.leftVertex == middle ? edge.rightVertex : edge.leftVertex;

		return angleBetweenOriented(a, middle, b);
	}

	/**
	 * Returns the oriented smallest angle between two vectors. The computed angle
	 * will be in the range (-Pi, Pi]. A positive result corresponds to a
	 * counterclockwise (CCW) rotation from v1 to v2; a negative result corresponds
	 * to a clockwise (CW) rotation; a zero result corresponds to no rotation.
	 *
	 * @param tip1 the tip of v1
	 * @param tail the tail of each vector
	 * @param tip2 the tip of v2
	 * @return the angle between v1 and v2, relative to v1
	 */
	private static double angleBetweenOriented(Vertex tip1, Vertex tail, Vertex tip2) {
		double a1 = angle(tail, tip1);
		double a2 = angle(tail, tip2);
		double angDel = a2 - a1;

		// normalize, maintaining orientation
		if (angDel <= -Math.PI) {
			return angDel + TWO_PI;
		}
		if (angDel > Math.PI) {
			return angDel - TWO_PI;
		}
		return angDel;
	}

	/**
	 * Returns the angle of the vector from p0 to p1, relative to the positive
	 * X-axis. The angle is normalized to be in the range [ -Pi, Pi ].
	 *
	 * @param p0 the initial point of the vector
	 * @param p1 the terminal point of the vector
	 * @return the normalized angle (in radians) that p0-p1 makes with the positive
	 *         x-axis.
	 */
	private static double angle(Vertex p0, Vertex p1) {
		double dx = p1.x - p0.x;
		double dy = p1.y - p0.y;
		return atan2Quick(dy, dx);
	}

	private static double atan2Quick(final double y, final double x) {
		final double THREE_QRTR_PI = Math.PI * 0.75;
		final double QRTR_PI = Math.PI * 0.25;
	
		double r, angle;
		final double abs_y = Math.abs(y) + 1e-10f; // kludge to prevent 0/0 condition
	
		if (x < 0.0f) {
			r = (x + abs_y) / (abs_y - x); // (3)
			angle = THREE_QRTR_PI; // (4)
		} else {
			r = (x - abs_y) / (x + abs_y); // (1)
			angle = QRTR_PI; // (2)
		}
		angle += (0.1963f * r * r - 0.9817f) * r; // (2 | 4)
		if (y < 0.0f) {
			return (-angle); // negate if in quad III or IV
		} else {
			return (angle);
		}
	}

	private static double acos(double x) {
	//		return Math.acos(x);
			return acosQuick(x);
		}

	private static double acosQuick(double x) {
		return atan2Quick(Math.sqrt((1.0 + x) * (1.0 - x)), x);
	}

	/**
	 * Computes the unoriented smallest difference between two angles. The angles
	 * are assumed to be normalized to the range [-Pi, Pi]. The result will be in
	 * the range [0, Pi].
	 *
	 * @param ang1 the angle of one vector (in [-Pi, Pi] )
	 * @param ang2 the angle of the other vector (in range [-Pi, Pi] )
	 * @return the angle (in radians) between the two vectors (in range [0, Pi] )
	 */
	public static double diff(double ang1, double ang2) {
		double delAngle;

		if (ang1 < ang2) {
			delAngle = ang2 - ang1;
		} else {
			delAngle = ang1 - ang2;
		}

		if (delAngle > Math.PI) {
			delAngle = (2 * Math.PI) - delAngle;
		}

		return delAngle;
	}

	// Return a common Vertex for edges this and e
	public Vertex commonVertex(Edge e) {
		if (hasVertex(e.leftVertex)) {
			return e.leftVertex;
		} else if (hasVertex(e.rightVertex)) {
			return e.rightVertex;
		} else {
			return null;
		}
	}

	// Return a common element for edges this and e
	public Element commonElement(Edge e) {
		if (hasElement(e.element1)) {
			return e.element1;
		} else if (e.element2 != null && hasElement(e.element2)) {
			return e.element2;
		} else {
			return null;
		}
	}

	// Add this Edge to the vertices' edgeLists. Careful, there's no safety checks!
	public void connectVertices() {
		leftVertex.edgeList.add(this);
		rightVertex.edgeList.add(this);
	}

	// Remove this Edge from the vertices' edgeLists. Careful, there's no safety
	// checks!
	public void disconnectVertices() {
		leftVertex.edgeList.remove(leftVertex.edgeList.indexOf(this));
		rightVertex.edgeList.remove(rightVertex.edgeList.indexOf(this));
	}

	// Remove this Edge from the vertices' edgeLists. Safety checks...
	public void tryToDisconnectVertices() {
		int i;
		i = leftVertex.edgeList.indexOf(this);
		if (i != -1) {
			leftVertex.edgeList.remove(i);
		}
		i = rightVertex.edgeList.indexOf(this);
		if (i != -1) {
			rightVertex.edgeList.remove(i);
		}
	}

	public void connectToTriangle(Triangle triangle) {
		if (hasElement(triangle)) {
			return;
		}
		if (element1 == null) {
			element1 = triangle;
		} else if (element2 == null) {
			element2 = triangle;
		} else {
			Msg.error("Edge.connectToTriangle(..): An edge cannot be connected to more than two elements. edge= " + descr());
		}
	}

	public void connectToQuad(Quad q) {
		if (hasElement(q)) {
			return;
		}
		if (element1 == null) {
			element1 = q;
		} else if (element2 == null) {
			element2 = q;
		} else {
			Msg.error("Edge.connectToQuad(..):An edge cannot be connected to more than two elements.");
		}
	}

	public void connectToElement(Element elem) {
		if (hasElement(elem)) {
			return;
		}
		if (element1 == null) {
			element1 = elem;
		} else if (element2 == null) {
			element2 = elem;
		} else {
			Msg.error("Edge.connectToElement(..):An edge cannot be connected to more than two elements.");
		}
	}

	// element1 should never be null:
	public void disconnectFromElement(Element elem) {
		if (element1 == elem) {
			element1 = element2;
			element2 = null;
		} else if (element2 == elem) {
			element2 = null;
		} else {
			Msg.error("Edge " + descr() + " is not connected to element " + elem.descr() + ".");
		}
	}

	/**
	 * @param wrongVertex a Vertex that we don't want returned
	 * @return a Vertex opposite to this edge in an adjacent triangle
	 */
	public Vertex oppositeVertex(Vertex wrongVertex) {
		Vertex candidate;
		Edge otherEdge;

		// Pick one of the other edges in element1
		int ind = element1.indexOf(this);
		if (ind == 0 || ind == 1) {
			otherEdge = element1.edgeList[2];
		} else {
			otherEdge = element1.edgeList[1];
		}

		// This edge contains an opposite Vertex... get this Vertex....
		if (otherEdge.leftVertex != leftVertex && otherEdge.leftVertex != rightVertex) {
			candidate = otherEdge.leftVertex;
		} else {
			candidate = otherEdge.rightVertex;
		}

		// Damn, it's the wrong Vertex! Then we must go look in element2.
		if (candidate == wrongVertex) {

			// Pick one of the other edges in element2
			ind = element2.indexOf(this);
			if (ind == 0 || ind == 1) {
				otherEdge = element2.edgeList[2];
			} else if (ind == 2) {
				otherEdge = element2.edgeList[1];
			}

			// This edge contains an opposite Vertex
			// get this Vertex....
			if (otherEdge.leftVertex != leftVertex && otherEdge.leftVertex != rightVertex) {
				candidate = otherEdge.leftVertex;
			} else if (otherEdge.rightVertex != leftVertex && otherEdge.rightVertex != rightVertex) {
				candidate = otherEdge.rightVertex;
			}
		}
		return candidate;
	}

	// Construct an Edge that is a unit normal to this Edge.
	// (Remember to add the new Vertex to vertexList if you want to keep it.)
	// nB: one of the vertices on this edge (leftVertex or rightVertex)

	/*
	 * C b ____----x ___---- | b=1 ___---- | x----------------------x A c B
	 * 
	 * angle(B)= PI/2
	 * 
	 */
	// xdiff= xB - xA
	// ydiff= yB - yA
	// a= 1, c= sqrt(xdiff^2 + ydiff^2)
	// xC = xA + b* cos(alpha +ang.A)
	// = xA + xdiff - ydiff*a/c
	// = xB - ydiff*a/c
	// = xB - ydiff/c
	//
	// yC = yA + b* sin(alpha + ang.A)
	// = yA + ydiff + xdiff*a/c
	// = yB + xdiff*a/c
	// = yB + xdiff/c
	//
	public Edge unitNormalAt(Vertex n) {
		double xdiff = rightVertex.x - leftVertex.x;
		double ydiff = rightVertex.y - leftVertex.y;

		double c = Math.sqrt(xdiff * xdiff + ydiff * ydiff);

		double xn = n.x - ydiff / c;
		double yn = n.y + xdiff / c;

		Vertex newVertex = new Vertex(xn, yn);

		return new Edge(n, newVertex);
	}

	public Edge getSwappedEdge() {
		if (element2 == null) {
			Msg.warning("getSwappedEdge: Cannot swap a boundary edge.");
			return null;
		}
		if (element1 instanceof Quad || element2 instanceof Quad) {
			Msg.warning("getSwappedEdge: Edge must lie between two triangles.");
			return null;
		}

		Vertex n = this.oppositeVertex(null);
		Vertex m = this.oppositeVertex(n);

		return new Edge(n, m);
	}

	/**
	 * Swap diagonal between the edge's two triangles and update locally (To be used
	 * with getSwappedEdge())
	 */
	public void swapToAndSetElementsFor(Edge e) {
		if (element1 == null || element2 == null) {
			Msg.error("Edge.swapToAndSetElementsFor(..): both elements not set");
		}

		// extract the outer edges

		Edge e1 = element1.neighborEdge(leftVertex, this);
		Edge e2 = element1.neighborEdge(e1.otherVertex(leftVertex), e1);
		Edge e3 = element2.neighborEdge(rightVertex, this);
		Edge e4 = element2.neighborEdge(e3.otherVertex(rightVertex), e3);

		element2.disconnectEdges(); // important: element2 *first*, then element1
		element1.disconnectEdges();

		Triangle t1 = new Triangle(e, e2, e3);
		Triangle t2 = new Triangle(e, e4, e1);

		t1.connectEdges();
		t2.connectEdges();

		// Update edgeLists at this.leftVertex and this.rightVertex
		// and at e.leftVertex and e.rightVertex:
		disconnectVertices();
		e.connectVertices();
	}

	public MyVector getVector() {
		return new MyVector(leftVertex, rightVertex);
	}

	public MyVector getVector(Vertex origin) {
		if (origin.equals(leftVertex)) {
			return new MyVector(leftVertex, rightVertex);
		} else if (origin.equals(rightVertex)) {
			return new MyVector(rightVertex, leftVertex);
		} else {
			Msg.error("Edge::getVector(Vertex): Vertex not an endpoint in this edge.");
			return null;
		}
	}

	public boolean bordersToTriangle() {
		if (element1 instanceof Triangle) {
			return true;
		} else if (element2 != null && element2 instanceof Triangle) {
			return true;
		} else {
			return false;
		}
	}

	public boolean boundaryEdge() {
		if (element1 == null || element2 == null) {
			return true;
		} else {
			return false;
		}
	}

	public boolean boundaryOrTriangleEdge() {
		if (element1 == null || element2 == null || element1 instanceof Triangle || element2 instanceof Triangle) {
			return true;
		} else {
			return false;
		}
	}

	public boolean hasVertex(Vertex n) {
		if (leftVertex == n || rightVertex == n) {
			return true;
		} else {
			return false;
		}
	}

	public boolean hasElement(Element elem) {
		if (element1 == elem || element2 == elem) {
			return true;
		} else {
			return false;
		}
	}

	public boolean hasFalseFrontNeighbor() {
		if (leftFrontNeighbor == null || !leftFrontNeighbor.frontEdge || rightFrontNeighbor == null || !rightFrontNeighbor.frontEdge) {
			return true;
		} else {
			return false;
		}
	}

	public boolean hasFrontNeighbor(Edge e) {
		if (leftFrontNeighbor == e || rightFrontNeighbor == e) {
			return true;
		} else {
			return false;
		}
	}

	public Vertex otherVertex(Vertex n) {
		if (n.equals(leftVertex)) {
			return rightVertex;
		} else if (n.equals(rightVertex)) {
			return leftVertex;
		} else {
			Msg.error("Edge.otherVertex(Vertex): n is not on this edge");
			return null;
		}
	}

	/**
	 * Extend this edge at a given Vertex and to a given lenth.
	 * 
	 * @param length the new length of this edge
	 * @param nJ     the Vertex from which the edge is extended
	 * @return the other Vertex on the new edge
	 */
	public Vertex otherVertexGivenNewLength(double length, Vertex nJ) {
		// First find the angle between the existing edge and the x-axis:
		// Use a triangle containing one 90 degrees angle:
		MyVector v = new MyVector(nJ, this.otherVertex(nJ));
		v.setLengthAndAngle(length, v.angle());

		// Use this to create the new Vertex:
		return new Vertex(v.origin.x + v.x, v.origin.y + v.y);
	}

	// Prefers the left Vertex if they are at equal y positions
	public Vertex upperVertex() {
		if (leftVertex.y >= rightVertex.y) {
			return leftVertex;
		} else {
			return rightVertex;
		}
	}

	// Prefers the right Vertex if they are at equal y positions
	public Vertex lowerVertex() {
		if (rightVertex.y <= leftVertex.y) {
			return rightVertex;
		} else {
			return leftVertex;
		}
	}

	/**
	 * Return true if the 1-orbit around this.commonVertex(e) through quad startQ
	 * from edge this to edge e doesn't contain any triangle elements. If Vertex n
	 * lies on the boundary, and the orbit contains a hole, the orbit simply skips
	 * the hole and continues on the other side.
	 */
	public boolean noTrianglesInOrbit(Edge e, Quad startQ) {
		Edge curEdge = this;
		Element curElem = startQ;
		Vertex n = commonVertex(e);
		if (n == null) {
			return false;
		}
		if (curEdge.boundaryEdge()) {
			curEdge = n.anotherBoundaryEdge(curEdge);
			curElem = curEdge.element1;
		}
		do {
			if (curElem instanceof Triangle) {
				return false;
			}
			curEdge = curElem.neighborEdge(n, curEdge);
			if (curEdge.boundaryEdge()) {
				curEdge = n.anotherBoundaryEdge(curEdge);
				curElem = curEdge.element1;
			} else {
				curElem = curElem.neighbor(curEdge);
			}
		} while (curEdge != e);

		return true;
	}

	public Edge findLeftFrontNeighbor(List<Edge> frontList) {
		List<Edge> list = new ArrayList<>();
		Edge candidate = null;
		double candAng = Double.POSITIVE_INFINITY, curAng;
		Triangle t;

		for (int j = 0; j < leftVertex.edgeList.size(); j++) {
			Edge leftEdge = leftVertex.edgeList.get(j);
			if (leftEdge != this && leftEdge.isFrontEdge() /* leftEdge.frontEdge */) {
				list.add(leftEdge);
			}
		}
		if (list.size() == 1) {
			return list.get(0);
		} else if (!list.isEmpty()) {

			// Choose the front edge with the smallest angle
			t = getTriangleElement();

			for (Edge leftEdge : list) {
				curAng = sumAngle(t, leftVertex, leftEdge);
				if (curAng < candAng) {
					candAng = curAng;
					candidate = leftEdge;
				}

				// if (noTrianglesInOrbit(leftEdge, q))
				// return leftEdge;
			}
			return candidate;
		}
		Msg.warning("findLeftFrontNeighbor(..): Returning null");
		return null;
	}

	public Edge findRightFrontNeighbor(List<Edge> frontList) {
		List<Edge> list = new ArrayList<>();
		Edge candidate = null;
		double candAng = Double.POSITIVE_INFINITY, curAng;
		Triangle t;

		for (int j = 0; j < rightVertex.edgeList.size(); j++) {
			Edge rightEdge = rightVertex.edgeList.get(j);
			if (rightEdge != this && rightEdge.isFrontEdge()/* frontList.contains(rightEdge) */) {
				list.add(rightEdge);
			}
		}
		if (list.size() == 1) {
			return list.get(0);
		} else if (!list.isEmpty()) {
			t = getTriangleElement();
			for (Edge rightEdge : list) {
				curAng = sumAngle(t, rightVertex, rightEdge);

				if (curAng < candAng) {
					candAng = curAng;
					candidate = rightEdge;
				}

				// if (noTrianglesInOrbit(rightEdge, q))
				// return rightEdge;
			}
			return candidate;
		}
		Msg.warning("findRightFrontNeighbor(..): List.size== " + list.size() + ". Returning null");
		return null;
	}

	/** Set the appropriate front neighbor to edge e. */
	public void setFrontNeighbor(Edge e) {
		if (e.hasVertex(leftVertex)) {
			leftFrontNeighbor = e;
		} else if (e.hasVertex(rightVertex)) {
			rightFrontNeighbor = e;
		} else {
			Msg.warning("Edge.setFrontNeighbor(..): Could not set.");
		}
	}

	/** Returns true if the frontEdgeNeighbors are changed. */
	public boolean setFrontNeighbors(List<Edge> frontList) {
		Edge lFront = findLeftFrontNeighbor(frontList);
		Edge rFront = findRightFrontNeighbor(frontList);
		boolean res = false;
		if (lFront != leftFrontNeighbor || rFront != rightFrontNeighbor) {
			res = true;
		}
		leftFrontNeighbor = lFront;
		rightFrontNeighbor = rFront;

		if (lFront != null && !lFront.hasFrontNeighbor(this)) {
			// res= true;
			lFront.setFrontNeighbor(this);
		}
		if (rFront != null && !rFront.hasFrontNeighbor(this)) {
			// res= true;
			rFront.setFrontNeighbor(this);
		}
		return res;
	}

	public void promoteToFront(int level, List<Edge> frontList) {
		if (!frontEdge) {
			frontList.add(this);
			this.level = level;
			frontEdge = true;
		}
	}

	public boolean removeFromFront(List<Edge> frontList) {
		int i = frontList.indexOf(this);
		frontEdge = false;
		if (i != -1) {
			frontList.remove(i);
			return true;
		} else {
			return false;
		}
	}

	/**
	 * Halve this Edge by introducing a new Vertex at the midpoint, and create two
	 * Edges from this midpoint to the each of the two opposite Vertices of Edge
	 * this: one in element1 and one in element2. Also create two new Edges from
	 * Vertex mid to the two Vertices of Edge this. Create four new Triangles.
	 * Update everything (also remove this Edge from edgeList and disconnect the
	 * vertices).
	 * 
	 * @return the new Edge incident with Vertex ben.
	 */
	public Edge splitTrianglesAt(Vertex nN, Vertex ben, List<Triangle> triangleList, List<Edge> edgeList, List<Vertex> vertexList) {
		Edge eK1 = new Edge(leftVertex, nN);
		Edge eK2 = new Edge(rightVertex, nN);

		Triangle tri1 = (Triangle) element1;
		Triangle tri2 = (Triangle) element2;

		Vertex n1 = tri1.oppositeOfEdge(this);
		Vertex n2 = tri2.oppositeOfEdge(this);
		Edge diagonal1 = new Edge(nN, n1);
		Edge diagonal2 = new Edge(nN, n2);

		Edge e12 = tri1.neighborEdge(leftVertex, this);
		Edge e13 = tri1.neighborEdge(rightVertex, this);
		Edge e22 = tri2.neighborEdge(leftVertex, this);
		Edge e23 = tri2.neighborEdge(rightVertex, this);

		Triangle t11 = new Triangle(diagonal1, e12, eK1);
		Triangle t12 = new Triangle(diagonal1, e13, eK2);
		Triangle t21 = new Triangle(diagonal2, e22, eK1);
		Triangle t22 = new Triangle(diagonal2, e23, eK2);

		// Update the vertices' edgeLists
		disconnectVertices();
		eK1.connectVertices();
		eK2.connectVertices();
		diagonal1.connectVertices();
		diagonal2.connectVertices();

		// Disconnect old Triangles
		tri1.disconnectEdges();
		tri2.disconnectEdges();

		// Connect Edges to new Triangles
		t11.connectEdges();
		t12.connectEdges();
		t21.connectEdges();
		t22.connectEdges();

		// Update "global" lists
		edgeList.remove(edgeList.indexOf(this));
		edgeList.add(eK1);
		edgeList.add(eK2);
		edgeList.add(diagonal1);
		edgeList.add(diagonal2);

		triangleList.remove(triangleList.indexOf(tri1));
		triangleList.remove(triangleList.indexOf(tri2));
		triangleList.add(t11);
		triangleList.add(t12);
		triangleList.add(t21);
		triangleList.add(t22);

		if (eK1.hasVertex(ben)) {
			return eK1;
		} else if (eK2.hasVertex(ben)) {
			return eK2;
		} else {
			Msg.error("");
			return null;
		}
	}

	/**
	 * Make new triangles by introducing new Edges at this' midpoint.
	 * 
	 * @return the "lower" (the one incident with the baseEdge) of the two edges
	 *         created from splitting this edge.
	 */
	public Edge splitTrianglesAtMyMidPoint(List<Triangle> triangleList, List<Edge> edgeList, List<Vertex> vertexList, Edge baseEdge) {
		Edge lowerEdge;
		Vertex ben = baseEdge.commonVertex(this);
		Vertex mid = this.midPoint();
		vertexList.add(mid);
		mid.color = java.awt.Color.blue;

		lowerEdge = splitTrianglesAt(mid, ben, triangleList, edgeList, vertexList);

		return lowerEdge;
	}

	/**
	 * Find the next edge adjacent a quad element, starting at this edge which is
	 * part of a given element and which is adjacent a given Vertex. Note that the
	 * method stops if the boundary is encountered.
	 * 
	 * @param n         the Vertex
	 * @param startElem
	 * @return the first edge of a quad found when parsing around Vertex n, starting
	 *         at edge e in element startElem and moving in the direction from e to
	 *         e's neighbor edge at n in startElem. If startElem happens to be a
	 *         quad, the method won't consider that particular quad. If
	 *         unsuccessful, the method returns null.
	 */
	public Edge nextQuadEdgeAt(Vertex n, Element startElem) {
		Element elem;
		Edge e;

		e = startElem.neighborEdge(n, this);
		elem = startElem.neighbor(e);

		while (elem != null && !(elem instanceof Quad) && elem != startElem) {
			e = elem.neighborEdge(n, e);
			elem = elem.neighbor(e);
		}
		if (elem != null && elem instanceof Quad && elem != startElem) {
			return e;
		} else {
			return null;
		}
	}

	// Returns a neighboring element that is a quad. When applied to inner front
	// edges, there are, of course, only one possible quad to return.
	public Quad getQuadElement() {
		if (element1 instanceof Quad) {
			return (Quad) element1;
		} else if (element2 instanceof Quad) {
			return (Quad) element2;
		} else {
			return null;
		}
	}

	// Returns a neighboring element that is a triangle. The method should work as
	// long as it is applied to a front edge.
	public Triangle getTriangleElement() {
		if (element1 instanceof Triangle) {
			return (Triangle) element1;
		} else if (element2 instanceof Triangle) {
			return (Triangle) element2;
		} else {
			return null;
		}
	}

	/** Return the neighboring quad that is also a neighbor of edge e. */
	public Quad getQuadWithEdge(Edge e) {
		if (element1 instanceof Quad && element1.hasEdge(e)) {
			return (Quad) element1;
		} else if (element2 instanceof Quad && element2.hasEdge(e)) {
			return (Quad) element2;
		} else {
			return null;
		}
	}

	/** Return the front neighbor edge at Vertex n. */
	public Edge frontNeighborAt(Vertex n) {
		if (leftFrontNeighbor != null && commonVertex(leftFrontNeighbor) == n) {
			return leftFrontNeighbor;
		} else if (rightFrontNeighbor != null && commonVertex(rightFrontNeighbor) == n) {
			return rightFrontNeighbor;
		} else {
			return null;
		}
	}

	/** @return the front neighbor next to this (not the prev edge). */
	public Edge nextFrontNeighbor(Edge prev) {
		if (leftFrontNeighbor != prev) {
			return leftFrontNeighbor;
		} else if (rightFrontNeighbor != prev) {
			return rightFrontNeighbor;
		} else {
			Msg.error("Edge.nextFrontNeighbor(Edge): Cannot find a suitable next edge.");
			return null;
		}
	}

	// Return a neighbor edge at Vertex n, that is front edge according to the
	// definition,
	// and that is part of the same loop as this edge.
	// Assumes that this is a true front edge.
	public Edge trueFrontNeighborAt(Vertex n) {
		Element curElem = getTriangleElement();
		Edge curEdge = this;

		if (!hasVertex(n)) {
			Msg.error("trueFrontNeighborAt(..): this Edge hasn't got Vertex " + n.descr());
		}

		do {
			curEdge = curElem.neighborEdge(n, curEdge);
			curElem = curElem.neighbor(curEdge);
		} while (!curEdge.isFrontEdge());

		return curEdge;
	}
}
