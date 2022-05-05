package meshditor;

/**
 * A class holding information for triangles, and with methods for handling issues regarding triangles.
 */
public class Triangle extends Element {
	
	public Triangle(double... vertices) {
		Vertex v1 = new Vertex(vertices[0], vertices[1]);
		Vertex v2 = new Vertex(vertices[2], vertices[3]);
		Vertex v3 = new Vertex(vertices[4], vertices[5]);
		Edge e1 = new Edge(v1, v2);
		Edge e2 = new Edge(v2, v3);
		Edge e3 = new Edge(v3, v1);

		edgeList = new Edge[3];
		edgeList[0] = e1;
		edgeList[1] = e2;
		edgeList[2] = e3;
		connectEdges();
		updateAngles();
	}

	public Triangle(Edge edge1, Edge edge2, Edge edge3) {
		edgeList = new Edge[3];
	
		if (edge1 == null || edge2 == null || edge3 == null) {
			Msg.error("Triangle: cannot create Triangle with null Edge.");
		}
	
		edgeList[0] = edge1;
		edgeList[1] = edge2;
		edgeList[2] = edge3;
	
		// Make a pointer to the base Vertex that is the origin of vector(edgeList[a])
		// so that the cross product vector(edgeList[a]) x vector(edgeList[b]) >= 0
		// where a,b in {1,2}
		firstVertex = edgeList[0].leftVertex;
		if (inverted()) {
			firstVertex = edgeList[0].rightVertex;
		}
	
		Edge temp;
		if (firstVertex == edgeList[0].commonVertex(edgeList[2])) {
			temp = edgeList[1];
			edgeList[1] = edgeList[2];
			edgeList[2] = temp;
		}
	
		ang = new double[3];
		updateAngles();
	}

	public Triangle(Edge edge1, Edge edge2, Edge edge3, double len1, double len2, double len3, double ang1, double ang2, double ang3,
			boolean lengthsOpt, boolean anglesOpt) {
		edgeList = new Edge[3];

		edgeList[0] = edge1;
		edgeList[1] = edge2;
		edgeList[2] = edge3;

		edgeList[0].len = len1;
		edgeList[1].len = len2;
		edgeList[2].len = len3;

		// Make a pointer to the base Vertex that is the origin of vector(edgeList[a])
		// so that the cross product vector(edgeList[a]) x vector(edgeList[b]) >= 0
		// where a,b in {1,2}
		firstVertex = edgeList[0].leftVertex;
		if (inverted()) {
			firstVertex = edgeList[0].rightVertex;
		}

		Edge temp;
		if (firstVertex == edgeList[0].commonVertex(edgeList[2])) {
			temp = edgeList[1];
			edgeList[1] = edgeList[2];
			edgeList[2] = temp;
		}

		ang = new double[3];
		ang[0] = ang1;
		ang[1] = ang2;
		ang[2] = ang3;

		if (!lengthsOpt) {
			updateLengths();
		}
		if (!anglesOpt) {
			updateAngles();
		}
	}

	// Makes a copy of the given triangle
	private Triangle(Triangle t) {
		edgeList = new Edge[3];

		edgeList[0] = new Edge(t.edgeList[0].leftVertex, t.edgeList[0].rightVertex);
		// t.edgeList[0].copy();
		edgeList[1] = new Edge(t.edgeList[1].leftVertex, t.edgeList[1].rightVertex);
		// t.edgeList[1].copy();
		edgeList[2] = new Edge(t.edgeList[2].leftVertex, t.edgeList[2].rightVertex);
		// t.edgeList[2].copy();

		ang = new double[3];
		ang[0] = t.ang[0];
		ang[1] = t.ang[1];
		ang[2] = t.ang[2];

		firstVertex = t.firstVertex;
	}

	// Create a simple triangle for testing purposes only
	// (constrainedLaplacianSmooth()
	// and optBasedSmooth(..))
	// The triangle created is a copy of this, and each of the edges are also
	// copies.
	// One of the vertices in the triangle has been replaced.
	// The Edge.repleceVertex(..) method updates the edge length.
	// *Not tested!!!*

	@Override
	public Element elementWithExchangedVertices(Vertex original, Vertex replacement) {
		Vertex Vertex1 = edgeList[0].leftVertex;
		Vertex Vertex2 = edgeList[0].rightVertex;
		Vertex Vertex3 = edgeList[1].rightVertex;
		if (Vertex3 == Vertex1 || Vertex3 == Vertex2) {
			Vertex3 = edgeList[1].leftVertex;
		}
		Vertex common1;

		// Make a copy of the original triangle...
		Triangle t = new Triangle(this);
		Edge edge1 = t.edgeList[0], edge2 = t.edgeList[1], edge3 = t.edgeList[2];

		// ... and then replace the Vertex
		if (Vertex1 == original) {
			common1 = edge1.commonVertex(edge2);
			if (original == common1) {
				edge2.replaceVertex(original, replacement);
			} else {
				edge3.replaceVertex(original, replacement);
			}
			edge1.replaceVertex(original, replacement);
		} else if (Vertex2 == original) {
			common1 = edge1.commonVertex(edge2);
			if (original == common1) {
				edge2.replaceVertex(original, replacement);
			} else {
				edge3.replaceVertex(original, replacement);
			}
			edge1.replaceVertex(original, replacement);
		} else if (Vertex3 == original) {
			common1 = edge2.commonVertex(edge1);
			if (original == common1) {
				edge1.replaceVertex(original, replacement);
			} else {
				edge3.replaceVertex(original, replacement);
			}
			edge2.replaceVertex(original, replacement);
		} else {
			return null;
		}

		if (original == firstVertex) {
			t.firstVertex = replacement;
		} else {
			t.firstVertex = firstVertex;
		}
		return t;
	}

	// Return true if the quad becomes inverted when Vertex n is relocated to pos.
	// (x,y).
	// Else return false.
	@Override
	public boolean invertedWhenVertexRelocated(Vertex n1, Vertex n2) {
		Msg.debug("Entering Triangle.invertedWhenVertexRelocated(..)");
		Vertex a, b, c;
		MyVector ac, bc;
		a = firstVertex;
		b = edgeList[0].otherVertex(a);
		c = edgeList[1].rightVertex;
		if (c == a || c == b) {
			c = edgeList[1].leftVertex;
		}

		if (a == n1) {
			a = n2;
		} else if (b == n1) {
			b = n2;
		} else if (c == n1) {
			c = n2;
		}

		ac = new MyVector(a, c);
		bc = new MyVector(b, c);

		Msg.debug("Leaving Triangle.invertedWhenVertexRelocated(..)");
		if (ac.cross(bc) <= 0) {
			return true;
		} else {
			return false;
		}
	}

	@Override
	public boolean equals(Object o) {
		if (o instanceof Triangle) {
			Triangle t = (Triangle) o;
			Vertex Vertex1, Vertex2, Vertex3;
			Vertex tVertex1, tVertex2, tVertex3;
			Vertex1 = edgeList[0].leftVertex;
			Vertex2 = edgeList[0].rightVertex;
			Vertex3 = edgeList[1].rightVertex;
			if (Vertex3 == Vertex1 || Vertex3 == Vertex2) {
				Vertex3 = edgeList[1].leftVertex;
			}

			tVertex1 = t.edgeList[0].leftVertex;
			tVertex2 = t.edgeList[0].rightVertex;
			tVertex3 = t.edgeList[1].rightVertex;
			if (tVertex3 == tVertex1 || tVertex3 == tVertex2) {
				tVertex3 = t.edgeList[1].leftVertex;
			}

			if (Vertex1 == tVertex1 && Vertex2 == tVertex2 && Vertex3 == tVertex3) {
				return true;
			} else {
				return false;
			}
		} else {
			return false;
		}
	}

	@Override
	public double angle(Edge e, Vertex n) {
		return ang[angleIndex(e, neighborEdge(n, e))];
	}

	@Override
	public double angle(Edge e1, Edge e2) {
		return ang[angleIndex(e1, e2)];
	}

	public void updateAttributes() {
		updateLengths();
		updateAngles();
	}

	public void updateLengths() {
		edgeList[0].len = edgeList[0].computeLength();
		edgeList[1].len = edgeList[1].computeLength();
		edgeList[2].len = edgeList[2].computeLength();
	}

	@Override
	public void updateAngle(Vertex n) {
		int j = angleIndex(n), i;

		i = angleIndex(edgeList[0], edgeList[1]);
		if (j == i) {
			ang[i] = edgeList[0].computeCCWAngle(edgeList[1]);
		} else {
			i = angleIndex(edgeList[1], edgeList[2]);
			if (j == i) {
				ang[i] = edgeList[1].computeCCWAngle(edgeList[2]);
			} else {
				i = angleIndex(edgeList[0], edgeList[2]);
				if (j == i) {
					ang[i] = edgeList[2].computeCCWAngle(edgeList[0]);
				} else {
					Msg.error("Quad.updateAngle(Vertex): Vertex not found.");
				}
			}
		}
	}

	@Override
	public void updateAngles() {
		int i;
		i = angleIndex(edgeList[0], edgeList[1]);
		ang[i] = edgeList[0].computeCCWAngle(edgeList[1]);

		i = angleIndex(edgeList[1], edgeList[2]);
		ang[i] = edgeList[1].computeCCWAngle(edgeList[2]);

		i = angleIndex(edgeList[0], edgeList[2]);
		ang[i] = edgeList[2].computeCCWAngle(edgeList[0]);
	}

	@Override
	public boolean hasEdge(Edge e) {
		if (edgeList[0] == e || edgeList[1] == e || edgeList[2] == e) {
			return true;
		} else {
			return false;
		}
	}

	@Override
	public boolean hasVertex(Vertex n) {
		if (edgeList[0].hasVertex(n) || edgeList[1].hasVertex(n) || edgeList[2].hasVertex(n)) {
			return true;
		} else {
			return false;
		}
	}

	// Return an Edge of this Triangle that is not Edge e.
	public Edge otherEdge(Edge e) {
		if (edgeList[0] != e) {
			return edgeList[0];
		} else if (edgeList[1] != e) {
			return edgeList[1];
		}
		if (edgeList[2] != e) {
			return edgeList[2];
		} else {
			Msg.error("Cannot find an other edge in triangle " + descr() + ".");
			return null;
		}
	}

	// Return the Edge of this Triangle that is not one of Edges e1 or e2.
	public Edge otherEdge(Edge e1, Edge e2) {
		if ((edgeList[0] == e1 && edgeList[1] == e2) || (edgeList[1] == e1 && edgeList[0] == e2)) {
			return edgeList[2];
		} else if ((edgeList[1] == e1 && edgeList[2] == e2) || (edgeList[2] == e1 && edgeList[1] == e2)) {
			return edgeList[0];
		} else if ((edgeList[0] == e1 && edgeList[2] == e2) || (edgeList[2] == e1 && edgeList[0] == e2)) {
			return edgeList[1];
		} else {
			Msg.error("Cannot find the other edge in triangle " + descr() + ".");
			return null;
		}
	}

	// Returns neighbor element sharing edge e
	@Override
	public Element neighbor(Edge e) {
		if (e.element1 == this) {
			return e.element2;
		} else if (e.element2 == this) {
			return e.element1;
		} else {
			Msg.warning("Triangle.neighbor(Edge): returning null");
			return null;
		}
	}

	// Returns the other edge belonging to this triangle that also share Vertex n
	// Returns null if not found
	@Override
	public Edge neighborEdge(Vertex n, Edge e) {
		if (edgeList[0] != e && edgeList[0].hasVertex(n)) {
			return edgeList[0];
		} else if (edgeList[1] != e && edgeList[1].hasVertex(n)) {
			return edgeList[1];
		} else if (edgeList[2] != e && edgeList[2].hasVertex(n)) {
			return edgeList[2];
		} else {
			Msg.error("Triangle.neighborEdge(Vertex, Edge): Neighbor not found.");
			return null; // e
		}
	}

	@Override
	public int indexOf(Edge e) {
		if (edgeList[0] == e) {
			return 0;
		} else if (edgeList[1] == e) {
			return 1;
		} else if (edgeList[2] == e) {
			return 2;
		} else {
			Msg.warning("Triangle.indexOf(Edge): Edge not part of triangle.");
			return -1;
		}
	}

	public int angleIndex(int e1Index, int e2Index) {
		// angle betw. edges 0 && 1
		if ((e1Index == 0 && e2Index == 1) || (e1Index == 1 && e2Index == 0)) {
			return 2;
		} else if ((e1Index == 1 && e2Index == 2) || (e1Index == 2 && e2Index == 1)) {
			return 0;
		} else if ((e1Index == 0 && e2Index == 2) || (e1Index == 2 && e2Index == 0)) {
			return 1;
		} else {
			return -1;
		}
	}

	@Override
	public int angleIndex(Vertex n) {
		// angle betw. edges 0 && 1
		if (edgeList[0].commonVertex(edgeList[1]) == n) {
			return 2;
		} else if (edgeList[1].commonVertex(edgeList[2]) == n) {
			return 0;
		} else if (edgeList[2].commonVertex(edgeList[0]) == n) {
			return 1;
		} else {
			Msg.error("Quad.angleIndex(Vertex): Vertex not found");
			return -1;
		}
	}

	@Override
	public int angleIndex(Edge e1, Edge e2) {
		return angleIndex(indexOf(e1), indexOf(e2));
	}

	public Edge oppositeOfVertex(Vertex n) {
		if (edgeList[0].hasVertex(n) && edgeList[1].hasVertex(n)) {
			return edgeList[2];
		} else if (edgeList[0].hasVertex(n) && edgeList[2].hasVertex(n)) {
			return edgeList[1];
		} else if (edgeList[1].hasVertex(n) && edgeList[2].hasVertex(n)) {
			return edgeList[0];
		} else {
			Msg.error("Cannot find opposide edge of Vertex.");
			return null;
		}
	}

	public Vertex oppositeOfEdge(Edge e) {
		if (edgeList[0] == e) {
			if (!e.hasVertex(edgeList[1].leftVertex)) {
				return edgeList[1].leftVertex;
			} else {
				return edgeList[1].rightVertex;
			}
		} else if (edgeList[1] == e) {
			if (!e.hasVertex(edgeList[2].leftVertex)) {
				return edgeList[2].leftVertex;
			} else {
				return edgeList[2].rightVertex;
			}
		} else if (edgeList[2] == e) {
			if (!e.hasVertex(edgeList[0].leftVertex)) {
				return edgeList[0].leftVertex;
			} else {
				return edgeList[0].rightVertex;
			}
		} else {
			Msg.error("oppositeOfEdge: Cannot find Vertex opposite of the supplied edge.");
			return null;
		}
	}

	// Make the Element pointers in each of the Edges in this Triangle's edgeList
	// point to this Triangle.
	@Override
	public void connectEdges() {
		edgeList[0].connectToTriangle(this);
		edgeList[1].connectToTriangle(this);
		edgeList[2].connectToTriangle(this);
	}

	// Release the element pointer of each Edge in edgeList that pointed to this
	// Element.
	@Override
	public void disconnectEdges() {
		edgeList[0].disconnectFromElement(this);
		edgeList[1].disconnectFromElement(this);
		edgeList[2].disconnectFromElement(this);
	}

	// Return the next ccw edge in this triangle.
	// Get the pos angles between e1 and each of the other two edges. Transform
	// these
	// angles into true angles. The largest angle indicates that this edge is ccw to
	// e1.
	public Edge nextCCWEdge(Edge e1) {
		Vertex e1commone2, e1commone3;
		Edge e2, e3;
		MyVector v1Forv2, v1Forv3, v2, v3;
		if (e1 == edgeList[0]) {
			e2 = edgeList[1];
			e3 = edgeList[2];
		} else if (e1 == edgeList[1]) {
			e2 = edgeList[0];
			e3 = edgeList[2];
		} else if (e1 == edgeList[2]) {
			e2 = edgeList[0];
			e3 = edgeList[1];
		} else {
			Msg.error("Edge " + e1.descr() + " is not part of triangle");
			return null;
		}

		e1commone2 = e1.commonVertex(e2);
		e1commone3 = e1.commonVertex(e3);

		v2 = new MyVector(e1commone2, e2.otherVertex(e1commone2));
		v3 = new MyVector(e1commone3, e3.otherVertex(e1commone3));
		v1Forv2 = new MyVector(e1commone2, e1.otherVertex(e1commone2));
		v1Forv3 = new MyVector(e1commone3, e1.otherVertex(e1commone3));

		// Positive angles between e1 and each of the other two edges:
		double ang1, ang2;
		ang1 = ang[angleIndex(e1, e2)];
		ang2 = ang[angleIndex(e1, e3)];

		if (v2.isCWto(v1Forv2)) {
			ang1 += Math.PI;
		}
		if (v3.isCWto(v1Forv3)) {
			ang2 += Math.PI;
		}

		if (ang2 > ang1) {
			return e3;
		} else {
			return e2;
		}
	}

	public Edge nextCWEdge(Edge e1) {
		Vertex e1commone2, e1commone3;
		Edge e2, e3;
		MyVector v1Forv2, v1Forv3, v2, v3;
		if (e1 == edgeList[0]) {
			e2 = edgeList[1];
			e3 = edgeList[2];
		} else if (e1 == edgeList[1]) {
			e2 = edgeList[0];
			e3 = edgeList[2];
		} else if (e1 == edgeList[2]) {
			e2 = edgeList[0];
			e3 = edgeList[1];
		} else {
			Msg.error("Edge " + e1.descr() + " is not part of triangle");
			return null;
		}

		e1commone2 = e1.commonVertex(e2);
		e1commone3 = e1.commonVertex(e3);

		v2 = new MyVector(e1commone2, e2.otherVertex(e1commone2));
		v3 = new MyVector(e1commone3, e3.otherVertex(e1commone3));
		v1Forv2 = new MyVector(e1commone2, e1.otherVertex(e1commone2));
		v1Forv3 = new MyVector(e1commone3, e1.otherVertex(e1commone3));

		// Positive angles between e1 and each of the other two edges:
		double ang1, ang2;
		ang1 = ang[angleIndex(e1, e2)];
		ang2 = ang[angleIndex(e1, e3)];

		// Transform into true angles:
		e1commone2 = e1.commonVertex(e2);
		e1commone3 = e1.commonVertex(e3);

		v2 = new MyVector(e1commone2, e2.otherVertex(e1commone2));
		v3 = new MyVector(e1commone3, e3.otherVertex(e1commone3));
		v1Forv2 = new MyVector(e1commone2, e1.otherVertex(e1commone2));
		v1Forv3 = new MyVector(e1commone3, e1.otherVertex(e1commone3));

		if (v2.isCWto(v1Forv2)) {
			ang1 += Math.PI;
		}
		if (v3.isCWto(v1Forv3)) {
			ang2 += Math.PI;
		}

		if (ang2 > ang1) {
			return e2;
		} else {
			return e3;
		}
	}

	// We simply check that the vertices of the element are not collinear.
	@Override
	public boolean areaLargerThan0() {
		Vertex na = edgeList[0].leftVertex;
		Vertex nb = edgeList[0].rightVertex;
		Vertex nc = oppositeOfEdge(edgeList[0]);

		if (cross(na, nc, nb, nc) != 0) {
			return true;
		} else {
			return false;
		}

		/*
		 * Vertex n= oppositeOfEdge(edgeList[0]); if (n.onLine(edgeList[0])) return false;
		 * else return true;
		 */
	}

	// Check if the old pos and the new pos of the Vertex are on different sides
	// of the vertices opposite edge.
	public boolean inverted(Vertex oldN, Vertex newN) {
		Edge e = oppositeOfVertex(newN);
		// Check with edge e:
		double oldN_e_det = (e.leftVertex.x - oldN.x) * (e.rightVertex.y - oldN.y) - (e.leftVertex.y - oldN.y) * (e.rightVertex.x - oldN.x);
		double newN_e_det = (e.leftVertex.x - newN.x) * (e.rightVertex.y - newN.y) - (e.leftVertex.y - newN.y) * (e.rightVertex.x - newN.x);

		// If different sign, or 0, they are inverted:
		if (oldN_e_det >= 0) {
			if (newN_e_det <= 0 || oldN_e_det == 0) {
				Msg.debug("Triangle.inverted(..): Triangle " + descr() + " is inverted");
				return true;
			}
		} else if (newN_e_det >= 0) {
			Msg.debug("Triangle.inverted(..): Triangle " + descr() + " is inverted");
			return true;
		}
		return false;
	}

	// Return true if the triangle has become inverted
	@Override
	public boolean inverted() {
		Vertex a, b, c;
		MyVector ac, bc;
		a = firstVertex;
		b = edgeList[0].otherVertex(a);
		c = edgeList[1].rightVertex;
		if (c == a || c == b) {
			c = edgeList[1].leftVertex;
		}

		if (cross(a, c, b, c) < 0) {
			// if (!ac.newcross(bc))
			return true;
		} else {
			return false;
		}
	}

	// Return true if the triangle has become inverted or the area is zero
	@Override
	public boolean invertedOrZeroArea() {
		Vertex a, b, c;
		MyVector ac, bc;
		a = firstVertex;
		b = edgeList[0].otherVertex(a);
		c = edgeList[1].rightVertex;
		if (c == a || c == b) {
			c = edgeList[1].leftVertex;
		}

		if (cross(a, c, b, c) <= 0) {
			// if (!ac.newcross(bc))
			return true;
		} else {
			return false;
		}
	}

	// Return true if the triangle area is zero.
	public boolean zeroArea() {
		Vertex a, b, c;
		MyVector ac, bc;
		a = firstVertex;
		b = edgeList[0].otherVertex(a);
		c = edgeList[1].rightVertex;
		if (c == a || c == b) {
			c = edgeList[1].leftVertex;
		}

		if (cross(a, c, b, c) == 0) {
			return true;
		} else {
			return false;
		}
	}

	// Triangles don't have concavities, so return false.
	@Override
	public boolean concavityAt(Vertex n) {
		return false;
	}

	@Override
	public void replaceEdge(Edge e, Edge replacement) {
		edgeList[indexOf(e)] = replacement;
	}

	// Returns an edge that is on the boundary. Returns null if not found.
	public Edge getBoundaryEdge() {
		if (edgeList[0].element2 == null) {
			return edgeList[0];
		} else if (edgeList[1].element2 == null) {
			return edgeList[1];
		} else if (edgeList[2].element2 == null) {
			return edgeList[2];
		} else {
			return null;
		}
	}

	/**
	 * Update the distortion metric according to the article "An approach to
	 * Combined Laplacian and Optimization-Based Smoothing for Triangular,
	 * Quadrilateral and Quad-Dominant Meshes" by by Cannan, Tristano, and Staten.
	 * 
	 * @return negative values for inverted triangles, else positive. Equilateral
	 *         triangles will return the maximum value of 1.
	 */
	@Override
	public void updateDistortionMetric() {
		updateDistortionMetric(sqrt3x2);
	}

	/**
	 * See updateDistortionMetric().
	 * 
	 * @return negative values for inverted triangles, else positive.
	 */
	public void updateDistortionMetric(double factor) {
		Msg.debug("Entering Triangle.updateDistortionMetric(..)");
		double AB = edgeList[0].len, CB = edgeList[1].len, CA = edgeList[2].len;
		Vertex a = edgeList[2].commonVertex(edgeList[0]), b = edgeList[0].commonVertex(edgeList[1]), c = edgeList[2].commonVertex(edgeList[1]);
		MyVector vCA = new MyVector(c, a), vCB = new MyVector(c, b);

		double temp = factor * Math.abs(vCA.cross(vCB)) / (CA * CA + AB * AB + CB * CB);
		if (inverted()) {
			distortionMetric = -temp;
		} else {
			distortionMetric = temp;
		}
		Msg.debug("Leaving Triangle.updateDistortionMetric(..): " + distortionMetric);
	}

	// Return the size of the largest angle
	@Override
	public double largestAngle() {
		double cand = ang[0];
		if (ang[1] > cand) {
			cand = ang[1];
		}
		if (ang[2] > cand) {
			cand = ang[2];
		}
		return cand;
	}

	// Return the Vertex at the largest interior angle
	@Override
	public Vertex VertexAtLargestAngle() {
		Vertex candVertex = edgeList[0].leftVertex;
		double cand = ang[angleIndex(candVertex)], temp;

		temp = ang[angleIndex(edgeList[0].rightVertex)];
		if (temp > cand) {
			candVertex = edgeList[0].rightVertex;
			cand = temp;
		}
		temp = ang[angleIndex(oppositeOfEdge(edgeList[0]))];
		if (temp > cand) {
			candVertex = oppositeOfEdge(edgeList[0]);
		}
		return candVertex;
	}

	// Return the length of the longest Edge
	@Override
	public double longestEdgeLength() {
		double temp = Math.max(edgeList[0].len, edgeList[1].len);
		return Math.max(temp, edgeList[2].len);
	}

	// Return the longest Edge
	public Edge longestEdge() {
		Edge temp;
		if (edgeList[0].len > edgeList[1].len) {
			temp = edgeList[0];
		} else {
			temp = edgeList[1];
		}

		if (edgeList[2].len > temp.len) {
			return edgeList[2];
		} else {
			return temp;
		}
	}

	// Set the color of the edges to green.
	@Override
	public void markEdgesLegal() {
		edgeList[0].color = java.awt.Color.green;
		edgeList[1].color = java.awt.Color.green;
		edgeList[2].color = java.awt.Color.green;
	}

	// Set the color of the edges to red.
	@Override
	public void markEdgesIllegal() {
		edgeList[0].color = java.awt.Color.red;
		edgeList[1].color = java.awt.Color.red;
		edgeList[2].color = java.awt.Color.red;
	}

	@Override
	public String descr() {
		return "";
//		Vertex Vertex1, Vertex2, Vertex3;
//		Vertex1 = edgeList[0].leftVertex;
//		Vertex2 = edgeList[0].rightVertex;
//		Vertex3 = edgeList[1].rightVertex;
//		if (Vertex3 == Vertex1 || Vertex3 == Vertex2) {
//			Vertex3 = edgeList[1].leftVertex;
//		}
//
//		return Vertex1.descr() + ", " + Vertex2.descr() + ", " + Vertex3.descr();
	}

	@Override
	public void printMe() {
		if (inverted()) {
			System.out.println(descr() + ", inverted");
		} else {
			System.out.println(descr() + ", not inverted");
		}

	}
}
