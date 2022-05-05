package meshditor;

import java.util.ArrayList;
import java.util.List;

/**
 * A class holding information for quadrilaterals, and with methods for the
 * handling of issues regarding quads.
 */

public class Quad extends Element {

	public boolean isFake;

	/** Create ordinary quad */
	public Quad(Edge baseEdge, Edge leftEdge, Edge rightEdge, Edge topEdge) {
		edgeList = new Edge[4];
		ang = new double[4];

		if (rightEdge != topEdge) {
			isFake = false;
		} else {
			isFake = true;
		}

		edgeList[base] = baseEdge;
		edgeList[left] = leftEdge;
		edgeList[right] = rightEdge;
		edgeList[top] = topEdge;

		edgeList[base].len = edgeList[base].computeLength();
		edgeList[left].len = edgeList[left].computeLength();
		edgeList[right].len = edgeList[right].computeLength();
		if (!isFake) {
			edgeList[top].len = edgeList[top].computeLength();
		}

		firstVertex = baseEdge.leftVertex;
		if (inverted()) {
			firstVertex = baseEdge.rightVertex;
		}

		updateAngles();
	}

	/**
	 * Create ordinary quad. Use Edge e's two elements (Triangles) as basis. Not
	 * tested thoroughly!!!
	 */
	public Quad(Edge e) {
		isFake = false;
		edgeList = new Edge[4];
		ang = new double[4];

		Triangle t1 = (Triangle) e.element1;
		Triangle t2 = (Triangle) e.element2;

		Edge c11 = t1.otherEdge(e);
		Edge c12 = t1.otherEdge(e, c11);
		Edge c21, c22, temp1 = t2.otherEdge(e);
		if (c11.commonVertex(temp1) != null) {
			c21 = temp1;
			c22 = t2.otherEdge(e, temp1);
		} else {
			c21 = t2.otherEdge(e, temp1);
			c22 = temp1;
		}

		// The following assignments to the array ang[] is based on a
		// geometric sketch. The correctness of this should be tested.
		if (c11.leftVertex.hasEdge(c12)) {
			edgeList[base] = c11;
			edgeList[left] = c12;
			ang[0] = t1.angle(c11, c12);

			if (c11.rightVertex.hasEdge(c21)) {
				edgeList[right] = c21;
				edgeList[top] = c22;
				ang[1] = t1.angle(c11, e) + t2.angle(e, c21);
				ang[2] = t1.angle(c12, e) + t2.angle(e, c22);
				ang[3] = t2.angle(c21, c22);
			} else {
				edgeList[right] = c22;
				edgeList[top] = c21;
				ang[1] = t1.angle(c11, e) + t2.angle(e, c22);
				ang[2] = t1.angle(c12, e) + t2.angle(e, c21);
				ang[3] = t2.angle(c21, c22);
			}
		} else { // c11.rightVertex.hasEdge(c12)
			edgeList[base] = c11;
			edgeList[right] = c12;
			ang[1] = t1.angle(c11, c12);

			if (c11.leftVertex.hasEdge(c21)) {
				edgeList[left] = c21;
				edgeList[top] = c22;
				ang[0] = t1.angle(c11, e) + t2.angle(e, c21);
				ang[2] = t2.angle(c21, c22);
				ang[3] = t1.angle(c12, e) + t2.angle(e, c22);
			} else { // c11.leftVertex.hasEdge(c22)
				edgeList[left] = c22;
				edgeList[top] = c21;
				ang[0] = t1.angle(c11, e) + t2.angle(e, c22);
				ang[2] = t2.angle(c21, c22);
				ang[3] = t1.angle(c12, e) + t2.angle(e, c21);
			}
		}

		firstVertex = edgeList[base].leftVertex;
		if (inverted()) {
			firstVertex = edgeList[base].rightVertex;
		}
	}

	/**
	 * Create a fake quad with 3 different vertices and
	 * edgeList[top]==edgeList[right]
	 */
	public Quad(Triangle t) {
		isFake = true;
		edgeList = new Edge[4];
		ang = new double[4];

		edgeList[base] = t.edgeList[0];
		if (t.edgeList[1].hasVertex(edgeList[base].leftVertex)) {
			edgeList[left] = t.edgeList[1];
			edgeList[right] = t.edgeList[2];
		} else {
			edgeList[left] = t.edgeList[2];
			edgeList[right] = t.edgeList[1];
		}

		edgeList[top] = edgeList[right];

		firstVertex = t.firstVertex;

		ang[0] = t.ang[0];
		ang[1] = t.ang[1];
		ang[2] = t.ang[2];
	}

	/**
	 * Create a simple quad for testing purposes only (nextCCWVertex(),
	 * isStrictlyconvex()). Not tested thoroughly!!!
	 * 
	 * @param e  is a diagonal edge
	 * @param n1 the first or the other two vertices in the quad
	 * @param n2 the second of the other two vertices in the quad
	 */
	public Quad(Edge e, Vertex n1, Vertex n2) {
		isFake = false;
		edgeList = new Edge[4];

		edgeList[base] = new Edge(e.leftVertex, n1);
		edgeList[top] = new Edge(n2, e.rightVertex);

		if (edgeList[base].leftVertex == e.leftVertex) {
			edgeList[right] = new Edge(edgeList[base].rightVertex, e.rightVertex);
			edgeList[left] = new Edge(edgeList[base].leftVertex, edgeList[top].otherVertex(e.rightVertex));
		} else {
			edgeList[right] = new Edge(edgeList[base].rightVertex, edgeList[top].otherVertex(e.rightVertex));
			edgeList[left] = new Edge(edgeList[base].leftVertex, e.rightVertex);
		}

		firstVertex = edgeList[base].leftVertex;
		if (inverted()) {
			firstVertex = edgeList[base].rightVertex;
		}
	}

	/**
	 * Constructor to make life easier for elementWithExchangedVertex(..) Create
	 * fake quad with only three vertices
	 */
	private Quad(Vertex n1, Vertex n2, Vertex n3, Vertex f) {
		isFake = true;
		edgeList = new Edge[4];
		ang = new double[4];

		edgeList[base] = new Edge(n1, n2);
		if (edgeList[base].leftVertex == n1) {
			edgeList[left] = new Edge(n1, n3);
			edgeList[right] = new Edge(n2, n3);
		} else {
			edgeList[left] = new Edge(n2, n3);
			edgeList[right] = new Edge(n1, n3);
		}
		edgeList[top] = edgeList[right];

		firstVertex = f;
		updateAngles();
	}

	/** Constructor to make life easier for elementWithExchangedVertex(..) */
	private Quad(Vertex n1, Vertex n2, Vertex n3, Vertex n4, Vertex f) {
		isFake = false;
		edgeList = new Edge[4];
		ang = new double[4];

		edgeList[base] = new Edge(n1, n2);
		edgeList[top] = new Edge(n3, n4);
		if (edgeList[base].leftVertex == n1) {
			edgeList[left] = new Edge(n1, n3);
			edgeList[right] = new Edge(n2, n4);
		} else {
			edgeList[left] = new Edge(n2, n4);
			edgeList[right] = new Edge(n1, n3);
		}

		firstVertex = f;
		updateAngles();
	}

	/**
	 * Create a simple quad for testing purposes only (constrainedLaplacianSmooth())
	 * Not tested thoroughly!!!
	 */
	@Override
	public Element elementWithExchangedVertices(Vertex original, Vertex replacement) {
		Vertex Vertex1 = edgeList[base].leftVertex;
		Vertex Vertex2 = edgeList[base].rightVertex;
		Vertex Vertex3 = edgeList[left].otherVertex(edgeList[base].leftVertex);

		if (isFake) {
			if (Vertex1 == original) {
				if (original == firstVertex) {
					return new Quad(replacement, Vertex2, Vertex3, replacement);
				} else {
					return new Quad(replacement, Vertex2, Vertex3, firstVertex);
				}
			} else if (Vertex2 == original) {
				if (original == firstVertex) {
					return new Quad(Vertex1, replacement, Vertex3, replacement);
				} else {
					return new Quad(Vertex1, replacement, Vertex3, firstVertex);
				}
			} else if (Vertex3 == original) {
				if (original == firstVertex) {
					return new Quad(Vertex1, Vertex2, replacement, replacement);
				} else {
					return new Quad(Vertex1, Vertex2, replacement, firstVertex);
				}
			} else {
				return null;
			}
		}

		Vertex Vertex4 = edgeList[right].otherVertex(edgeList[base].rightVertex);

		if (Vertex1 == original) {
			if (original == firstVertex) {
				return new Quad(replacement, Vertex2, Vertex3, Vertex4, replacement);
			} else {
				return new Quad(replacement, Vertex2, Vertex3, Vertex4, firstVertex);
			}
		} else if (Vertex2 == original) {
			if (original == firstVertex) {
				return new Quad(Vertex1, replacement, Vertex3, Vertex4, replacement);
			} else {
				return new Quad(Vertex1, replacement, Vertex3, Vertex4, firstVertex);
			}
		} else if (Vertex3 == original) {
			if (original == firstVertex) {
				return new Quad(Vertex1, Vertex2, replacement, Vertex4, replacement);
			} else {
				return new Quad(Vertex1, Vertex2, replacement, Vertex4, firstVertex);
			}
		} else if (Vertex4 == original) {
			if (original == firstVertex) {
				return new Quad(Vertex1, Vertex2, Vertex3, replacement, replacement);
			} else {
				return new Quad(Vertex1, Vertex2, Vertex3, replacement, firstVertex);
			}
		} else {
			return null;
		}
	}

	/**
	 * @return true if the quad becomes inverted when Vertex n1 is relocated to pos.
	 *         n2. Else return false.
	 */
	@Override
	public boolean invertedWhenVertexRelocated(Vertex n1, Vertex n2) {
		Vertex thisFirstVertex = firstVertex;
		Vertex a = edgeList[base].leftVertex, b = edgeList[base].rightVertex, c = edgeList[right].otherVertex(b),
				d = edgeList[left].otherVertex(a);

		if (a == n1) {
			a = n2;
		} else if (b == n1) {
			b = n2;
		} else if (c == n1) {
			c = n2;
		} else if (d == n1) {
			d = n2;
		}

		if (n1 == firstVertex) {
			thisFirstVertex = n2;
		}

		// We need at least 3 okays to be certain that this quad is not inverted
		int okays = 0;
		if (cross(a, d, b, d) > 0) {
			okays++;
		}
		if (cross(a, c, b, c) > 0) {
			okays++;
		}
		if (cross(c, a, d, a) > 0) {
			okays++;
		}
		if (cross(c, b, d, b) > 0) {
			okays++;
		}

		if (b == thisFirstVertex) {
			okays = 4 - okays;
		}

		if (okays >= 3) {
			return false;
		} else {
			return true;
		}
	}

	/**
	 * Test whether any neighboring elements becomes inverted if the quad is
	 * collapsed in a particular manner.
	 * 
	 * @param q     the quad to be collapsed
	 * @param n     a Vertex holding the position for where the joined vertices are
	 *              to be located
	 * @param n1    the Vertex in quad q that is to be joined with opposite Vertex
	 *              n2
	 * @param n2    the Vertex in quad q that is to be joined with opposite Vertex
	 *              n1
	 * @param list1 the list of elements adjacent n1
	 * @param list2 the list of elements adjacent n2
	 * @return true if any elements adjacent to quad q becomes inverted when
	 *         collapsing quad q, joining its two opposite vertices n1 and n2 to the
	 *         position held by Vertex n. Vertex n must be located somewhere inside
	 *         quad q.
	 */
	public boolean anyInvertedElementsWhenCollapsed(Vertex n, Vertex n1, Vertex n2, List<Element> list1, List<Element> list2) {
		for (Element elem : list1) {
			if (elem != this && elem.invertedWhenVertexRelocated(n1, n)) {
				return true;
			}
		}

		for (Element elem : list1) {
			if (elem != this && elem.invertedWhenVertexRelocated(n2, n)) {
				return true;
			}
		}

		return false;
	}
	
	@Override
	public int hashCode() {
		int hash = edgeList[base].leftVertex.hashCode();
		hash ^= edgeList[base].rightVertex.hashCode();
		hash ^= edgeList[top].leftVertex.hashCode();
		hash ^= edgeList[top].rightVertex.hashCode();
		return hash;
	}

	@Override
	public boolean equals(Object o) {
		if (o instanceof Quad) {
			Quad q = (Quad) o;
			Vertex Vertex1, Vertex2, Vertex3, Vertex4;
			Vertex qVertex1, qVertex2, qVertex3, qVertex4;
			Vertex1 = edgeList[base].leftVertex;
			Vertex2 = edgeList[base].rightVertex;
			Vertex3 = edgeList[left].otherVertex(edgeList[base].leftVertex);
			Vertex4 = edgeList[right].otherVertex(edgeList[base].rightVertex);

			qVertex1 = q.edgeList[base].leftVertex;
			qVertex2 = q.edgeList[base].rightVertex;
			qVertex3 = q.edgeList[left].otherVertex(q.edgeList[base].leftVertex);
			qVertex4 = q.edgeList[right].otherVertex(q.edgeList[base].rightVertex);

			if (Vertex1 == qVertex1 && Vertex2 == qVertex2 && Vertex3 == qVertex3 && Vertex4 == qVertex4) {
				return true;
			} else {
				return false;
			}
		} else {
			return false;
		}
	}

	/** @return edge's index in this quad's edgeList */
	@Override
	public int indexOf(Edge e) {
		if (edgeList[base] == e) {
			return base;
		} else if (edgeList[left] == e) {
			return left;
		} else if (edgeList[right] == e) {
			return right;
		} else if (edgeList[top] == e) {
			return top;
		} else {
			return -1;
		}
	}

	/**
	 * Collapse the quad by welding together two neighboring edges. For each edge in
	 * the edgeList of the otherVertex (otherE2) of edge e2: Provided that otherE1
	 * or the otherVertex of otherE2 is not found in any edge in the edgeList of the
	 * otherVertex of edge e1 (otherE1), then the edge is added to otherE1's
	 * edgeList. Assumes that the inner quad area is empty and that the edges of the
	 * quad don't have the quad as an element.
	 * 
	 * @param e1 first edge
	 * @param e2 second edge, and the neighbor to edge e1 in this quad
	 */
	public void closeQuad(Edge e1, Edge e2) {
		Vertex nK = e1.commonVertex(e2);
		Vertex nKp1 = e1.otherVertex(nK), nKm1 = e2.otherVertex(nK), other;
		Quad q;
		boolean found = false;
		Edge e, eI, eJ;
		ArrayList<Edge> addList = new ArrayList<>();
		ArrayList<Element> quadList = new ArrayList<>();

		for (int i = 0; i < nKm1.edgeList.size(); i++) {
			eI = nKm1.edgeList.get(i);
			other = eI.otherVertex(nKm1);
			for (int j = 0; j < nKp1.edgeList.size(); j++) {
				eJ = nKp1.edgeList.get(j);

				if (other == eJ.otherVertex(nKp1)) {
					found = true;
					other.edgeList.remove(other.edgeList.indexOf(eI));

					if (eI.element1 != null) {
						if (eI.element1.firstVertex == nKm1) {
							eI.element1.firstVertex = nKp1; // Don't forget firstVertex!!
						}
						eI.element1.replaceEdge(eI, eJ);
						eJ.connectToElement(eI.element1);
						if (eI.element1 instanceof Quad) {
							quadList.add(eI.element1); // but that must be done later
						}
					}
					break;
				}
			}

			if (!found) {
				if (eI.element1 != null && eI.element1.firstVertex == nKm1) {
					eI.element1.firstVertex = nKp1; // Don't forget firstVertex!!
				}
				if (eI.element2 != null && eI.element2.firstVertex == nKm1) {
					eI.element2.firstVertex = nKp1; // Don't forget firstVertex!!
				}

				nKm1.edgeList.set(i, null);
				eI.replaceVertex(nKm1, nKp1);
				addList.add(eI); // nKp1.edgeList.add(eI);
			} else {
				found = false;
			}
		}

		for (Edge element : addList) {
			e = element;
			nKp1.edgeList.add(e);
		}

		for (Object element : quadList) {
			q = (Quad) element;
			q.updateLR();
		}

		nKm1.edgeList.clear();
	}

	/**
	 * Create a new quad by combining this quad with quad q. The two quads must
	 * initially share two incident edges:
	 * 
	 * @param e1 first common edge
	 * @param e2 second common edge
	 */
	public Quad combine(Quad q, Edge e1, Edge e2) {
		Quad quad;
		Edge e, edges[] = new Edge[4];
		int i;

		edges[0] = null;
		for (i = 0; i < 4; i++) {
			e = edgeList[i];

			if (e != e1 && e != e2) {
				if (edges[0] == null) {
					edges[0] = e;
				} else if (e.hasVertex(edges[0].leftVertex)) {
					edges[1] = e;
				} else if (e.hasVertex(edges[0].rightVertex)) {
					edges[2] = e;
				}
			}
		}

		for (i = 0; i < 4; i++) {
			if (q.isFake && i == 3) {
				edges[3] = edges[2];
			} else {
				e = q.edgeList[i];

				if (e != e1 && e != e2) {
					if (e.hasVertex(edges[0].leftVertex)) {
						edges[1] = e;
					} else if (e.hasVertex(edges[0].rightVertex)) {
						edges[2] = e;
					} else {
						edges[3] = e;
					}
				}
			}
		}

		// A triangle (fake quad) will have edges[2]== edges[3].
		quad = new Quad(edges[0], edges[1], edges[2], edges[3]);
		return quad;
	}

	/**
	 * Create a new triangle by combining this quad with Triangle t. The two
	 * elements must initially share two incident edges:
	 * 
	 * @param e1 first common edge
	 * @param e2 second common edge
	 */
	public Triangle combine(Triangle t, Edge e1, Edge e2) {
		Triangle tri;
		Edge e, edges[] = new Edge[3];
		int i;

		edges[0] = null;
		for (i = 0; i < 4; i++) {
			e = edgeList[i];

			if (e != e1 && e != e2) {
				if (edges[0] == null) {
					edges[0] = e;
				} else if (e.hasVertex(edges[0].leftVertex)) {
					edges[1] = e;
				} else {
					edges[2] = e;
				}
			}
		}

		for (i = 0; i < 3; i++) {
			e = t.edgeList[i];

			if (e != e1 && e != e2) {
				if (e.hasVertex(edges[0].leftVertex)) {
					edges[1] = e;
				} else {
					edges[2] = e;
				}
			}
		}

		tri = new Triangle(edges[0], edges[1], edges[2]);
		return tri;
	}

	/** @return neighbor element sharing edge e */
	@Override
	public Element neighbor(Edge e) {
		if (e.element1 == this) {
			return e.element2;
		} else if (e.element2 == this) {
			return e.element1;
		} else {
			Msg.warning("Quad.neighbor(Edge): returning null");
			return null;
		}
	}

	/**
	 * @return the edge in this quad that is a neighbor of the given Vertex and
	 *         edge.
	 */
	@Override
	public Edge neighborEdge(Vertex n, Edge e) {
		if (edgeList[base].leftVertex == n) {
			if (edgeList[base] == e) {
				return edgeList[left];
			} else if (edgeList[left] == e) {
				return edgeList[base];
			} else {
				Msg.warning("Quad.neighborEdge(Vertex n, Edge e): neighbor not found:");
				Msg.error("quad: " + descr() + ", n: " + n.descr() + ", e: " + e.descr());
				return null;
			}
		} else if (edgeList[base].rightVertex == n) {
			if (edgeList[base] == e) {
				return edgeList[right];
			} else if (edgeList[right] == e) {
				return edgeList[base];
			} else {
				Msg.warning("Quad.neighborEdge(Vertex n, Edge e): neighbor not found:");
				Msg.error("quad: " + descr() + ", n: " + n.descr() + ", e: " + e.descr());
				return null;
			}
		} else if (edgeList[left].otherVertex(edgeList[base].leftVertex) == n) {
			if (isFake) {
				if (edgeList[right] == e) {
					return edgeList[left];
				} else if (edgeList[left] == e) {
					return edgeList[right];
				} else {
					Msg.warning("Quad.neighborEdge(Vertex n, Edge e): neighbor not found:");
					Msg.error("quad: " + descr() + ", n: " + n.descr() + ", e: " + e.descr());
					return null;
				}
			} else if (edgeList[top] == e) {
				return edgeList[left];
			} else if (edgeList[left] == e) {
				return edgeList[top];
			} else {
				Msg.warning("Quad.neighborEdge(Vertex n, Edge e): neighbor not found:");
				Msg.error("quad: " + descr() + ", n: " + n.descr() + ", e: " + e.descr());
				return null;
			}
		} else if (edgeList[right].otherVertex(edgeList[base].rightVertex) == n) {
			if (edgeList[top] == e) {
				return edgeList[right];
			} else if (edgeList[right] == e) {
				return edgeList[top];
			} else {
				Msg.warning("Quad.neighborEdge(Vertex n, Edge e): neighbor not found:");
				Msg.error("quad: " + descr() + ", n: " + n.descr() + ", e: " + e.descr());
				return null;
			}
		} else {
			Msg.warning("Quad.neighborEdge(Vertex n, Edge e): neighbor not found:");
			Msg.error("quad: " + descr() + ", n: " + n.descr() + ", e: " + e.descr());
			return null;
		}
	}

	/** @return an edge in this quad that is a neighbor of the given Vertex. */
	public Edge neighborEdge(Vertex n) {
		if (edgeList[base].leftVertex == n) {
			return edgeList[base];
		} else if (edgeList[base].rightVertex == n) {
			return edgeList[right];
		} else if (edgeList[left].otherVertex(edgeList[base].leftVertex) == n) {
			return edgeList[left];
		} else if (edgeList[right].otherVertex(edgeList[base].rightVertex) == n) {
			return edgeList[top];
		} else {
			Msg.error("Quad.neighborEdge(Vertex n): neighbor not found.");
			return null;
		}
	}

	@Override
	public double angle(Edge e, Vertex n) {
		// Find this edge's index
		int thisEdgeIndex = indexOf(e);

		// Find other edge's index
		int otherEdgeIndex;
		if (edgeList[base] != e && edgeList[base].hasVertex(n)) {
			otherEdgeIndex = base;
		} else if (edgeList[left] != e && edgeList[left].hasVertex(n)) {
			otherEdgeIndex = left;
		} else if (edgeList[right] != e && edgeList[right].hasVertex(n)) {
			otherEdgeIndex = right;
		} else {
			otherEdgeIndex = top;
		}

		// Return correct angle
		return ang[angleIndex(thisEdgeIndex, otherEdgeIndex)];
	}

	public int angleIndex(int e1Index, int e2Index) {
		if ((e1Index == base && e2Index == left) || // angle at base, left
				(e1Index == left && e2Index == base)) {
			return 0;
		} else if ((e1Index == base && e2Index == right) || // angle at base,right
				(e1Index == right && e2Index == base)) {
			return 1;
		} else if ((e1Index == top && e2Index == left) || // angle at top, left
				(e1Index == left && e2Index == top)) {
			return 2;
		} else {
			return 3;
		}
	}

	@Override
	public int angleIndex(Edge e1, Edge e2) {
		return angleIndex(indexOf(e1), indexOf(e2));
	}

	@Override
	public int angleIndex(Vertex n) {
		if (edgeList[base].leftVertex == n) {
			return 0;
		} else if (edgeList[base].rightVertex == n) {
			return 1;
		} else if (edgeList[left].otherVertex(edgeList[base].leftVertex) == n) {
			return 2;
		} else if (edgeList[right].otherVertex(edgeList[base].rightVertex) == n) {
			return 3;
		} else {
			Msg.error("Quad.angleIndex(Vertex): Vertex not found");
			return -1;
		}
	}

	@Override
	public double angle(Edge e1, Edge e2) {
		return ang[angleIndex(indexOf(e1), indexOf(e2))];
	}

	@Override
	public boolean hasVertex(Vertex n) {
		if (edgeList[base].leftVertex.equals(n) || edgeList[base].rightVertex.equals(n) || edgeList[top].leftVertex.equals(n)
				|| edgeList[top].rightVertex.equals(n)) {
			return true;
		} else {
			return false;
		}
	}

	@Override
	public boolean hasEdge(Edge e) {
		if (edgeList[base] == e || edgeList[left] == e || edgeList[right] == e || edgeList[top] == e) {
			return true;
		} else {
			return false;
		}
	}

	public boolean boundaryQuad() {
		if (neighbor(edgeList[base]) instanceof Quad || neighbor(edgeList[left]) instanceof Quad
				|| neighbor(edgeList[right]) instanceof Quad || neighbor(edgeList[top]) instanceof Quad) {
			return true;
		} else {
			return false;
		}
	}

	/**
	 * Method to verify that this quad is a boundary diamond. A boundary diamond is
	 * defined as a quad with only one Vertex on the boundary.
	 */
	public boolean boundaryDiamond() {
		int bVertices = 0;
		if (edgeList[base].leftVertex.boundaryVertex()) {
			bVertices++;
		}
		if (edgeList[base].rightVertex.boundaryVertex()) {
			bVertices++;
		}
		if (edgeList[top].leftVertex.boundaryVertex()) {
			bVertices++;
		}
		if (edgeList[top].rightVertex.boundaryVertex()) {
			bVertices++;
		}

		if (bVertices == 1) {
			return true;
		} else {
			return false;
		}
	}

	/** Get any Vertex on the boundary that belongs to this quad. */
	public Vertex getBoundaryVertex() {
		if (edgeList[base].leftVertex.boundaryVertex()) {
			return edgeList[base].leftVertex;
		} else if (edgeList[base].rightVertex.boundaryVertex()) {
			return edgeList[base].rightVertex;
		} else if (edgeList[top].leftVertex.boundaryVertex()) {
			return edgeList[top].leftVertex;
		} else if (edgeList[top].rightVertex.boundaryVertex()) {
			return edgeList[top].rightVertex;
		} else {
			return null;
		}
	}

	/**
	 * Method to verify that the quad has an area greater than 0. We simply check
	 * that the vertices of the element are not colinear.
	 */
	@Override
	public boolean areaLargerThan0() {
		Vertex Vertex3 = edgeList[top].leftVertex;
		Vertex Vertex4 = edgeList[top].rightVertex;

		if (!Vertex3.onLine(edgeList[base]) || !Vertex4.onLine(edgeList[base])) {
			return true;
		}
		return true;
	}

	/** Method to verify that the quad is convex. */
	public boolean isConvex() {
		Vertex n1 = edgeList[base].leftVertex;
		Vertex n2 = edgeList[base].rightVertex;
		Vertex n3 = edgeList[left].otherVertex(n1);
		Vertex n4 = edgeList[right].otherVertex(n2);

		MyVector d1 = new MyVector(n1, n4);
		MyVector d2 = new MyVector(n2, n3);

		return d1.pointIntersects(d2);
	}

	/**
	 * Method to verify that the quad is strictly convex, that is, convex in the
	 * common sense and in addition demanding that no three Vertices are colinear.
	 */
	public boolean isStrictlyConvex() {
		Vertex n1 = edgeList[base].leftVertex;
		Vertex n2 = edgeList[base].rightVertex;
		Vertex n3 = edgeList[left].otherVertex(n1);
		Vertex n4 = edgeList[right].otherVertex(n2);

		MyVector d1 = new MyVector(n1, n4);
		MyVector d2 = new MyVector(n2, n3);
		if (d1.innerpointIntersects(d2)) {
			return true;
		} else {
			return false;
		}
	}

	/**
	 * @return true if the quad is a bowtie, defined as a quad with two opposite
	 *         edges that intersect.
	 */
	public boolean isBowtie() {
		Vertex n1 = edgeList[base].leftVertex;
		Vertex n2 = edgeList[base].rightVertex;
		Vertex n3 = edgeList[left].otherVertex(n1);
		Vertex n4 = edgeList[right].otherVertex(n2);

		MyVector d1 = new MyVector(n1, n2);
		MyVector d2 = new MyVector(n3, n4);
		MyVector d3 = new MyVector(n1, n3);
		MyVector d4 = new MyVector(n2, n4);
		if (d1.pointIntersects(d2) || d3.pointIntersects(d4)) {
			return true;
		} else {
			return false;
		}
	}

	/**
	 * @return true if the quad is a chevron, defined as a quad with a greatest
	 *         angle that is greater than 200 degrees.
	 */
	public boolean isChevron() {
		if (largestAngle() >= CHEVRONMIN) {
			return true;
		} else {
			return false;
		}
	}

	/**
	 * @return true if the largest angle of the quad is greater than 180 degrees.
	 */
	public boolean largestAngleGT180() {
		if (largestAngle() > DEG_180) {
			return true;
		} else {
			return false;
		}
	}

	/**
	 * Not tested much yet, but should work very well in principle.
	 * 
	 * @return the next Vertex in the ccw direction around this quad.
	 */
	public Vertex nextCCWVertex(Vertex n) {
		Vertex n2, n3, n4;
		Edge e2, e3, e4;
		if (n == firstVertex) {
			return edgeList[base].otherVertex(firstVertex); // n2
		} else {
			n2 = edgeList[base].otherVertex(firstVertex);
			e2 = neighborEdge(n2, edgeList[base]);
			if (n == n2) {
				return e2.otherVertex(n2); // n3
			} else {
				n3 = e2.otherVertex(n2);
				e3 = neighborEdge(n3, e2);
				if (n == n3) {
					return e3.otherVertex(n3); // n4
				} else {
					n4 = e3.otherVertex(n3);
					if (n == n4) {
						return firstVertex;
					} else {
						return null;
					}
				}
			}
		}
	}

	/**
	 * Update so that the edge connected to edgeList[base].leftVertex is
	 * edgeList[left] and that the edge connected to edgeList[base].rightVertex is
	 * edgeList[right]. The angle between base and left is at pos 0 in the ang
	 * array. The angle between right and base is at pos 1 in the ang array. The
	 * angle between left and top is at pos 2 in the ang array. The angle between
	 * top and right is at pos 3 in the ang array.
	 */
	public void updateLR() {
		Edge temp;
		double dt0, dt1, dt2, dt3;
		if (!edgeList[left].hasVertex(edgeList[base].leftVertex)) {
			temp = edgeList[left];
			edgeList[left] = edgeList[right];
			edgeList[right] = temp;

			dt0 = ang[0];
			dt1 = ang[1];
			dt2 = ang[2];
			dt3 = ang[3];

			ang[0] = dt1;
			ang[1] = dt0;
			ang[2] = dt3;
			ang[3] = dt2;
		}
	}

	/**
	 * Update the values in the ang array. Works correctly only for uninverted
	 * quads.
	 */
	@Override
	public void updateAngles() {
		if (isFake) {
			if (firstVertex == edgeList[base].leftVertex) {
				ang[0] = edgeList[base].computeCCWAngle(edgeList[left]);
				ang[1] = edgeList[right].computeCCWAngle(edgeList[base]);
				ang[2] = edgeList[left].computeCCWAngle(edgeList[right]);
			} else {
				ang[0] = edgeList[left].computeCCWAngle(edgeList[base]);
				ang[1] = edgeList[base].computeCCWAngle(edgeList[right]);
				ang[2] = edgeList[right].computeCCWAngle(edgeList[left]);
			}
		} else if (firstVertex == edgeList[base].leftVertex) {
			ang[0] = edgeList[base].computeCCWAngle(edgeList[left]);
			ang[1] = edgeList[right].computeCCWAngle(edgeList[base]);
			ang[2] = edgeList[left].computeCCWAngle(edgeList[top]);
			ang[3] = edgeList[top].computeCCWAngle(edgeList[right]);
		} else {
			ang[0] = edgeList[left].computeCCWAngle(edgeList[base]);
			ang[1] = edgeList[base].computeCCWAngle(edgeList[right]);
			ang[2] = edgeList[top].computeCCWAngle(edgeList[left]);
			ang[3] = edgeList[right].computeCCWAngle(edgeList[top]);
		}
	}

	/** Update the values in the ang array except at the specified Vertex. */
	public void updateAnglesExcept(Vertex n) {
		int i = angleIndex(n);
		if (isFake) {
			if (firstVertex == edgeList[base].leftVertex) {
				if (i != 0) {
					ang[0] = edgeList[base].computeCCWAngle(edgeList[left]);
				}
				if (i != 1) {
					ang[1] = edgeList[right].computeCCWAngle(edgeList[base]);
				}
				if (i != 2) {
					ang[2] = edgeList[left].computeCCWAngle(edgeList[right]);
				}
			} else {
				if (i != 0) {
					ang[0] = edgeList[left].computeCCWAngle(edgeList[base]);
				}
				if (i != 1) {
					ang[1] = edgeList[base].computeCCWAngle(edgeList[right]);
				}
				if (i != 2) {
					ang[2] = edgeList[right].computeCCWAngle(edgeList[left]);
				}
			}
		} else if (firstVertex == edgeList[base].leftVertex) {
			if (i != 0) {
				ang[0] = edgeList[base].computeCCWAngle(edgeList[left]);
			}
			if (i != 1) {
				ang[1] = edgeList[right].computeCCWAngle(edgeList[base]);
			}
			if (i != 2) {
				ang[2] = edgeList[left].computeCCWAngle(edgeList[top]);
			}
			if (i != 3) {
				ang[3] = edgeList[top].computeCCWAngle(edgeList[right]);
			}
		} else {
			if (i != 0) {
				ang[0] = edgeList[left].computeCCWAngle(edgeList[base]);
			}
			if (i != 1) {
				ang[1] = edgeList[base].computeCCWAngle(edgeList[right]);
			}
			if (i != 2) {
				ang[2] = edgeList[top].computeCCWAngle(edgeList[left]);
			}
			if (i != 3) {
				ang[3] = edgeList[right].computeCCWAngle(edgeList[top]);
			}
		}
	}

	@Override
	public void updateAngle(Vertex n) {
		int i = angleIndex(n);

		if (isFake) {
			if (firstVertex == edgeList[base].leftVertex) {
				if (i == 0) {
					ang[0] = edgeList[base].computeCCWAngle(edgeList[left]);
				} else if (i == 1) {
					ang[1] = edgeList[right].computeCCWAngle(edgeList[base]);
				} else if (i == 2) {
					ang[2] = edgeList[left].computeCCWAngle(edgeList[right]);
				}
			} else if (i == 0) {
				ang[0] = edgeList[left].computeCCWAngle(edgeList[base]);
			} else if (i == 1) {
				ang[1] = edgeList[base].computeCCWAngle(edgeList[right]);
			} else if (i == 2) {
				ang[2] = edgeList[right].computeCCWAngle(edgeList[left]);
			}
		} else if (firstVertex == edgeList[base].leftVertex) {
			if (i == 0) {
				ang[0] = edgeList[base].computeCCWAngle(edgeList[left]);
			} else if (i == 1) {
				ang[1] = edgeList[right].computeCCWAngle(edgeList[base]);
			} else if (i == 2) {
				ang[2] = edgeList[left].computeCCWAngle(edgeList[top]);
			} else if (i == 3) {
				ang[3] = edgeList[top].computeCCWAngle(edgeList[right]);
			}
		} else if (i == 0) {
			ang[0] = edgeList[left].computeCCWAngle(edgeList[base]);
		} else if (i == 1) {
			ang[1] = edgeList[base].computeCCWAngle(edgeList[right]);
		} else if (i == 2) {
			ang[2] = edgeList[top].computeCCWAngle(edgeList[left]);
		} else if (i == 3) {
			ang[3] = edgeList[right].computeCCWAngle(edgeList[top]);
		}
	}

	/** Method to test whether the quad is inverted. */
	@Override
	public boolean inverted() {
		if (isFake) {

			Vertex a, b, c;
			a = firstVertex;
			b = edgeList[base].otherVertex(a);
			if (a == edgeList[base].leftVertex) {
				c = edgeList[left].otherVertex(a);
			} else {
				c = edgeList[right].otherVertex(a);
			}

			if (cross(a, c, b, c) < 0) {
				return true;
			} else {
				return false;
			}
		}

		Vertex a = edgeList[base].leftVertex, b = edgeList[base].rightVertex, c = edgeList[right].otherVertex(b),
				d = edgeList[left].otherVertex(a);

		// We need at least 3 okays to be certain that this quad is not inverted
		int okays = 0;
		if (cross(a, d, b, d) > 0) {
			okays++;
		}
		if (cross(a, c, b, c) > 0) {
			okays++;
		}
		if (cross(c, a, d, a) > 0) {
			okays++;
		}
		if (cross(c, b, d, b) > 0) {
			okays++;
		}

		if (b == firstVertex) {
			okays = 4 - okays;
		}

		if (okays >= 3) {
			return false;
		} else {
			return true;
		}
	}

	/** Method to test whether the quad is inverted or its area is zero. */
	@Override
	public boolean invertedOrZeroArea() {
		if (isFake) {

			Vertex a, b, c;
			a = firstVertex;
			b = edgeList[base].otherVertex(a);
			if (a == edgeList[base].leftVertex) {
				c = edgeList[left].otherVertex(a);
			} else {
				c = edgeList[right].otherVertex(a);
			}

			if (cross(a, c, b, c) <= 0) {
				return true;
			} else {
				return false;
			}
		}

		Vertex a = edgeList[base].leftVertex, b = edgeList[base].rightVertex, c = edgeList[right].otherVertex(b),
				d = edgeList[left].otherVertex(a);

		// We need at least 3 okays to be certain that this quad is not inverted
		int okays = 0;
		if (cross(a, d, b, d) >= 0) {
			okays++;
		}
		if (cross(a, c, b, c) >= 0) {
			okays++;
		}
		if (cross(c, a, d, a) >= 0) {
			okays++;
		}
		if (cross(c, b, d, b) >= 0) {
			okays++;
		}

		if (b == firstVertex) {
			okays = 4 - okays;
		}

		if (okays >= 3) {
			return false;
		} else {
			return true;
		}
	}

	/**
	 * Determines whether there is a concavity (angle > 180 degrees) at the
	 * specified Vertex.
	 * 
	 * @param n the Vertex at the angle to investigate
	 * @return true if the element has a concavity at the specified Vertex.
	 */
	@Override
	public boolean concavityAt(Vertex n) {
		if (ang[angleIndex(n)] >= Math.PI) {
			return true;
		} else {
			return false;
		}
	}

	/**
	 * @return the centroid of this quad.... or at least a point *inside* the
	 *         quad... Assumes that the quad is not inverted.
	 */
	public Vertex centroid() {
		double x = 0, y = 0;
		Vertex bL = edgeList[base].leftVertex, bR = edgeList[base].rightVertex;
		Vertex uL = edgeList[left].otherVertex(bL), uR = edgeList[right].otherVertex(bR);

		if (isFake) {
			x = (bL.x + uL.x + bR.x) / 3.0;
			y = (bL.y + uL.y + bR.y) / 3.0;
			return new Vertex(x, y);
		} else {
			if (concavityAt(bL)) {
				x = (bL.x + uR.x) / 2.0;
				y = (bL.y + uR.y) / 2.0;
			} else if (concavityAt(bR)) {
				x = (bR.x + uL.x) / 2.0;
				y = (bR.y + uL.y) / 2.0;
			} else if (concavityAt(uL)) {
				x = (uL.x + bR.x) / 2.0;
				y = (uL.y + bR.y) / 2.0;
			} else if (concavityAt(uR)) {
				x = (uR.x + bL.x) / 2.0;
				y = (uR.y + bL.y) / 2.0;
			} else {
				x = (bL.x + uL.x + uR.x + bR.x) / 4.0;
				y = (bL.y + uL.y + uR.y + bR.y) / 4.0;

			}
			return new Vertex(x, y);
		}
	}

	/**
	 * @param n a Vertex in this quad
	 * @return the Vertex on the opposite side of Vertex n in the quad
	 */
	public Vertex oppositeVertex(Vertex n) {
		// 2 out of 4 edges has Vertex n, so at least 1 out of 3 edge must have it, too:
		Edge startEdge;
		if (edgeList[0].hasVertex(n)) {
			startEdge = edgeList[0];
		} else if (edgeList[1].hasVertex(n)) {
			startEdge = edgeList[1];
		} else if (edgeList[2].hasVertex(n)) {
			startEdge = edgeList[2];
		} else {
			return null; // Most likely, Vertex n is not part of this Quad.
		}

		Vertex n2 = startEdge.otherVertex(n);
		Edge e = neighborEdge(n2, startEdge);
		return e.otherVertex(n2);
	}

	/**
	 * @return the opposite Edge of Vertex n that is cw to the other opposite Edge
	 */
	public Edge cwOppositeEdge(Vertex n) {
		if (n == edgeList[base].leftVertex) {
			return edgeList[right];
		} else if (n == edgeList[base].rightVertex) {
			return edgeList[top];
		} else if (n == edgeList[left].otherVertex(edgeList[base].leftVertex)) {
			return edgeList[base];
		} else if (n == edgeList[right].otherVertex(edgeList[base].rightVertex)) {
			return edgeList[left];
		} else {
			return null;
		}
	}

	/**
	 * @return the opposite Edge of Vertex n that is ccw to the other opposite Edge
	 */
	public Edge ccwOppositeEdge(Vertex n) {
		if (n == edgeList[base].leftVertex) {
			return edgeList[top];
		} else if (n == edgeList[base].rightVertex) {
			return edgeList[left];
		} else if (n == edgeList[left].otherVertex(edgeList[base].leftVertex)) {
			return edgeList[right];
		} else if (n == edgeList[right].otherVertex(edgeList[base].rightVertex)) {
			return edgeList[base];
		} else {
			return null;
		}
	}

	public Edge oppositeEdge(Edge e) {
		if (e == edgeList[base]) {
			return edgeList[top];
		} else if (e == edgeList[left]) {
			return edgeList[right];
		} else if (e == edgeList[right]) {
			return edgeList[left];
		} else if (e == edgeList[top]) {
			return edgeList[base];
		} else {
			return null;
		}
	}

	/**
	 * Check to see if any of the neighboring quad elements have become inverted
	 * NOTE 1: I might not need to check those elements that lies behind the front.
	 */
	public boolean invertedNeighbors() {
		Vertex uLVertex = edgeList[left].otherVertex(edgeList[base].leftVertex);
		Vertex uRVertex = edgeList[right].otherVertex(edgeList[base].rightVertex);
		Element curElem;
		Edge curEdge;

		// Parse all adjacent elements at upper left Vertex from neigbor of left edge
		// to,
		// but not including, neighbor of top edge.
		curElem = neighbor(edgeList[left]);
		curEdge = edgeList[left];
		while (curElem != null && curEdge != edgeList[top]) {
			if (curElem.inverted()) {
				return true;
			}
			curEdge = curElem.neighborEdge(uLVertex, curEdge);
			curElem = curElem.neighbor(curEdge);
		}

		// Parse all adjacent elements at upper right Vertex from neigbor of top edge
		// to,
		// but not including, neighbor of right edge.
		curElem = neighbor(edgeList[top]);
		curEdge = edgeList[top];
		while (curElem != null && curEdge != edgeList[right]) {
			if (curElem.inverted()) {
				return true;
			}
			curEdge = curElem.neighborEdge(uRVertex, curEdge);
			curElem = curElem.neighbor(curEdge);
		}

		// Parse all adjacent elements at lower right Vertex from neigbor of right edge
		// to,
		// but not including, neighbor of base edge.
		curElem = neighbor(edgeList[right]);
		curEdge = edgeList[right];
		while (curElem != null && curEdge != edgeList[base]) {
			if (curElem.inverted()) {
				return true;
			}
			curEdge = curElem.neighborEdge(edgeList[base].rightVertex, curEdge);
			curElem = curElem.neighbor(curEdge);
		}

		// Parse all adjacent elements at lower left Vertex from neigbor of base edge
		// to,
		// but not including, neighbor of left edge.
		curElem = neighbor(edgeList[base]);
		curEdge = edgeList[base];
		while (curElem != null && curEdge != edgeList[left]) {
			if (curElem.inverted()) {
				return true;
			}
			curEdge = curElem.neighborEdge(edgeList[base].leftVertex, curEdge);
			curElem = curElem.neighbor(curEdge);
		}
		return false;
	}

	/** @return a list of all triangles adjacent to this quad. */
	public List<Element> getAdjTriangles() {
		List<Element> triangleList;
		Vertex uLVertex = edgeList[left].otherVertex(edgeList[base].leftVertex);
		Vertex uRVertex = edgeList[right].otherVertex(edgeList[base].rightVertex);
		Vertex bLVertex = edgeList[base].leftVertex;
		Vertex bRVertex = edgeList[base].rightVertex;

		triangleList = bLVertex.adjTriangles();

		// Parse all adjacent elements at upper left Vertex from, but not not including,
		// neigbor of left edge to, but not including, neighbor of top edge.
		Element curElem = neighbor(edgeList[left]);
		Edge curEdge = edgeList[left];
		if (curElem != null) {
			curEdge = curElem.neighborEdge(uLVertex, curEdge);
			curElem = curElem.neighbor(curEdge);

			while (curElem != null && curEdge != edgeList[top]) {
				if (curElem instanceof Triangle) {
					triangleList.add(curElem);
				}
				curEdge = curElem.neighborEdge(uLVertex, curEdge);
				curElem = curElem.neighbor(curEdge);
			}
		}
		// Parse all adjacent elements at upper right Vertex from neigbor of top edge
		// to,
		// but not including, neighbor of right edge.
		curElem = neighbor(edgeList[top]);
		curEdge = edgeList[top];
		while (curElem != null && curEdge != edgeList[right]) {
			if (curElem instanceof Triangle && !triangleList.contains(curElem)) {
				triangleList.add(curElem);
			}
			curEdge = curElem.neighborEdge(uRVertex, curEdge);
			curElem = curElem.neighbor(curEdge);
		}

		triangleList.addAll(bRVertex.adjTriangles());

		return triangleList;
	}

	/** @return a list of all vertices adjacent to this quad. */
	public List<Vertex> getAdjVertices() {
		ArrayList<Vertex> vertexList = new ArrayList<>();
		Edge e;
		Vertex n;
		Vertex bLVertex = edgeList[base].leftVertex;
		Vertex bRVertex = edgeList[base].rightVertex;
		Vertex uLVertex = edgeList[left].otherVertex(bLVertex);

		int i;

		for (i = 0; i < bLVertex.edgeList.size(); i++) {
			e = bLVertex.edgeList.get(i);
			if (e != edgeList[base] && e != edgeList[left] && e != edgeList[right] && e != edgeList[top]) {
				vertexList.add(e.otherVertex(bLVertex));
			}
		}

		for (i = 0; i < bRVertex.edgeList.size(); i++) {
			e = bRVertex.edgeList.get(i);
			if (e != edgeList[base] && e != edgeList[left] && e != edgeList[right] && e != edgeList[top]) {
				n = e.otherVertex(bRVertex);
				if (!vertexList.contains(n)) {
					vertexList.add(n);
				}
			}
		}

		for (i = 0; i < uLVertex.edgeList.size(); i++) {
			e = uLVertex.edgeList.get(i);
			if (e != edgeList[base] && e != edgeList[left] && e != edgeList[right] && e != edgeList[top]) {
				n = e.otherVertex(uLVertex);
				if (!vertexList.contains(n)) {
					vertexList.add(n);
				}
			}
		}

		if (!isFake) {
			Vertex uRVertex = edgeList[right].otherVertex(bRVertex);
			for (i = 0; i < uRVertex.edgeList.size(); i++) {
				e = uRVertex.edgeList.get(i);
				if (e != edgeList[base] && e != edgeList[left] && e != edgeList[right] && e != edgeList[top]) {
					n = e.otherVertex(uRVertex);
					if (!vertexList.contains(n)) {
						vertexList.add(n);
					}
				}
			}
		}
		return vertexList;
	}

	@Override
	public void replaceEdge(Edge e, Edge replacement) {
		edgeList[indexOf(e)] = replacement;
	}

	/**
	 * Make the Element pointers in each of the Edges in this Quad point to this
	 * Quad.
	 */
	@Override
	public void connectEdges() {
		edgeList[base].connectToQuad(this);
		edgeList[left].connectToQuad(this);
		edgeList[right].connectToQuad(this);
		if (!isFake) {
			edgeList[top].connectToQuad(this);
		}
	}

	/**
	 * Release the element pointer of the edges in edgeList that pointed to this
	 * Quad.
	 */
	@Override
	public void disconnectEdges() {
		edgeList[base].disconnectFromElement(this);
		edgeList[left].disconnectFromElement(this);
		edgeList[right].disconnectFromElement(this);
		if (!isFake) {
			edgeList[top].disconnectFromElement(this);
		}
	}

	/**
	 * @return an Edge that is common to both this Quad and Quad q at Vertex n.
	 *         Return null if none exists.
	 */
	public Edge commonEdgeAt(Vertex n, Quad q) {
		for (Edge e : n.edgeList) {
			if (hasEdge(e) && q.hasEdge(e)) {
				return e;
			}
		}
		return null;
	}

	/**
	 * @param q a neighbor quad sharing an edge with this quad.
	 * @return an edge that is common to both this quad and quad q. Return null if
	 *         none exists.
	 */
	public Edge commonEdge(Quad q) {
		if (q == neighbor(edgeList[base])) {
			return edgeList[base];
		} else if (q == neighbor(edgeList[left])) {
			return edgeList[left];
		} else if (q == neighbor(edgeList[right])) {
			return edgeList[right];
		} else if (q == neighbor(edgeList[top])) {
			return edgeList[top];
		} else {
			return null;
		}
	}

	/**
	 * @return true if at least one of the edges connected to Vertex n is a front
	 *         edge.
	 */
	public boolean hasFrontEdgeAt(Vertex n) {
		if (edgeList[left].hasVertex(n)) {
			if (edgeList[base].hasVertex(n)) {
				if (edgeList[left].isFrontEdge() || edgeList[base].isFrontEdge()) {
					return true;
				} else {
					return false;
				}
			} else if (edgeList[top].hasVertex(n)) {
				if (edgeList[left].isFrontEdge() || edgeList[top].isFrontEdge()) {
					return true;
				} else {
					return false;
				}
			}
		} else if (edgeList[right].hasVertex(n)) {
			if (edgeList[base].hasVertex(n)) {
				if (edgeList[right].isFrontEdge() || edgeList[base].isFrontEdge()) {
					return true;
				} else {
					return false;
				}
			} else if (edgeList[top].hasVertex(n)) {
				if (edgeList[right].isFrontEdge() || edgeList[top].isFrontEdge()) {
					return true;
				} else {
					return false;
				}
			}
		}
		return false;
	}

	/**
	 * @return the number of quad neighbors sharing an edge with this quad at Vertex
	 *         n. This quad is not counted. Values are 0, 1, or 2.
	 */
	public int nrOfQuadsSharingAnEdgeAt(Vertex n) {
		int count = 0;

		if (edgeList[left].hasVertex(n)) {
			if (neighbor(edgeList[left]) instanceof Quad) {
				count++;
			}
			if (edgeList[base].hasVertex(n)) {
				if (neighbor(edgeList[base]) instanceof Quad) {
					count++;
				}
			} else if (neighbor(edgeList[top]) instanceof Quad) {
				count++;
			}
			return count;
		} else if (edgeList[right].hasVertex(n)) {
			if (neighbor(edgeList[right]) instanceof Quad) {
				count++;
			}
			if (edgeList[base].hasVertex(n)) {
				if (neighbor(edgeList[base]) instanceof Quad) {
					count++;
				}
			} else if (neighbor(edgeList[top]) instanceof Quad) {
				count++;
			}

			return count;
		}
		return count;
	}

	/**
	 * Update the distortion metric according to the paper "An approach to Combined
	 * Laplacian and Optimization-Based Smoothing for Triangular, Quadrilateral and
	 * Quad-Dominant Meshes" by by Cannan, Tristano, and Staten
	 * 
	 * @return negative values for inverted quadrilaterals, else positive.
	 *         Equilateral quadrilaterals should return the maximum value of 1.
	 */
	//
	// This is a simple sketch of the quadrilateral with vertices and divided
	// into four triangles:
	//
	// n3__________n4
	// |\ /|
	// | \ t4 / |
	// | \ / |
	// | t2 X t3 |
	// | / \ |
	// | / t1 \ |
	// |/___________\|
	// n1 n2
	//
	// Also, I tried to sketch the case where the quad has an angle > than 180
	// degrees
	// Note that t3 is part of t1 and that t4 is part of t2 in the sketch.
	// t3 and t4 are inverted.
	//
	// n3
	// |\ \
	// | \ \
	// | \t4 X
	// |t2 \ / \
	// | /\ t3 \
	// | /t1 ---__n2
	// |/______-----=
	// n1
	@Override
	public void updateDistortionMetric() {
		if (isFake) {
			double AB = edgeList[base].len, CB = edgeList[left].len, CA = edgeList[right].len;

			Vertex a = edgeList[base].commonVertex(edgeList[right]), b = edgeList[base].commonVertex(edgeList[left]),
					c = edgeList[left].commonVertex(edgeList[right]);
			MyVector vCA = new MyVector(c, a), vCB = new MyVector(c, b);

			double temp = sqrt3x2 * Math.abs(vCA.cross(vCB)) / (CA * CA + AB * AB + CB * CB);
			if (inverted()) {
				distortionMetric = -temp;
			} else {
				distortionMetric = temp;
			}

			return;
		}

		Vertex n1 = edgeList[base].leftVertex;
		Vertex n2 = edgeList[base].rightVertex;
		Vertex n3 = edgeList[left].otherVertex(n1);
		Vertex n4 = edgeList[right].otherVertex(n2);

		// The two diagonals
		Edge e1 = new Edge(n1, n4);
		Edge e2 = new Edge(n2, n3);

		// The four triangles
		Triangle t1 = new Triangle(edgeList[base], edgeList[left], e2);
		Triangle t2 = new Triangle(edgeList[base], e1, edgeList[right]);
		Triangle t3 = new Triangle(edgeList[top], edgeList[right], e2);
		Triangle t4 = new Triangle(edgeList[top], e1, edgeList[left]);

		// Place the firstVertices correctly
		t1.firstVertex = firstVertex;
		if (firstVertex == n1) {
			t2.firstVertex = n1;
			t3.firstVertex = n4;
			t4.firstVertex = n4;
		} else {
			t2.firstVertex = n2;
			t3.firstVertex = n3;
			t4.firstVertex = n3;
		}

		// Compute and get alpha values for each triangle
		t1.updateDistortionMetric(4.0);
		t2.updateDistortionMetric(4.0);
		t3.updateDistortionMetric(4.0);
		t4.updateDistortionMetric(4.0);

		double alpha1 = t1.distortionMetric, alpha2 = t2.distortionMetric, alpha3 = t3.distortionMetric, alpha4 = t4.distortionMetric;

		int invCount = 0;
		if (alpha1 < 0) {
			invCount++;
		}
		if (alpha2 < 0) {
			invCount++;
		}
		if (alpha3 < 0) {
			invCount++;
		}
		if (alpha4 < 0) {
			invCount++;
		}

		double temp12 = Math.min(alpha1, alpha2);
		double temp34 = Math.min(alpha3, alpha4);
		double alphaMin = Math.min(temp12, temp34);
		double negval = 0;

		if (invCount >= 3) {
			if (invCount == 3) {
				negval = 2.0;
			} else {
				negval = 3.0;
			}
		} else if (ang[0] < DEG_6 || ang[1] < DEG_6 || ang[2] < DEG_6 || ang[3] < DEG_6 || coincidentVertices(n1, n2, n3, n4)
				|| invCount == 2) {
			negval = 1.0;
		}

		distortionMetric = alphaMin - negval;
	}

	/** Test whether any vertices of the quad are coincident. */
	private boolean coincidentVertices(Vertex n1, Vertex n2, Vertex n3, Vertex n4) {
		double x12diff = n2.x - n1.x;
		double y12diff = n2.y - n1.y;
		double x13diff = n3.x - n1.x;
		double y13diff = n3.y - n1.y;
		double x14diff = n4.x - n1.x;
		double y14diff = n4.y - n1.y;

		double x23diff = n3.x - n2.x;
		double y23diff = n3.y - n2.y;
		double x24diff = n4.x - n2.x;
		double y24diff = n4.y - n2.y;

		double x34diff = n4.x - n3.x;
		double y34diff = n4.y - n3.y;

		// Using Pythagoras: hyp^2= kat1^2 + kat2^2
		double l12 = Math.sqrt(x12diff * x12diff + y12diff * y12diff);
		double l13 = Math.sqrt(x13diff * x13diff + y13diff * y13diff);
		double l14 = Math.sqrt(x14diff * x14diff + y14diff * y14diff);
		double l23 = Math.sqrt(x23diff * x23diff + y23diff * y23diff);
		double l24 = Math.sqrt(x24diff * x24diff + y24diff * y24diff);
		double l34 = Math.sqrt(x34diff * x34diff + y34diff * y34diff);

		if (l12 < COINCTOL || l13 < COINCTOL || l14 < COINCTOL || l23 < COINCTOL || l24 < COINCTOL || l34 < COINCTOL) {
			return true;
		} else {
			return false;
		}
	}

	/** @return the size of the largest interior angle */
	@Override
	public double largestAngle() {
		double cand = ang[0];
		if (ang[1] > cand) {
			cand = ang[1];
		}
		if (ang[2] > cand) {
			cand = ang[2];
		}
		if (ang[3] > cand) {
			cand = ang[3];
		}
		return cand;
	}

	/** @return the Vertex at the largest interior angle */
	@Override
	public Vertex VertexAtLargestAngle() {
		Vertex candVertex = edgeList[base].leftVertex;
		double cand = ang[0];

		if (ang[1] > cand) {
			candVertex = edgeList[base].rightVertex;
			cand = ang[1];
		}
		if (ang[2] > cand) {
			candVertex = edgeList[left].otherVertex(edgeList[base].leftVertex);
			cand = ang[2];
		}
		if (ang[3] > cand) {
			candVertex = edgeList[right].otherVertex(edgeList[base].rightVertex);
		}
		return candVertex;
	}

	/** @return the length of the longest Edge in the quad */
	@Override
	public double longestEdgeLength() {
		double t1 = Math.max(edgeList[base].len, edgeList[left].len);
		double t2 = Math.max(t1, edgeList[right].len);
		return Math.max(t2, edgeList[top].len);
	}

	/**
	 * @param first a triangle that is located inside the quad
	 * @return a list of triangles contained within the four edges of this quad.
	 */
	public ArrayList<Triangle> trianglesContained(Triangle first) {
		ArrayList<Triangle> tris = new ArrayList<>();
		Element neighbor;
		Triangle cur;
		Edge e;

		tris.add(first);
		for (int j = 0; j < tris.size(); j++) {
			cur = tris.get(j);
			for (int i = 0; i < 3; i++) {
				e = cur.edgeList[i];
				if (!hasEdge(e)) {
					neighbor = cur.neighbor(e);
					if (neighbor != null && !tris.contains(neighbor)) {
						tris.add((Triangle) neighbor);
					}
				}
			}
		}
		return tris;
	}

	/**
	 * Test whether the quad contains a hole.
	 * 
	 * @param tris the interior triangles
	 * @return true if there are one or more holes present within the four edges
	 *         defining the quad.
	 */
	public boolean containsHole(List<Triangle> tris) {
		Triangle t;

		if (tris.size() == 0) {
			return true; // Corresponds to a quad defined on a quad-shaped hole
		}

		for (Object element : tris) {
			t = (Triangle) element;
			if (t.edgeList[0].boundaryEdge() && !hasEdge(t.edgeList[0]) || t.edgeList[1].boundaryEdge() && !hasEdge(t.edgeList[1])
					|| t.edgeList[2].boundaryEdge() && !hasEdge(t.edgeList[2])) {
				return true;
			}
		}
		return false;
	}

	/** Set the color of the edges to green. */
	@Override
	public void markEdgesLegal() {
		edgeList[base].color = java.awt.Color.green;
		edgeList[left].color = java.awt.Color.green;
		edgeList[right].color = java.awt.Color.green;
		edgeList[top].color = java.awt.Color.green;
	}

	/** Set the color of the edges to red. */
	@Override
	public void markEdgesIllegal() {
		edgeList[base].color = java.awt.Color.red;
		edgeList[left].color = java.awt.Color.red;
		edgeList[right].color = java.awt.Color.red;
		edgeList[top].color = java.awt.Color.red;
	}

	/**
	 * Give a string representation of the quad.
	 * 
	 * @return a string representation of the quad.
	 */
	@Override
	public String descr() {
		Vertex Vertex1, Vertex2, Vertex3, Vertex4;
		Vertex1 = edgeList[base].leftVertex;
		Vertex2 = edgeList[base].rightVertex;
		Vertex3 = edgeList[left].otherVertex(Vertex1);
		Vertex4 = edgeList[right].otherVertex(Vertex2);
		String v3Desc = Vertex3 == null ? "" : Vertex3.descr();
		String v4Desc = Vertex4 == null ? "" : Vertex4.descr();
		return Vertex1.descr() + ", " + Vertex2.descr() + ", " + v3Desc + ", " + v4Desc;
	}

	/** Output a string representation of the quad. */
	@Override
	public void printMe() {
		System.out.println(descr() + ", inverted(): " + inverted() + ", ang[0]: " + Math.toDegrees(ang[0]) + ", ang[1]: "
				+ Math.toDegrees(ang[1]) + ", ang[2]: " + Math.toDegrees(ang[2]) + ", ang[3]: " + Math.toDegrees(ang[3])
				+ ", firstVertex is " + firstVertex.descr());
	}
}
