package meshditor;

/**
 * This class declares methods and variables that are common to quads and
 * triangles.
 */

public abstract class Element extends Constants {
	
	/** An array of interior angles */
	public double[] ang;
	/** An array of edges */
	public Edge[] edgeList;
	/** Vertex used for determining inversion, amonst other things. */
	Vertex firstVertex;
	
	/** @return neighbor element sharing edge e */
	public abstract Element neighbor(Edge e);

	/** @return local angle inside element at Vertex n */
	public abstract double angle(Edge e, Vertex n);

	/** Compute & set the angles at the vertices of the element. */
	public abstract void updateAngles();

	/**
	 * Compute & set the angle at this particular Vertex incident with this Element
	 * Edge
	 */
	public abstract void updateAngle(Vertex n);

	/** @return description string for element (list of Vertex coordinates) */
	public abstract String descr();

	/** Output description string for element (list of Vertex coordinates) */
	public abstract void printMe();

	/** Verify that the quad has the specified edge. */
	public abstract boolean hasEdge(Edge e);

	/** Verify that the quad has the specified Vertex. */
	public abstract boolean hasVertex(Vertex n);

	/** Verify that the area of the quad is greater than 0. */
	public abstract boolean areaLargerThan0();

	/** Return local neighboring edge at Vertex n. */
	public abstract Edge neighborEdge(Vertex n, Edge e);

	/** Return the index to this edge in this element's edgeList */
	public abstract int indexOf(Edge e);

	/** Return the index to this angle in this element's ang array */
	public abstract int angleIndex(Edge e1, Edge e2);

	public abstract int angleIndex(Vertex n);

	/** Return the angle between this Element's Edges e1 and e2. */
	public abstract double angle(Edge e1, Edge e2);

	/** @return true if the element has become inverted */
	public abstract boolean inverted();

	/** @return true if the element has become inverted or its area is zero. */
	public abstract boolean invertedOrZeroArea();

	/** @return true if the element has a concavity at its Vertex n. */
	public abstract boolean concavityAt(Vertex n);

	/** Replace one of the specified edges e with a replacement edge. */
	public abstract void replaceEdge(Edge e, Edge replacement);

	/** Make one element pointer of each Edge in edgeList point to this Element. */
	public abstract void connectEdges();

	/**
	 * Point the element pointer of each Edge in edgeList that previously pointed to
	 * this Element to point to null.
	 */
	public abstract void disconnectEdges();

	/** Create a simple element for testing purposes only. */
	public abstract Element elementWithExchangedVertices(Vertex original, Vertex replacement);

	/**
	 * @return true if the quad becomes inverted when Vertex n1 is relocated to pos.
	 *         n2. Else return false.
	 */
	public abstract boolean invertedWhenVertexRelocated(Vertex n1, Vertex n2);

	/**
	 * Update the distortion metric according to the article "An approach to
	 * Combined Laplacian and Optimization-Based Smoothing for Triangular,
	 * Quadrilateral and Quad-Dominant Meshes" by by Cannan, Tristano, and Staten
	 */
	public abstract void updateDistortionMetric();

	/**
	 * Doubles to hold the cur. distortion metric & the metric after perturbation
	 */
	double distortionMetric, newDistortionMetric;
	/** Doubles to hold the gradient vector */
	double gX, gY;

	/** Return the length of the longest Edge. */
	public abstract double longestEdgeLength();

	/** Return the size of the largest angle. */
	public abstract double largestAngle();

	/** Return the Vertex at the largest interior angle. */
	public abstract Vertex VertexAtLargestAngle();

	/** Set the color of the edges to red. */
	public abstract void markEdgesIllegal();

	/** Set the color of the edges to green. */
	public abstract void markEdgesLegal();

	/**
	 * A method for fast computation of the cross product of two vectors.
	 * 
	 * @param o1 origin of first vector
	 * @param p1 endpoint of first vector
	 * @param o2 origin of second vector
	 * @param p2 endpoint of second vector
	 * @return the cross product of the two vectors
	 */
	protected double cross(Vertex o1, Vertex p1, Vertex o2, Vertex p2) {
		double x1 = p1.x - o1.x;
		double x2 = p2.x - o2.x;
		double y1 = p1.y - o1.y;
		double y2 = p2.y - o2.y;
		return x1 * y2 - x2 * y1;
	}
}
