package meshditor;

/** A very simple implementation of a dart. */
public class Dart {
	
	public Vertex n;
	public Edge e;
	public Element elem;
	
	public Dart() {
		n = null;
		e = null;
		elem = null;
	}

	public Dart(Vertex n, Edge e, Element elem) {
		this.n = n;
		this.e = e;
		this.elem = elem;
	}

	public String descr() {
		return "(elem: " + elem.descr() + ", e: " + e.descr() + ", n: " + n.descr() + ")";
	}

}
