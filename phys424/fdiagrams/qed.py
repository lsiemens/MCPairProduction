"""Generate QED feynman diagrams
"""

class QEDError(Exception):
    pass

class GraphError(QEDError):
    pass

class VertexError(QEDError):
    pass

def get_id(obj):
    if obj is not None:
        return obj.id
    else:
        return "None"

class ExternalParticle:
    _ID = 0

    def __init__(self, incoming):
        self._id = ExternalParticle._ID
        ExternalParticle._ID += 1

        self.incoming = incoming
        if self.incoming:
            self.id = f"i{self._id}"
        else:
            self.id = f"o{self._id}"

class Vertex:
    """QED Vertex

    Define a vertex in a Feynman diagram. Note a incoming fermion (fin)
    must link to the outgoin fermion (fout) of another vertex, and photons
    must always link together.

    Which is to say for a pair of vertices `self` and `other` if `self.fin`
    points to the vertex `other`, then it must point back with `other.fout`.
    """
    _ID = 0

    def __init__(self):
        self._id = Vertex._ID
        Vertex._ID += 1

        self.id = f"v{self._id}"
        self.fin = None
        self.fout = None
        self.photon = None

    def validate(self):
        """Check if fin connects to fout, photon to photon ...
        """
        if isinstance(self.photon, Vertex):
            if (self.photon).photon != self:
                raise VertexError("Photon not connected to photon.")

        if isinstance(self.fin, Vertex):
            if (self.fin).fout != self:
                raise VertexError("Incomming fermion not connected to outgoing fermio.")

        if isinstance(self.fout, Vertex):
            if (self.fout).fin != self:
                raise VertexError("Outgoing fermion not connected to incomming fermio.")

    def __str__(self):
        string = f"Vertex(id:{self._id}, fin:{get_id(self.fin)}, fout:{get_id(self.fout)}, photon:{get_id(self.photon)})"
        return string

class InternalGraph:
    def __init__(self, vertices):
        self.vertices = vertices

    def validate(self):
        """Validate all interconnections of the vertices
        """

        for vertex in self.vertices:
            vertex.validate()

            if isinstance(vertex.fin, Vertex):
                if vertex.fin not in self.vertices:
                    raise GraphError("Incoming internal fermion not in graph.")

            if isinstance(vertex.fout, Vertex):
                if vertex.fout not in self.vertices:
                    raise GraphError("Outgoing internal fermion not in graph.")

            if isinstance(vertex.photon, Vertex):
                if vertex.photon not in self.vertices:
                    raise GraphError("Internal photon not in graph.")

    def get_free_edges(self):
        """Get lists of fin, fout and photon connections that are not
        internal to the graph.
        """
        photons = []
        fins = []
        fouts = []

        self.validate()

        for vertex in self.vertices:
            if not isinstance(vertex.fin, Vertex):
                fins.append(vertex)

            if not isinstance(vertex.fout, Vertex):
                fouts.append(vertex)

            if not isinstance(vertex.photon, Vertex):
                photons.append(vertex)

        return fins, fouts, photons

    def __str__(self):
        string = "Vertices[" + ",".join([str(vertex) for vertex in self.vertices]) + "]"
        return string

ein = ExternalParticle(True)
v1 = Vertex()
v2 = Vertex()
v1.fin = ein
v1.fout = v2
v1.photon = v2
v2.fin = v1
v2.photon = v1

v3 = Vertex()
g2 = InternalGraph([v3])

graph = InternalGraph([v1, v2])

print(graph.get_free_edges())
