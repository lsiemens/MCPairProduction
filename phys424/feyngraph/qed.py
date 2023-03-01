"""Generate QED feynman diagrams
"""

class QEDError(Exception):
    pass

class GraphError(QEDError):
    pass

class VertexError(QEDError):
    pass

class ExternalParticleError(QEDError):
    pass


def get_id(obj):
    """Helper function for printing the unique index of an object
    """
    if obj is not None:
        return obj.id
    else:
        return "None"


class ExternalParticle:
    _ID = 0

    def __init__(self):
        self._id = ExternalParticle._ID
        ExternalParticle._ID += 1

        self.id = "ExternalParticle"

    def validate(self):
        raise NotImplementedError("ExternalParticle is class prototype")

    def get_connection(self):
        raise NotImplementedError("ExternalParticle is class prototype")

    def set_connection(self, connection):
        raise NotImplementedError("ExternalParticle is class prototype")

    def __str__(self):
        string = f"ExternalParticle(id:{get_id(self)}, connection:{get_id(self.get_connection())})"
        return string

class ExternalFin(ExternalParticle):
    def __init__(self):
        super().__init__()
        self.fin = None

        self.id = "ExternalFin"

    def get_connection(self):
        return self.fin

    def set_connection(self, connection):
        self.fin = connection

    def validate(self):
        if isinstance(self.fin, Vertex) or isinstance(self.fin, ExternalFout):
            if (self.fin).fout != self:
                raise ExternalParticleError("External incomming particle is not connected to outgoing particle or external outgoing particle.")
        if isinstance(self.fin, ExternalPhoton) or isinstance(self.fin, ExternalFin):
            raise ExternalParticleError("External incomming particle is not connected to outgoing particle or external outgoing particle.")

class ExternalFout(ExternalParticle):
    def __init__(self):
        super().__init__()
        self.fout = None

        self.id = "ExternalFout"

    def get_connection(self):
        return self.fout

    def set_connection(self, connection):
        self.fout = connection

    def validate(self):
        if isinstance(self.fout, Vertex) or isinstance(self.fout, ExternalFin):
            if (self.fout).fin != self:
                raise ExternalParticleError("External outgoing particle is not connected to incomming particle or external incomming particle.")
        if isinstance(self.fout, ExternalPhoton) or isinstance(self.fout, ExternalFout):
            raise ExternalParticleError("External outgoing particle is not connected to incomming particle or external incomming particle.")

class ExternalPhoton(ExternalParticle):
    def __init__(self):
        super().__init__()
        self.photon = None

        self.id = "ExternalPhoton"

    def get_connection(self):
        return self.photon

    def set_connection(self, connection):
        self.photon = connection

    def validate(self):
        if isinstance(self.photon, Vertex) or isinstance(self.photon, ExternalPhoton):
            if (self.photon).photon != self:
                raise ExternalParticleError("External photon is not connected to photon or external photon.")
        if isinstance(self.photon, ExternalFin) or isinstance(self.photon, ExternalFout):
            raise ExternalParticleError("External photon is not connected to photon or external photon.")

class InitalFermion(ExternalFout):
    def __init__(self):
        super().__init__()

        self.id = f"iF{self._id}"

class InitalAFermion(ExternalFin):
    def __init__(self):
        super().__init__()

        self.id = f"iaF{self._id}"

class InitalPhoton(ExternalPhoton):
    def __init__(self):
        super().__init__()

        self.id = f"iP{self._id}"
    pass

class FinalFermion(ExternalFin):
    def __init__(self):
        super().__init__()

        self.id = f"oF{self._id}"

class FinalAFermion(ExternalFout):
    def __init__(self):
        super().__init__()

        self.id = f"oaF{self._id}"

class FinalPhoton(ExternalPhoton):
    def __init__(self):
        super().__init__()

        self.id = f"oP{self._id}"


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
        string = "InternalGraph\n\t Vertices[" + ",".join(["\n\t\t" + str(vertex) for vertex in self.vertices]) + "]"
        return string


class FeynmanDiagram:
    def __init__(self, vertices, externalParticles):
        self.vertices = vertices
        self.externalParticles = externalParticles

    def validate(self):
        """Validate the Feynman Diagram, all vertices and propegators
        """

        for vertex in self.vertices:
            vertex.validate()

            if (vertex.fin not in self.vertices) and (vertex.fin not in self.externalParticles):
                raise GraphError("Incoming internal fermion not in Feynman diagram.")

            if (vertex.fout not in self.vertices) and (vertex.fout not in self.externalParticles):
                raise GraphError("Outgoing internal fermion not in Feynman diagram.")

            if (vertex.photon not in self.vertices) and (vertex.photon not in self.externalParticles):
                raise GraphError("Internal photon not in Feynman diagram.")

        for externalParticle in self.externalParticles:
            externalParticle.validate()

            if (externalParticle.get_connection() not in self.vertices) and (externalParticle.get_connection() not in self.externalParticles):
                raise GraphError("External particle not in Feynman diagram.")

    def __str__(self):
        string = "FeynmanDiagram\n\t Vertices[" + ",".join(["\n\t\t" + str(vertex) for vertex in self.vertices]) + "]\n\t ExternalParticles[" + ",".join(["\n\t\t" + str(externalParticle) for externalParticle in self.externalParticles]) + "]"
        return string


# define internal graph

ein = ExternalParticle()
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
graph.validate()
print(graph)

# define Feynman diagram

i1 = InitalPhoton()
i2 = InitalFermion()
o1 = FinalPhoton()
o2 = FinalFermion()

v1 = Vertex()
v2 = Vertex()

v1.photon = i1
i1.set_connection(v1)

v1.fin = i2
i2.set_connection(v1)

v1.fout = v2
v2.fin = v1

v2.photon = o1
o1.set_connection(v2)

v2.fout = o2
o2.set_connection(v2)

diagram = FeynmanDiagram([v1, v2], [i1, i2, o1, o2])
diagram.validate()
print(diagram)
