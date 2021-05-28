import numpy as np
import matplotlib.pyplot as plt
import gmsh
np.set_printoptions(precision = 4, linewidth = 150)


gmsh.initialize()
gmsh.model.add('test')

lc = 20
L = 10

#Primero puntos
p1 = gmsh.model.geo.addPoint(0, 0, 0, lc)
p2 = gmsh.model.geo.addPoint(2*L, 0, 0, lc)
p3 = gmsh.model.geo.addPoint(2*L, L, 0, lc)
p4 = gmsh.model.geo.addPoint(0, L, 0, lc)

#Luego líneas
l1 = gmsh.model.geo.addLine(p1,p2)
l2 = gmsh.model.geo.addLine(p2,p3)
l3 = gmsh.model.geo.addLine(p3,p4)
l4 = gmsh.model.geo.addLine(p4,p1)

#Luego bordes
C1 = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
S1 = gmsh.model.geo.addPlaneSurface([C1])

gmsh.model.geo.synchronize()
#Los physical groups son categorías, tengo que darle la
#dimensión a cada grupo y por qué está compuesto
Empotrado = gmsh.model.addPhysicalGroup(1, [l4])
gmsh.model.setPhysicalName(1, Empotrado, 'Empotrado')
Traccionado = gmsh.model.addPhysicalGroup(1, [l2])
gmsh.model.setPhysicalName(1, Traccionado, 'Traccionado')
Superficie = gmsh.model.addPhysicalGroup(2, [S1])
gmsh.model.setPhysicalName(2, Superficie, 'Superficie')

gmsh.model.geo.synchronize()

#Genero el mallado, y quiero que sea 2D
gmsh.model.mesh.generate(2)
#Veo como va quedando la cosa
gmsh.model.geo.synchronize()
#gmsh.fltk.run()

#Obtengo info
NodeInfo = gmsh.model.mesh.get_nodes()
NumeroNodos = NodeInfo[0].shape[0]

MatrizNodos = NodeInfo[1].reshape(NumeroNodos, 3)

#plt.rc('figure', figsize=(10,10))
#plt.plot(MatrizNodos[:,0], MatrizNodos[:,1], 'ok')
#plt.show()

ElementInfo = gmsh.model.mesh.get_elements()
ETYPES = ElementInfo[0]

ETAGS, ELEMENTS = gmsh.model.mesh.get_elements_by_type(2)
MatrizConectividad = ELEMENTS.reshape([ETAGS.shape[0],3])

NodosEmpotrados = gmsh.model.mesh.get_nodes_for_physical_group(1,Empotrado)
NodosTraccionados = gmsh.model.mesh.get_nodes_for_physical_group(1, Traccionado)

print(NodosTraccionados)
print(NodosEmpotrados)

np.savetxt('MatrizNodos2D.dat', MatrizNodos, fmt='%1.3f')
np.savetxt('MatrizConectividad2D.dat', MatrizConectividad, fmt = '%d')
np.savetxt('TraccionadosX.dat', NodosTraccionados[0], fmt = '%d')
np.savetxt('Empotrados.dat', NodosEmpotrados[0], fmt = '%d')

gmsh.model.geo.synchronize()
gmsh.fltk.run()
gmsh.finalize()