## Esta función requiere un universo creado con MDAnalysis (trr y tpr file)
## y especificar si quieres la monocapa superior (up) o inferior (down)
## Para usar con distintos lípidos cambiar sel1,2,3, dependiendo de la situación.


def up_or_down(universe, side):
    """Función para separar a la bicapa en dos monocapas.
       c/u de ellas como un AtomGroup"""

    sel1 = "resname POPC and name C32"
    sel2 = "resname PSM and name C2F"
    sel3 = "resname CHL1 and name C23 "
    u.trajectory[0] 
    # Select atoms within this particular frame
    lpd1_atoms = u.select_atoms(sel1) 
    #print lpd1_atoms, u.select_atoms
    lpd2_atoms = u.select_atoms(sel2) 
    lpd3_atoms = u.select_atoms(sel3) 

    zmean = np.mean(np.concatenate([lpd1_atoms.positions[:,2], lpd2_atoms.positions[:,2], lpd3_atoms.positions[:,2]]))
    #print zmean

    num_lpd1 = lpd1_atoms.n_atoms # Número de átomos de POPC
    num_lpd2 = lpd2_atoms.n_atoms # Número de átomos de PSM 
    num_lpd3 = lpd3_atoms.n_atoms # Número de átomos de CHO
    
    if side == "up":
        lpd1i = []
        zpos = lpd1_atoms.positions[:,2] # Matriz con las posiciones de lpd1_atoms [:,2]==(z coord) 
        for i in range(num_lpd1): 
            if zpos[i] > zmean: # seleccionando valores por encima de zmean
                lpd1i.append(lpd1_atoms[i]) # Se agrega a la matriz lpd1i, todos los valores que indican que estan en la bicapa superior
        lpd1i = np.sum(np.array(lpd1i)) # AtomGroup con todos los átomos del lípido seleccionado en up
        print(type(lpd1i))

        lpd2i = []
        zpos = lpd2_atoms.positions[:,2]
        for i in range(num_lpd2):
            if zpos[i] > zmean:
                lpd2i.append(lpd2_atoms[i])
        lpd2i = np.sum(np.array(lpd2i)) 
        print(type(lpd2i))
        lpd3i = []
        zpos = lpd3_atoms.positions[:,2]
        for i in range(num_lpd3):
            if zpos[i] > zmean:
                lpd3i.append(lpd3_atoms[i])
        lpd3i = np.sum(np.array(lpd3i)) 
        print(type(lpd3i))

    elif side == 'down':
        lpd1i = []
        zpos = lpd1_atoms.positions[:,2]
        for i in range(num_lpd1):
            if zpos[i] < zmean:
                lpd1i.append(lpd1_atoms[i])
        lpd1i = np.sum(np.array(lpd1i))

        lpd2i = []
        zpos = lpd2_atoms.positions[:,2]
        for i in range(num_lpd2):
            if zpos[i] < zmean:
                lpd2i.append(lpd2_atoms[i])
        lpd2i = np.sum(np.array(lpd2i))

        lpd3i = []
        zpos = lpd3_atoms.positions[:,2]
        for i in range(num_lpd3):
            if zpos[i] > zmean:
                lpd3i.append(lpd3_atoms[i])
        lpd3i = np.sum(np.array(lpd3i)) 

    lpd_atms = lpd1i + lpd2i + lpd3i


    return lpd_atms

monocapa = up_or_down(u, "up")
