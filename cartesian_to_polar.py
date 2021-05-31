######### Programa para acomodar data (trayectoría DM) tal que extraemos las posiciones (x,y) de los átomos que queramos utilizar. 
######### Construimos un array con posiciones R_i(t^k) = x_i(t^k),y_i(t^k) de los i-átomos seleccionados a los distintos tiempos t^k 
######### Con i = 1,2,...N y k = 1,...,n respectivamente número de moléculas y número de frames en la simulación. 
######### Después se divide cada elemento y_i(t^k)/x_i(t^k) = arg, para usar ángulo polar = tan^-1(arg).
######## SE NECESITA LA FUNCIÓN up_or_down PARA LAS SELECCIONES INICIALES.

# Selecciones
atom1 = monocapa.select_atoms("name C32")
atom2 = monocapa.select_atoms("name C2F")
atom3 = monocapa.select_atoms("name C23")

lipid_name1 = atom1.resnames[0]
lipid_name2 = atom2.resnames[0]
lipid_name3 = atom3.resnames[0]

num_atom1 = atom1[:30]
num_atom2 = atom2[:30]
num_atom3 = atom3[:30]

# Posiciones x,y de c/u de las selecciones evolucionando en el tiempo
atom1_pos_t = np.array([num_atom1.positions[:, :2] for ts in u.trajectory[::100]]).reshape(51, 60)
atom2_pos_t = np.array([num_atom2.positions[:, :2] for ts in u.trajectory[::100]]).reshape(51, 60)
#atom3_pos_t = np.array([num_atom3.positions[:, :2] for ts in u.trajectory[::100]]).reshape(51,60)


atom1_df = pd.DataFrame(atom1_pos_t)
print(atom1_df)
atom2_df = pd.DataFrame(atom2_pos_t)
#print(atom2_df)

def cartesian_to_polars(dataframe, lpd_name):
    """Función para transformar un dataframe.shape = (ts, 2N) con (x, y) coords para cada átomo 
    a theta = tan^-1(y/x) con theta_df.shape = (ts, N)"""

    dimensions = dataframe.shape
    num_rows = dimensions[0]
    num_cols = dimensions[1]
    print("Las dimensiones del dataframe son (%d,%d)" %(num_rows, num_cols))

    polar_angle_arr = np.zeros(shape = (num_rows, num_cols//2))
    if lpd_name == lipid_name1:
        for col in range(num_cols - 30):
            for row in range(num_rows):
                #y_over_x_arr[row, col] = dataframe[(2*col) + 1][row] / dataframe[2*col][row]
                #y_over_x_df = pd.DataFrame(y_over_x_arr)
                y = dataframe[(2*col) + 1][row]
                x = dataframe[2*col][row]
                polar_angle_arr[row,col] = np.arctan2(y, x)
                polar_angle_df = pd.DataFrame(polar_angle_arr).add_prefix(lpd_name + "theta_")
                mean_angle = polar_angle_df.mean(axis = 0)
                std_polar_angle = polar_angle_df - mean_angle
                #print(polar_angle_df)
                #print(mean_angle[0])
                #print(std_polar_angle)
    elif lpd_name == lipid_name2:
        for col in range(num_cols - 30):
            for row in range(num_rows):
                y = dataframe[(2*col) + 1][row]
                x = dataframe[2*col][row]
                polar_angle_arr[row,col] = np.arctan2(y, x)
                polar_angle_df = pd.DataFrame(polar_angle_arr).add_prefix(lpd_name + "theta_")
                mean_angle = polar_angle_df.mean(axis = 0)
                std_polar_angle = polar_angle_df - mean_angle
        #print(std_polar_angle)

    else:
        lpd_name = "CHO"
        for col in range(num_cols - 3):
            for row in range(num_rows):
                y = dataframe[(2*col) + 1][row]
                x = dataframe[2*col][row]
                polar_angle_arr[row,col] = np.arctan2(y, x)
                polar_angle_df = pd.DataFrame(polar_angle_arr).add_prefix(lpd_name + "theta_")
                mean_angle = polar_angle_df.mean(axis = 0)
                std_polar_angle = polar_angle_df - mean_angle
        #print(std_polar_angle)
        
    return std_polar_angle

popc = cartesian_to_polars(atom1_df, "POPC")
#print(popc)
psm = cartesian_to_polars(atom2_df, "PSM")

dataframe_final = pd.concat([popc, psm], axis = 1)
print(dataframe_final)
dataframe_final.to_csv("atan2_polarstd_angle_PSMPOPCup_cho35.csv")
