def acomodar_data(monocapa):

    """ Función para acomodar data (trayectoría DM) tal que extraemos las posiciones de los átomos que queramos utilizar. 
    Construimos un array con posiciones R_i(t^k) = x_i(t^k),y_i(t^k),z_i(t^k) de los i-átomos seleccionados a los distintos tiempos t^k 
    Con i = 1,2,...N y k = 1,...,n respectivamente número de moléculas y número de frames en la simulación. 
    Después se calcula la media (m_i) de cada columna y se le resta a c/u de los elementos de la columna que le corresponde
    Sacamos el cuadrado a todos los elementos y sumamos deltaR_i(t^k)= [x_i(t^k)]^2 + [y_i(t^k)]^2 + [z_i(t^k)]^2 
    Y obtenemos raíz cuadrada de todos los elementos de la matriz contraída por la suma, tal que obtenemos la magnitud.
    El valor del parámetro monocapa proviene de la función up_or_down"""

    seleccion1 = monocapa.select_atoms("name C32")
    seleccion2 = monocapa.select_atoms("name C2F")
    seleccion3 = monocapa.select_atoms("name C23")

    num_sel1 = seleccion1[:5] # Seleccionar el número de átomos pertenecientes a la selección1,2o3
    num_sel2 = seleccion2[:5]
    num_sel3 = seleccion3[:5]

    sel1_pos_t = np.array([num_sel1.positions for ts in u.trajectory[::10]]).reshape(501,15)
    sel2_pos_t = np.array([num_sel2.positions for ts in u.trajectory[::10]]).reshape(501,15)
    sel3_pos_t = np.array([num_sel3.positions for ts in u.trajectory[::10]]).reshape(501,15)

    sel1_mean_t = sel1_pos_t.mean(axis = 0).reshape(1, 15) # La media de todas las columnas y se cambia de forma de (30,) a (1,30)
    sel2_mean_t = sel2_pos_t.mean(axis = 0).reshape(1, 15)
    sel3_mean_t = sel3_pos_t.mean(axis = 0).reshape(1, 15)
    #print(sel1_mean_t)
    #print(sel1_pos_t)
    sel1_minus_mean = sel1_pos_t - sel1_mean_t  # Se resta la media de esa columna a cada uno de los frames y así en todas las columas (broadcast)
    sel2_minus_mean = sel2_pos_t - sel2_mean_t
    sel3_minus_mean = sel3_pos_t - sel3_mean_t

    sel1_df = pd.DataFrame(sel1_minus_mean)
    sel1_df2 = sel1_df**2
    sel2_df = pd.DataFrame(sel2_minus_mean)
    sel2_df2 = sel2_df**2
    sel3_df = pd.DataFrame(sel3_minus_mean)
    sel3_df2 = sel3_df**2
    # La siguiente linea suma 3 columnas renglón x renglón. Se genera una nueva columna que se llena con los resultados
    deltaR_tk_sel1 = sel1_df2.groupby((np.arange(len(sel1_df2.columns)) // 3) + 1, axis=1).sum().add_prefix('POPC_') 
    #print(deltaR_tk)
    squared_deltaR_tk_sel1 = deltaR_tk_sel1.transform(lambda x : x**0.5) #Raíz cuadrada a todos los elementos.
    deltaR_tk_sel2 = sel2_df2.groupby((np.arange(len(sel2_df2.columns)) // 3) + 1, axis=1).sum().add_prefix('PSM_')
    squared_deltaR_tk_sel2 = deltaR_tk_sel2.transform(lambda x : x**0.5)
    deltaR_tk_sel3 = sel3_df2.groupby((np.arange(len(sel3_df2.columns)) // 3) + 1, axis=1).sum().add_prefix('CHO_')
    squared_deltaR_tk_sel3 = deltaR_tk_sel3.transform(lambda x : x**0.5) 

    final_df = pd.concat([squared_deltaR_tk_sel1, squared_deltaR_tk_sel2], axis = 1)

    print(squared_deltaR_tk_sel1)
    print(squared_deltaR_tk_sel2)
    print(final_df)

    final_df.to_csv("squared_deltaR_PSMPOPCup_pruebade5cho35.csv")
    return final_df
