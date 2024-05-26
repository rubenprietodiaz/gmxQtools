import os
import pandas as pd

def read_xvg(file_path):
    """Función para leer un archivo .xvg, ignorando las líneas de comentario."""
    with open(file_path, 'r') as file:
        lines = file.readlines()
    # Filtrar solo las líneas con datos numéricos
    data_lines = [line for line in lines if not line.startswith(('@', '#'))]
    # Convertir a DataFrame
    return pd.DataFrame([list(map(float, line.split())) for line in data_lines], columns=['Time', 'RMSD'])

# Ruta del directorio base donde se encuentran las carpetas 1.md, 2.md, 3.md
base_path = '.'

# Lista de los directorios de estudio
study_dirs = ['1.md', '2.md', '3.md']

# Directorio para guardar los resultados
output_dir = os.path.join(base_path, 'rmsd_results')
os.makedirs(output_dir, exist_ok=True)  # Crear el directorio si no existe

# Buscar subdirectorios comunes en los directorios de estudio
common_subdirs = None
for study_dir in study_dirs:
    subdirs = {subdir for subdir in os.listdir(os.path.join(base_path, study_dir)) if os.path.isdir(os.path.join(base_path, study_dir, subdir))}
    if common_subdirs is None:
        common_subdirs = subdirs
    else:
        common_subdirs &= subdirs

# Procesar cada subdirectorio común
for subdir in common_subdirs:
    # DataFrame para almacenar los datos fusionados
    merged_data = None
    
    # Leer y fusionar los archivos .xvg correspondientes
    for study_dir in study_dirs:
        file_path = os.path.join(base_path, study_dir, subdir, f"{subdir}_rmsd.xvg")
        # Leer los datos
        data = read_xvg(file_path)
        
        if merged_data is None:
            # Si es el primer archivo, tomar la columna 'Time' y la 'RMSD'
            merged_data = data.copy()
            merged_data.columns = ['Time', f'RMSD_{study_dir}']
        else:
            # Solo añadir la columna 'RMSD' renombrada
            merged_data[f'RMSD_{study_dir}'] = data['RMSD']
    
    # Guardar el archivo fusionado en la carpeta de resultados
    output_file_path = os.path.join(output_dir, f"merged_{subdir}_rmsd.csv")
    merged_data.to_csv(output_file_path, index=False)
    print(f"Archivo fusionado guardado en: {output_file_path}")

print("Proceso completado.")

