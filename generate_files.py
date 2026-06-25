import numpy as np
import pandas as pd

# Funzione per generare posizioni lungo una circonferenza attorno alle ancore
def generate_circle_positions(anchors, radius, num_points):
    center_x, center_y = anchors['AnchorX'].mean(), anchors['AnchorY'].mean()
    angles = np.linspace(0, 2 * np.pi, num_points, endpoint=False)
    positions = pd.DataFrame({
        'True_T_x': center_x + radius * np.cos(angles),
        'True_T_y': center_y + radius * np.sin(angles),
        'True_T_z': np.zeros(num_points)
    })
    return positions

# Funzione per generare posizioni lungo un lato specifico delle ancore
def generate_side_positions(anchors, side, num_points):
    min_x, max_x = anchors['AnchorX'].min(), anchors['AnchorX'].max()
    min_y, max_y = anchors['AnchorY'].min(), anchors['AnchorY'].max()

    if side == 'left':
        x = np.full(num_points, min_x)
        y = np.linspace(min_y, max_y, num_points)
    elif side == 'right':
        x = np.full(num_points, max_x)
        y = np.linspace(min_y, max_y, num_points)
    elif side == 'bottom':
        x = np.linspace(min_x, max_x, num_points)
        y = np.full(num_points, min_y)
    elif side == 'top':
        x = np.linspace(min_x, max_x, num_points)
        y = np.full(num_points, max_y)
    else:
        raise ValueError("Side must be 'left', 'right', 'bottom', or 'top'")

    positions = pd.DataFrame({'True_T_x': x, 'True_T_y': y, 'True_T_z': np.zeros(num_points)})
    return positions

# Funzione per aggiungere rumore gaussiano
def add_gaussian_noise(positions, std_dev=0.05):
    noisy_positions = positions.copy()
    noisy_positions.iloc[:, 0] += np.random.normal(0, std_dev, size=len(positions))
    noisy_positions.iloc[:, 1] += np.random.normal(0, std_dev, size=len(positions))
    noisy_positions.iloc[:, 2] += np.random.normal(0, std_dev, size=len(positions))
    noisy_positions.columns = ['InputPoint_x', 'InputPoint_y', 'InputPoint_z']
    return noisy_positions

# Funzione per generare posizioni ancore in rettangolo con dimensioni specificate, centrato in 0,0,0
def generate_rectangle_anchors(width, height):
    positions = pd.DataFrame({
        'AnchorX': [-width/2, width/2, width/2, -width/2],
        'AnchorY': [-height/2, -height/2, height/2, height/2],
        'AnchorZ': [0, 0, 0, 0]
    })
    return positions

# Salvataggio CSV
def save_positions_to_csv(positions, file_path):
    positions.to_csv(file_path, index=False)


# Esempio d'uso parametrizzato
if __name__ == "__main__":
    anchors = generate_rectangle_anchors(width=0.4, height=0.6)

    position_type = 'circle'  # 'circle' o 'side'

    if position_type == 'circle':
        true_positions = generate_circle_positions(anchors, radius=3, num_points=50)
    elif position_type == 'side':
        true_positions = generate_side_positions(anchors, side='left', num_points=20)
    else:
        raise ValueError("position_type deve essere 'circle' o 'side'")

    # Salva posizioni vere del tag
    save_positions_to_csv(true_positions, 'TrueTagPositions.csv')

    # Aggiunge rumore e salva come InputPoint
    noisy_positions = add_gaussian_noise(true_positions)
    save_positions_to_csv(noisy_positions, 'InputPoint.csv')

    # Salva posizioni ancore
    save_positions_to_csv(anchors, 'Anchors.csv')
