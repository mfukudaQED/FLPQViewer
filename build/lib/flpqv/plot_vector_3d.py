import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot_vector_field():
    """
    Plot a 3D vector field using matplotlib.
    """
    # Create a grid
    x, y, z = np.meshgrid(np.linspace(-2, 2, 10),
                           np.linspace(-2, 2, 10),
                           np.linspace(-2, 2, 5))
    
    # Define the vector field (example: rotating field)
    u = -y  # x-component
    v = x   # y-component
    w = z   # z-component
    
    # Plot
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.quiver(x, y, z, u, v, w, length=0.3, normalize=True)
    
    # Labels
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('3D Vector Field')
    
    plt.show()

# Execute
plot_vector_field()
