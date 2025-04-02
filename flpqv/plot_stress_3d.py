import os
import numpy as np
import pprint
import copy
import trimesh
import argparse

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as colors
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

from scipy.spatial.distance import cdist

BohrR = 0.529177249 # Angstrom

class Data_atoms():
    """ class of color list of atoms """
    def __init__(self,
                 name = None,
                 color = None,
                 ):

        self.name =\
        ["E", "H",  "He", "Li", "Be",  "B",  "C",  "N",  "O",  "F",\
        "Ne", "Na", "Mg", "Al", "Si",  "P",  "S", "Cl", "Ar",  "K",\
        "Ca", "Sc", "Ti",  "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu",\
        "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr",  "Y",\
        "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",\
        "Sn", "Sb", "Te",  "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr",\
        "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm",\
        "Yb", "Lu", "Hf", "Ta",  "W", "Re", "Os", "Ir", "Pt", "Au",\
        "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac",\
        "Th", "Pa",  "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es",\
        "Fm", "Md", "No", "Lr"]

        # number = List_element.index("H")                                                                                                                                                                                                                                                                                         
        self.color =\
        ["#FFFFFF", "#FAF0E6", "#F0E68C", "#EE82EE", "#DB7093", "#56256e", "#778899", "#0033ff", "#ee0011", "#7fff7f",\
         "#FF00FF", "#48D1CC", "#8B008B", "#D8BFD8", "#00CED1", "#BC043C", "#FFFF55", "#66CDAA", "#FF69B4", "#ADD8E6",\
         "#5c9291", "#483D8B", "#8B708B", "#008080", "#B22222", "#a033a0", "#BC8F8F", "#50C878", "#66CDAA", "#D2691E",\
         "#BC8F8F", "#F08080", "#E6E6FA", "#DAA520", "#e0ebaf", "#b44c97", "#00a497", "#b77b57", "#006e54", "#d69090",\
         "#634950", "#44617b", "#E9967A", "#FF7F50", "#FFA07A", "#A0522D", "#FFE4C4", "#C0C0C0", "#BDB76B", "#6B8E23",\
         "#556B2F", "#ADFF2F", "#008000", "#7a4171", "#7FFFD4", "#008080", "#00aaaa", "#999900", "#4682B4", "#4169E1",\
         "#7B68EE", "#8A2BE2", "#BA55D3", "#800080", "#FF1493", "#c89932", "#5c9291", "#9d5b8b", "#eeeaec", "#e29676",\
         "#5f6527", "#c70067", "#e9dacb", "#28ff93", "#0582ff", "#baff75", "#43676b", "#47585c", "#fde8d0", "#FFD700",\
         "#dcd6d9", "#bf794e", "#f5b1aa", "#cd5e3c", "#95859c", "#71686c", "#203744", "#ec6d71", "#b55233", "#a19361",\
         "#cc3399", "#3399cc", "#339966", "#ffff00", "#ccff33", "#cc3300", "#cc6600", "#ff0033", "#660066", "#006666",\
         "#ffcc00", "#33cccc", "#99ffff", "#ff6666"]

        self.size =\
        [0.20, 0.41, 0.41, 1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57,
         0.58, 1.66, 1.41, 1.21, 1.11, 1.07, 1.05, 1.02, 1.06, 2.03,
         1.76, 1.70, 1.60, 1.53, 1.39, 1.39, 1.32, 1.26, 1.24, 1.32,
         1.22, 1.22, 1.20, 1.19, 1.20, 1.20, 1.16, 2.20, 1.95, 1.90,
         1.75, 1.64, 1.54, 1.47, 1.46, 1.42, 1.39, 1.45, 1.42, 1.42,
         1.39, 1.39, 1.38, 1.39, 1.40, 2.44, 2.15, 2.07, 2.04, 2.03,
         2.01, 1.99, 1.98, 1.98, 1.96, 1.94, 1.92, 1.92, 1.89, 1.90,
         1.87, 1.87, 1.75, 1.70, 1.62, 1.51, 1.44, 1.41, 1.36, 1.36,
         1.32, 1.45, 1.46, 1.48, 1.40, 1.50, 1.50, 2.60, 2.21, 2.15,
         2.06, 2.00, 1.96, 1.90, 1.87, 1.80, 1.69, 1.60, 1.60, 1.60,
         1.60, 1.60, 1.60, 1.60]

    def hex_to_rgb(hex_code):
        hex_code = hex_code.lstrip("#")  # "#" を削除
        return tuple(int(hex_code[i:i+2], 16) for i in (0, 2, 4))
        # print(hex_to_rgb("#FF5733"))  # (255, 87, 51)

    def __repr__(self):
        return self.__class__.__name__ +"(" + pprint.pformat(self.__dict__) +")"


class Data_cube:
    """ class of data from md file """
    def __init__(self,
                 atomnum = None,
                 numbers = None,   # int
                 valence = None,   # int
                 positions_bohr = None,  # bohr, coordinate of atoms
                 positions_ang = None,  # Ang, coordinate of atoms
                 origin = None, # Ang
                 num_mesh = None,
                 mesh_vec = None, # Ang primitive vector for mesh
                 mesh_grid = None, # Ang
                 val = None,
                 small_mesh_grid = None, # Ang
                 small_val = None):

        pass

    def __repr__(self):
        return self.__class__.__name__ +"(" + pprint.pformat(self.__dict__) +")"

    def read(self,filename):
        with open(filename, 'r') as file:

            # Skip the header lines
            for _ in range(2):
                next(file)

            # Read the number of atoms and the origin
            line = next(file).split()
            self.atomnum = int(line[0])
            #self.origin = [float(coord) for coord in line[1:4]]
            self.origin = [float(coord)*BohrR for coord in line[1:4]]


            # Read the mesh info
            self.num_mesh = np.empty(3).tolist()
            self.mesh_vec = np.empty((3,3)).tolist()

            for i in range(3):
                line = next(file).split()
                self.num_mesh[i] = int(line[0])
                self.mesh_vec[i] = [float(coord) for coord in line[1:4]]

            ## if num_mesh<0 -> Ang, if num_mesh>0 -> Bohr
            ## convert bohr to Angstrom unit
            for i in range(3):
                if self.num_mesh[i]>0:
                    for j in range(3):
                        self.mesh_vec[i][j] = self.mesh_vec[i][j]*BohrR
                else:
                    self.num_mesh[i] = -self.num_mesh[i]


            # Read the atomic symbols and positions_bohr
            self.numbers = []
            self.valence = []
            self.positions_bohr = []
            for _ in range(self.atomnum):
                line = next(file).split()
                self.numbers.append(line[0])
                self.valence.append(line[1])
                self.positions_bohr.append([float(coord) for coord in line[2:5]])

                          ## convert bohr to Angstrom unit
            self.positions_ang = []
            for i in range(self.atomnum):
                self.positions_ang.append([])
                for j in range(3):
                    self.positions_ang[i].append(self.positions_bohr[i][j]*BohrR)

            self.val = np.empty((self.num_mesh[0]*self.num_mesh[1]*self.num_mesh[2])).tolist()

            # Read the file using the next function until the end of the file
            idx = 0
            while True:
                try:
                    line = next(file).strip().split()
                    for value in line:
                        self.val[idx] = float(value)
                        #print(self.val[idx])
                        idx += 1

                except StopIteration:
                    # End of file reached, break out of the loop
                    break

            ### Data structure ###
            #idx = 0
            #for i in range(self.num_mesh[0]):
            #    x = self.origin[0] + i*self.mesh_vec[0]
            #    for j in range(self.num_mesh[1]):
            #        y = self.origin[1] + j*self.mesh_vec[1]
            #        for k in range(num_mesh[2]):
            #            z = self.origin[2] + k*self.mesh_vec[2]
            #            print(x,y,z,self.val[idx])
            #            idx += 1
        return

    def read_header(self,filename):
        with open(filename, 'r') as file:

            # Skip the header lines
            for _ in range(2):
                next(file)

            # Read the number of atoms and the origin
            line = next(file).split()
            self.atomnum = int(line[0])
            self.origin = [float(coord)*BohrR for coord in line[1:4]]

            # Read the mesh info
            self.num_mesh = np.empty(3).tolist()
            self.mesh_vec = np.empty((3,3)).tolist()

            for i in range(3):
                line = next(file).split()
                self.num_mesh[i] = int(line[0])
                self.mesh_vec[i] = [float(coord) for coord in line[1:4]]

            ## if num_mesh<0 -> Ang, if num_mesh>0 -> Bohr
            ## convert bohr to Angstrom unit
            for i in range(3):
                if self.num_mesh[i]>0:
                    for j in range(3):
                        self.mesh_vec[i][j] = self.mesh_vec[i][j]*BohrR
                else:
                    self.num_mesh[i] = -self.num_mesh[i]

            # Read the atomic symbols and positions_bohr
            self.numbers = []
            self.valence = []
            self.positions_bohr = []
            for _ in range(self.atomnum):
                line = next(file).split()
                self.numbers.append(line[0])
                self.valence.append(line[1])
                self.positions_bohr.append([float(coord) for coord in line[2:5]])

            ## convert bohr to Angstrom unit
            self.positions_ang = []
            for i in range(self.atomnum):
                self.positions_ang.append([])
                for j in range(3):
                    self.positions_ang[i].append(self.positions_bohr[i][j]*BohrR)

        return

class Data_tensor():
    """ class of tensor data """
    def __init__(self,
                 filename = None,
                 list_center = None,
                 list_eigenvalues = None,
                 list_eigenvectors = None):

        pass

    def __repr__(self):
        return self.__class__.__name__ +"(" + pprint.pformat(self.__dict__) +")"

class Data_tensor_diag:
    """ class of tensor data """
    def __init__(self,
                 coord = None,
                 val = None,
                 val_dev = None,
                 tensor = None,
                 mesh = None):

        self.comment = []
        self.blank_line =[]

        pass



    def __repr__(self):
        return self.__class__.__name__ +"(" + pprint.pformat(self.__dict__) +")"


    def read(self,filename):
        """
            Read tensor data.
        """
        self.coord  = []
        self.val = []
        self.val_dev = []
        self.tensor = []
        self.mesh = []

        with open(filename.strip(), 'r') as f2:
        
            f = f2.readlines()

            self.mesh.append(np.int64(f[0].strip().split()[2]))
            self.mesh.append(np.int64(f[0].strip().split()[3]))
            self.mesh.append(np.int64(f[0].strip().split()[4]))

            for index, line in enumerate(f):
                if line.startswith("#"):
                    self.comment.append(line)
                    continue

                if line.startswith("\n"):
                    self.blank_line.append(index)
                    continue

                #li = line.strip().split()
                li = list(np.float64(line.strip().split()))

                self.coord.append(np.array([li[0],li[1],li[2]]))
                self.val.append(np.array([li[3],li[4],li[5]]))
                self.tensor.append(np.array([[li[6],li[7],li[8]],
                                 [li[9],li[10],li[11]],
                                 [li[12],li[13],li[14]]]))

                val = np.float64([li[3],li[4],li[5]])
                val_average = (val[0] + val[1] + val[2])/3.0
                self.val_dev.append([val[0]-val_average,val[1]-val_average,val[2]-val_average])

                li.clear()

        return

    def write(self,filename):
        """
            Write tensor data.
        """
        with open(filename.strip(), 'w') as f:
        
            count = 0
            for line in self.comment:
                f.write(line)
                count+=1

            for i in range(len(self.coord)):
                if(count in self.blank_line):
                    f.write("\n")
                else:
                    table = [ self.coord[i] + self.val[i] + self.tensor[i][0] + self.tensor[i][1] + self.tensor[i][2] ]
                    f.write(tabulate(table, tablefmt="plain", floatfmt='.8f'))
                    f.write("\n")
                    table.clear()

                count+=1

def read_toml(toml_file):
    
    #  Python バージョンによる `tomllib` or `toml` の切り替え
    try:
        if sys.version_info >= (3, 11):
            import tomllib
            def load_toml(toml_file):
                with open(toml_file, "rb") as f:
                    return tomllib.load(f)
        else:
            import toml
            def load_toml(toml_file):
                with open(toml_file, "r", encoding="utf-8") as f:
                    return toml.load(f)
    except ImportError:
        print("`toml` module cannot be found\n `pip install toml` ")
        sys.exit(1)
    
    #  TOML ファイルを読み込む
    try:
        #print("start read toml_file")
        data_toml = load_toml(toml_file)
        #print("Succeeded reading TOML file")
        #print(data_toml)
        return data_toml
    except FileNotFoundError:
        print(f" cannot find '{toml_file}' ")
        sys.exit(1)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)

    #  TOML ファイルの読み込み終了



def plot_structure(ax, data_toml, data_cube):
    """
    原子構造をプロットするコード
    """
    # 原子座標を取得
    positions = data_cube.positions_ang
    numbers = np.int64(data_cube.numbers)

    centroid = np.zeros(3)
    if(data_toml["view"]["centering"]):
        # **重心を求める**
        centroid = np.mean(data_cube.positions_ang, axis=0)

    positions = positions - centroid

    # 結合の閾値（Å単位）
    bond_factor = float(data_toml["view"]["bond_factor"]) # 原子間距離がこの値以下なら結合とみなす
    bond_size_ratio = float(data_toml["view"]["bond_size_ratio"])
    atom_size_ratio = float(data_toml["view"]["atom_size_ratio"]) 
    
    # **結合リストを作る**
    distances = cdist(positions, positions)
    bonds = [(i, j) for i in range(len(positions)) for j in range(i+1, len(positions)) if distances[i, j] < bond_factor]

    data_atoms = Data_atoms()
    if not (data_toml["color_list"]["color_atoms"]):
        data_atoms.color[:] = ["#778899"] * len(data_atoms.color)
    
    # 色の設定
    #unique_numbers = list(set(numbers))
    #colors = {el: plt.cm.jet(i / len(unique_numbers)) for i, el in enumerate(unique_numbers)}
    
    # 球を描く関数
    def plot_sphere(ax, center, radius, color):
        u = np.linspace(0, 2 * np.pi, 30)
        v = np.linspace(0, np.pi, 20)
        x = radius * np.outer(np.cos(u), np.sin(v)) + center[0]
        y = radius * np.outer(np.sin(u), np.sin(v)) + center[1]
        z = radius * np.outer(np.ones(np.size(u)), np.cos(v)) + center[2]
        ax.plot_surface(x, y, z, color=color, alpha=1.0, edgecolor="k", linewidth=0.0, antialiased=True)

        vertices = np.array([x.ravel(), y.ravel(), z.ravel()]).T
        faces = []
        rows, cols = x.shape
        for i in range(rows - 1):
            for j in range(cols - 1):
                q1 = i * cols + j
                q2 = q1 + 1
                q3 = (i + 1) * cols + j
                q4 = q3 + 1
                faces.append([q1, q2, q3])  # 三角形1
                faces.append([q2, q4, q3])  # 三角形2

        return vertices, faces


    # **円柱（bond）をプロットする関数**
    def plot_cylinder(ax, p1, p2, radius=0.2, n=20, r_atom1=0.3, r_atom2=0.3):
        """ p1, p2 を結ぶ円柱を描画（原子半径を考慮） """
        v = np.array(p2) - np.array(p1)  # ベクトル p1 → p2
        length = np.linalg.norm(v)  # 円柱の高さ

        # **bond の端点を原子表面にシフト**
        delta_r1 = np.sqrt(float(r_atom1)*float(r_atom1) - float(radius)*float(radius))
        delta_r2 = np.sqrt(float(r_atom2)*float(r_atom2) - float(radius)*float(radius))
        p1_new = p1 + (delta_r1 / length) * v  # p1 から r_atom だけ v 方向へシフト
        p2_new = p2 - (delta_r2 / length) * v  # p2 から r_atom だけ v 方向へシフト
        v_new = p2_new - p1_new
        length_new = np.linalg.norm(v_new)  # 調整後の bond の長さ
    
        # 円柱の基本メッシュ (z 軸方向に伸びる円柱)
        theta = np.linspace(0, 2 * np.pi, n)
        z = np.linspace(0, length_new, 2)  # 円柱の高さ
        theta, z = np.meshgrid(theta, z)
        x = radius * np.cos(theta)
        y = radius * np.sin(theta)
    
        # 方向を合わせる（回転）
        def rotation_matrix(v1, v2):
            """ v1 から v2 に回転する行列を求める """
            v1, v2 = v1 / np.linalg.norm(v1), v2 / np.linalg.norm(v2)
            cross = np.cross(v1, v2)
            dot = np.dot(v1, v2)
            if np.allclose(cross, 0):  # 既に同じ方向なら回転不要
                return np.eye(3)
            skew = np.array([[0, -cross[2], cross[1]], [cross[2], 0, -cross[0]], [-cross[1], cross[0], 0]])
            return np.eye(3) + skew + (skew @ skew) * (1 / (1 + dot))
    
        # z軸基準の円柱を p1 → p2 に回転
        R = rotation_matrix(np.array([0, 0, 1]), v)
        xyz = np.dot(R, np.array([x.ravel(), y.ravel(), z.ravel()]))  # 座標変換
        x, y, z = xyz[0].reshape(x.shape) + p1_new[0], xyz[1].reshape(y.shape) + p1_new[1], xyz[2].reshape(z.shape) + p1_new[2]
    
        # プロット
        ax.plot_surface(x, y, z, color="white", alpha=1.0, linewidth=0.0, antialiased=True)

        vertices = np.array([x.ravel(), y.ravel(), z.ravel()]).T
        faces = []
        rows, cols = x.shape
        for i in range(rows - 1):
            for j in range(cols - 1):
                q1 = i * cols + j
                q2 = q1 + 1
                q3 = (i + 1) * cols + j
                q4 = q3 + 1
                faces.append([q1, q2, q3])  # 三角形1
                faces.append([q2, q4, q3])  # 三角形2

        return vertices, faces


    # **すべての結合を円柱で描画**
    for i, j in bonds:
        plot_cylinder(ax, positions[i], positions[j], radius=0.1*bond_size_ratio,
                                    r_atom1=data_atoms.size[numbers[i]]*atom_size_ratio,
                                    r_atom2=data_atoms.size[numbers[j]]*atom_size_ratio)

    # 各原子を球でプロット
    for pos, num in zip(positions, numbers):
        plot_sphere(ax, pos, data_atoms.size[num]*atom_size_ratio, data_atoms.color[num])

    ## **すべての結合を円柱で描画**
    #list_vertices_faces_bond = []
    #for i, j in bonds:
    #    list_vertices_faces_bond.append(plot_cylinder(ax, positions[i], positions[j], radius=0.1, r_atom=atom_radius))

    #mesh_list_bond = []
    #for tmp in list_vertices_faces_bond:
    #    mesh_list_bond.append(trimesh.Trimesh(vetices=tmp[0], faces=tmp[1]))


    ## 各原子を球でプロット
    #list_vertices_faces_atom= []
    #for pos, elem in zip(positions, numbers):
    #    list_vertices_faces_atom.append(plot_sphere(ax, pos, atom_radius, colors[elem]))

    #mesh_list_atom = []
    #for tmp in list_vertices_faces_atom:
    #    mesh_list_atom.append(trimesh.Trimesh(vetices=tmp[0], faces=tmp[1]))

    return

def plot_plane(ax,data_cube):
    # 1️⃣ 四角形の頂点を定義
    #verts = [
    #    [[0, -1, -1], [0, 1, -1], [0, 1, 1], [0, -1, 1]]  # 四角形の4点
    #]
    #    #[[-1, -1, 0], [1, -1, 0], [1, 1, 0], [-1, 1, 0]]  # 四角形の4点

    #scale = 3
    #verts = scale*verts
    #
    ## 2️⃣ 四角形を描画
    #plane = Poly3DCollection(verts, color='lightblue', edgecolor='k', alpha=0.5)
    #ax.add_collection3d(plane)

    scale = 1
    x = scale*np.linspace(-1, 1, 10)
    y = scale*np.linspace(-1, 1, 10)
    X, Y = np.meshgrid(x, y)
    Z = np.zeros_like(X)  # Z=0 の平面

    ax.plot_surface(X, Y, Z, color='lightblue', alpha=0.5, zorder=20)

def plot_tensor_3d(data_tensor,data_cube,data_toml):

    nskip = 1

    list_coord        = BohrR*np.array(data_tensor.coord).reshape(data_tensor.mesh[0],data_tensor.mesh[1],data_tensor.mesh[2], 3)[::nskip, ::nskip, ::nskip, :].reshape(-1,3)
    list_eigenvalues  = np.array(data_tensor.val).reshape(data_tensor.mesh[0],data_tensor.mesh[1],data_tensor.mesh[2], 3)[::nskip, ::nskip, ::nskip, :].reshape(-1,3)
    list_eigenvectors = np.array(data_tensor.tensor).reshape(data_tensor.mesh[0],data_tensor.mesh[1],data_tensor.mesh[2], 3, 3)[::nskip, ::nskip, ::nskip, :, :].reshape(-1,3,3)

    #print(list_coord)

    #list_coord       = [data_tensor.coord[i] for i in range(0, len(data_tensor.coord), nskip)]
    #list_eigenvalues  = [data_tensor.val[i] for i in range(0, len(data_tensor.val), nskip)]
    #list_eigenvectors = [data_tensor.tensor[i] for i in range(0, len(data_tensor.tensor), nskip)]

    max_val = np.max(list_eigenvalues.flatten())
    min_val = np.min(list_eigenvalues.flatten())

    color_values = np.zeros(len(list_eigenvalues))
    for i, eigenvalues in enumerate(list_eigenvalues):
        color_values[i] = eigenvalues[2]

    mask = np.logical_or(color_values < -0.0001, color_values > 0.0001)
    color_values= color_values[mask]
    list_coord       = list_coord[mask]
    list_eigenvalues  = list_eigenvalues[mask]
    list_eigenvectors = list_eigenvectors[mask]

    # ⭐ 数値データをRGBAカラーに変換する ⭐
    #norm = colors.Normalize(vmin=min_val, vmax=max_val)
    norm = colors.TwoSlopeNorm(vmin=-0.005, vcenter=0, vmax=0.001) 
    #cmap = cm.get_cmap('RdYlBu')  # カラーマップを選択'coolwarm'
    cmap = colors.LinearSegmentedColormap.from_list( "rb_cmap", [(0, "blue"), (0.5, "white"), (1, "red")] )
    list_color = cmap(norm(color_values))  # RGBAに変換！

    list_alpha = np.interp(np.abs(color_values), [0, 0.001], [0.0, 1.0])

    ## 色を最大固有値で決定
    #list_color = []
    #for eigenvalues in list_eigenvalues:
    #    max_eigenvalue = np.max(eigenvalues)
    #    list_color.append( "r" if max_eigenvalue > 0 else "b" ) # 赤（正）or 青（負）

    # 球のメッシュ作成
    u = np.linspace(0, 2 * np.pi, 10)
    v = np.linspace(0, np.pi, 5)
    #x = np.outer(np.cos(u), np.sin(v))
    #y = np.outer(np.sin(u), np.sin(v))
    #z = np.outer(np.ones_like(u), np.cos(v))
    x = (np.outer(np.cos(u), np.sin(v)))**3
    y = (np.outer(np.sin(u), np.sin(v)))**3
    z = (np.outer(np.ones_like(u), np.cos(v)))**3
    
    # 楕円球への変換
    size_point = 1.0
    scale_point = 0.2*size_point
    xyz = scale_point*np.stack([x, y, z], axis=0)

    
    list_translated_xyz = []
    for i, eigenvalues in enumerate(list_eigenvalues):
        scaled_xyz = (0.3+30*np.clip(eigenvalues, -0.005, 0.001))[:, np.newaxis, np.newaxis] * xyz  # スケール適用
        transformed_xyz = np.einsum("ij, jkl -> ikl", list_eigenvectors[i], scaled_xyz)  # 回転適用
        translated_xyz = transformed_xyz + list_coord[i][:, np.newaxis, np.newaxis]  # 中心座標を適用
        list_translated_xyz.append(copy.deepcopy(translated_xyz))
    

    # プロット
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection="3d")
    ax.set_box_aspect([1, 1, 1])  # アスペクト比を統一
    
    for i, translated_xyz in enumerate(list_translated_xyz):
        ax.plot_surface(translated_xyz[0], translated_xyz[1], translated_xyz[2], color=list_color[i], alpha=list_alpha[i], zorder=50)
    
    # 軸の設定
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    
    # 軸範囲を自動調整
    #max_radius = np.max(np.sqrt(eigenvalues))  # 最大半径
    #ax.set_xlim([coord[0] - max_radius, coord[0] + max_radius])
    #ax.set_ylim([coord[1] - max_radius, coord[1] + max_radius])
    #ax.set_zlim([coord[2] - max_radius, coord[2] + max_radius])

    ############  set plot range  ##############
    # 各軸の範囲を取得
    x_min, x_max = ax.get_xlim()
    y_min, y_max = ax.get_ylim()
    z_min, z_max = ax.get_zlim()
    
    # ★ 最大の範囲を求めて、すべての軸に適用
    ratio_zoom = 1.5
    max_range = (max(x_max - x_min, y_max - y_min, z_max - z_min) / 2)/ratio_zoom
    
    # ★ 各軸の中心を求める
    x_mid = (x_max + x_min) / 2
    y_mid = (y_max + y_min) / 2
    z_mid = (z_max + z_min) / 2
    
    # ★ 各軸の範囲を統一する
    ax.set_xlim([x_mid - max_range, x_mid + max_range])
    ax.set_ylim([y_mid - max_range, y_mid + max_range])
    ax.set_zlim([z_mid - max_range, z_mid + max_range])

    ############  set color bar  ##############
    # カラーバーを追加
    sm = cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, shrink=0.5, aspect=20)
    cbar.outline.set_linewidth(0.2)
    #cbar.set_ticks([min_val, max_val])
    cbar.set_ticks([-0.005, 0.0, 0.001])
    #cbar.locator = MaxNLocator(nbins=2)  # 目盛りを5個に設定
    #cbar.update_ticks()  # 目盛りを更新
    # 目盛り線（tick marks）を非表示にする
    cbar.ax.tick_params(axis='both', which='both', length=0)

    #plot_plane(ax,data_cube)
    plot_structure(ax,data_toml,data_cube)

    plt.show()
    #plt.savefig("stress.pdf", bbox_inches='tight', pad_inches=0)

if __name__ == '__main__':

    import sys

    parser = argparse.ArgumentParser(description="Generate BibTeX for cited Python modules.")
    parser.add_argument("filename", nargs="?", help="File to read (optional).")
    parser.add_argument("-cite", action="store_true", help="Print the BibTeX citations.")
    args = parser.parse_args()
    
    filename = args.filename  # コマンドライン引数からファイル名を取得

    data_cube = Data_cube()
    data_cube.read_header("header.cube")
    print(data_cube)
    data_toml = read_toml("input.toml")

    data_tensor = Data_tensor_diag()
    data_tensor.read(filename)
    print("finish reading tensor file")
    plot_tensor_3d(data_tensor,data_cube,data_toml)

    
## 固有値と固有ベクトル（例）
#list_eigenvalues = [np.array([4, 2, 1]), np.array([1, 2, 4]) ] # 軸の長さ（平方根が半径）
#list_eigenvectors = [np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]),  # 正規直交基底
#                     np.array([[1/np.sqrt(2), 0, 0], [0, 1/np.sqrt(2), 0], [0, 0, 1]])]
#
## 楕円球の中心座標を設定
#list_center = [np.array([10.0, -1.0, 2.0]),  # (cx, cy, cz)
#               np.array([0.0, -1.0, 2.0])]  # (cx, cy, cz)



